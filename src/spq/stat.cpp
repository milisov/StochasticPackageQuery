#include <boost/algorithm/string.hpp>
#include <fmt/core.h>
#include <libpq-fe.h>
#include <omp.h>
#include <iostream>

#include "stat.hpp"
#include "util/uconfig.hpp"
#include "util/udebug.hpp"
#include "util/udeclare.hpp"
#include "core/kde.hpp"

using boost::algorithm::to_lower_copy;
using std::make_unique;
using std::cerr;
using std::to_string;

bool Stat::rebuildStat = to_lower_copy(Config::getInstance()->pt.get<string>("build.rebuild_stat")) == "true";
const string Stat::statTable = Config::getInstance()->pt.get<string>("pgmanager.stat_table");

Stat::Stat(){
    pg = make_unique<PgManager>();
    if (!pg->existTable(statTable) || rebuildStat){
        rebuildStat = false;
        if (pg->existTable(statTable)) pg->dropTable(statTable);
        string sql = fmt::format("CREATE TABLE \"{}\" (\
            table_name VARCHAR(50), \
            column_name VARCHAR(50), \
            sum NUMERIC, \
            m2 NUMERIC, \
            count BIGINT)", statTable);
        auto res = PQexec(pg->conn.get(), sql.c_str());
        ck(pg->conn, res);
        PQclear(res);
    }
}

bool Stat::isAnalyzed(const string& tableName, const string& columnName){
    string sql = fmt::format("SELECT count FROM \"{}\" WHERE table_name='{}' AND column_name='{}'", statTable, tableName, columnName);
    auto res = PQexec(pg->conn.get(), sql.c_str());
    ck(pg->conn, res);
    bool analyzed = PQntuples(res) > 0;
    PQclear(res);
    return analyzed;
}

bool Stat::analyze(const string& tableName){
    if (!pg->existTable(tableName)){
        cerr << fmt::format("Table '{}' does not exist\n", tableName);
        return false;
    }
    auto strColumns = Config::getInstance()->pt.get_optional<string>("partition."+tableName);
    if (strColumns){
        vector<string> columnNames;
        boost::split(columnNames, *strColumns, boost::is_any_of(","));
        auto columns = pg->getColumns(tableName);
        for (auto columnName : columnNames){
            if (!columns.count(columnName)){
                cerr << fmt::format("Column '{}' does not exist in table '{}'\n", columnName, tableName);
                return false;
            }
            if (columns[columnName] != Column::array_type && columns[columnName] != Column::numeric_type){
                cerr << fmt::format("Column '{}' is not supported for partitioning in table '{}'\n", columnName, tableName);
                return false;
            }
        }

        vector<string> dColumnNames;
        for (auto columnName : columnNames){
            if (!isAnalyzed(tableName, columnName)){
                if (columns[columnName] == Column::numeric_type){
                    dColumnNames.push_back(columnName);
                } else if (columns[columnName] == Column::array_type){
                    analyzeStochastic(tableName, columnName);
                }
            }
        }
        analyzeDeterministic(tableName, dColumnNames);
        return true;
    } else{
        cerr << fmt::format("No partitioning columns defined for table '{}'\n", tableName);
        return false;
    }
}

void Stat::analyzeStochastic(const string& tableName, const string& columnName){
    long long tableSize = pg->getTableSize(tableName);
    auto columns = pg->getColumns(tableName);
    if (columns[columnName] == Column::array_type){
        string summaryTable = fmt::format("{}_summary", tableName);
        if (!pg->existTable(summaryTable)){
            string sql = fmt::format("CREATE TABLE \"{}\" ({} SERIAL PRIMARY KEY)", summaryTable, PgManager::id);
            auto res = PQexec(pg->conn.get(), sql.c_str());
            ck(pg->conn, res);
            PQclear(res);
            sql = fmt::format("INSERT INTO \"{}\" ({}) SELECT GENERATE_SERIES(1, {})", summaryTable, PgManager::id, tableSize);
            res = PQexec(pg->conn.get(), sql.c_str());
            ck(pg->conn, res);
            PQclear(res);
        }
        auto summaryColumns = pg->getColumns(summaryTable);
        vector<string> summaryColumnNames = {columnName + "_mean", columnName + "_variance", columnName + "_qseq"};
        vector<string> summaryColumnTypes = {"DOUBLE PRECISION", "DOUBLE PRECISION", "REAL[]"};
        vector<string> updateColumnNames (summaryColumnNames.size());
        for (size_t i = 0; i < summaryColumnNames.size(); ++i){
            if (!summaryColumns.count(summaryColumnNames[i])){
                pg->addColumn(summaryTable, summaryColumnNames[i], summaryColumnTypes[i]);
            }
            updateColumnNames[i] = fmt::format("{}=${}", summaryColumnNames[i], i+2);
        }
        auto nCores = Config::getInstance()->nPhysicalCores;
        auto intervals = divideInterval(1, tableSize, nCores);
        int count = 0, N = pg->getColumnLength(tableName, columnName);
        vector<float> quantiles (2*N+1, 0);
        quantiles.back() = 1;
        for (size_t i=1; i < quantiles.size()-1; ++i){
            quantiles[i] = static_cast<float>(i)/(2*N);
        }
        AccAggregator agg;
        #pragma omp parallel num_threads(nCores)
        {
            auto coreIndex = omp_get_thread_num();
            AsyncUpdate au = AsyncUpdate(fmt::format("UPDATE \"{}\" SET {} WHERE {}=$1", summaryTable, boost::join(updateColumnNames, ","), PgManager::id));
            SingleRow sr = SingleRow(fmt::format("SELECT {},{} FROM \"{}\" WHERE {} BETWEEN {} AND {}", PgManager::id, columnName, tableName, PgManager::id, intervals[coreIndex], intervals[coreIndex+1]-1));
            while (sr.fetchRow()){
                long long id = sr.getBigInt(0);
                vector<double> samples;
                sr.getDoubleArray(1, samples);
                KDE kde (samples, true);
                au.set(0, id);
                au.set(1, kde.getMean());
                au.set(2, kde.getVariance());
                vector<float> qseq = quantiles;
                kde.getQuantiles(qseq);
                au.set(3, qseq);
                au.send();
                #pragma omp critical
                {
                    agg.add(kde.getSum(), kde.getM2(), samples.size());
                }
            }
        }
        addStat(tableName, columnName, agg.getSum(), agg.getM2(), agg.getCount());
    }
}

void Stat::analyzeDeterministic(const string& tableName, const vector<string>& columnNames){
    auto columns = pg->getColumns(tableName);
    auto nCores = Config::getInstance()->nPhysicalCores;
    auto intervals = divideInterval(1, pg->getTableSize(tableName), nCores);
    string selectColumns = boost::join(columnNames, ",");
    vector<AccAggregator> aggs (columnNames.size());
    #pragma omp parallel num_threads(nCores)
    {
        auto coreIndex = omp_get_thread_num();
        string sql = fmt::format("SELECT {} FROM \"{}\" WHERE {} BETWEEN {} AND {}", selectColumns, tableName, PgManager::id, intervals[coreIndex], intervals[coreIndex+1]-1);
        SingleRow sr = SingleRow(sql);
        vector<AccSet> accs (columnNames.size());
        long long count = 0;
        while (sr.fetchRow()){
            for (size_t j = 0; j < columnNames.size(); ++j){
                if (columns[columnNames[j]] == Column::numeric_type){
                    accs[j](sr.getNumeric(j));
                }
            }
            count ++;
        }
        #pragma omp critical
        {
            for (size_t j = 0; j < columnNames.size(); ++j){
                aggs[j].add(ba::sum(accs[j]), count*ba::variance(accs[j]), count);
            }
        }
    }
    for (size_t j = 0; j < columnNames.size(); ++j){
        if (columns[columnNames[j]] == Column::numeric_type){
            addStat(tableName, columnNames[j], aggs[j].getSum(), aggs[j].getM2(), aggs[j].getCount());
        }
    }
}

void Stat::addStat(const string& tableName, const string& columnName, const long_double sum, const long_double m2, const long long count){
    string sql;
    if (!isAnalyzed(tableName, columnName)){
        sql = fmt::format("INSERT INTO \"{}\" VALUES ('{}', '{}', {}, {}, {})", statTable, tableName, columnName, strf(sum), strf(m2), count);
    } else{
        long_double avg = sum / count;
        sql = fmt::format("UPDATE \"{}\" SET \
            m2=m2+{}+(sum/count+{})*(sum/count+{})*count/(count+{})*{}\
            sum=sum+{}, count=count+{} \
            WHERE table_name='{}' AND column_name='{}'",
            statTable, strf(m2), strf(avg), strf(avg), count, count, strf(sum), count, tableName, columnName);
    }
    auto res = PQexec(pg->conn.get(), sql.c_str());
    ck(pg->conn, res);
    PQclear(res);
}

void Stat::getStochasticMeanVar(const string& tableName, const string& columnName, unordered_map<long long, pair<double, double>>& storages){
    if (!isAnalyzed(tableName, columnName)){
        cerr << fmt::format("Column '{}' has not been analyzed in table '{}'\n", columnName, tableName);
        exit(1);
    }
    vector<string> strIds; strIds.reserve(storages.size());
    for (const auto& p : storages) {
        strIds.push_back(to_string(p.first));
    }
    string sql = fmt::format("SELECT {},{}_mean,{}_variance FROM \"{}_summary\" WHERE {} IN ({})", PgManager::id, columnName, columnName, tableName, PgManager::id, boost::join(strIds, ","));
    auto res = PQexec(pg->conn.get(), sql.c_str());
    ck(pg->conn, res);
    for (int i = 0; i < PQntuples(res); ++i){
        storages[atoll(PQgetvalue(res, i, 0))] = {atof(PQgetvalue(res, i, 1)), atof(PQgetvalue(res, i, 2))};
    }
    PQclear(res);
}