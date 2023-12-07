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

        vector<string> dColumnNames, sColumnNames;
        bool rebuildSummary = false;
        for (auto columnName : columnNames){
            if (!isAnalyzed(tableName, columnName)){
                if (columns[columnName] == Column::numeric_type){
                    dColumnNames.push_back(columnName);
                } else if (columns[columnName] == Column::array_type){
                    rebuildSummary = true;
                }
            }
            if (columns[columnName] == Column::array_type){
                sColumnNames.push_back(columnName);
            }
        }
        analyzeDeterministic(tableName, dColumnNames);
        if (rebuildSummary) analyzeStochastic(tableName, sColumnNames);
        return true;
    } else{
        cerr << fmt::format("No partitioning columns defined for table '{}'\n", tableName);
        return false;
    }
}

void Stat::analyzeStochastic(const string& tableName, const vector<string>& columnNames){
    long long tableSize = pg->getTableSize(tableName);
    string summaryTable = fmt::format("{}_summary", tableName);
    if (pg->existTable(summaryTable)) pg->dropTable(summaryTable);
    vector<string> summaryColumnPostfixes = {"_mean", "_variance", "_quantiles"};
    vector<string> summaryColumnTypes = {"DOUBLE PRECISION", "DOUBLE PRECISION", "REAL[]"};
    auto perCol = summaryColumnPostfixes.size();
    auto nColumns = columnNames.size();
    vector<string> summaryColumns (perCol*nColumns);
    vector<string> createColumns (perCol*nColumns);
    for (size_t i = 0; i < nColumns; ++i){
        for (size_t j = 0; j < perCol; ++j){
            auto ind = i*perCol+j;
            summaryColumns[ind] = columnNames[i] + summaryColumnPostfixes[j];
            createColumns[ind] = fmt::format("{} {}", summaryColumns[ind], summaryColumnTypes[j]);
        }
    }
    string sql = fmt::format("CREATE TABLE \"{}\" ({} BIGINT,{})", summaryTable, PgManager::id, boost::join(createColumns, ","));
    auto res = PQexec(pg->conn.get(), sql.c_str());
    ck(pg->conn, res);
    PQclear(res);
    vector<vector<float>> quantiles (nColumns);
    for (size_t i = 0; i < nColumns; ++i){
        int N = pg->getColumnLength(tableName, columnNames[i]);
        vector<float> quantile (2*N+1, 0);
        quantile.back() = 1;
        for (size_t i=1; i < quantile.size()-1; ++i){
            quantile[i] = static_cast<float>(i)/(2*N);
        }
        quantiles[i] = quantile;
    }
    vector<AccAggregator> aggs (nColumns);
    vector<omp_lock_t> locks (nColumns);
    for (size_t i = 0; i < nColumns; ++i) omp_init_lock(&(locks[i]));
    auto nCores = Config::getInstance()->nPhysicalCores;
    auto intervals = divideInterval(1, tableSize, nCores);
    string selectColumns = boost::join(columnNames, ",");
    #pragma omp parallel num_threads(nCores)
    {
        auto coreIndex = omp_get_thread_num();
        SingleRow sr = SingleRow(fmt::format("SELECT {},{} FROM \"{}\" WHERE {} BETWEEN {} AND {}", PgManager::id, selectColumns, tableName, PgManager::id, intervals[coreIndex], intervals[coreIndex+1]-1));
        BulkCopy bc = BulkCopy(summaryTable);
        while (sr.fetchRow()){
            long long id = sr.getBigInt(0);
            bc.add(id);
            for (size_t i = 0; i < nColumns; ++i){
                vector<double> samples;
                sr.getArray(i+1, samples);
                KDE kde (samples, true);
                vector<float> qseq = quantiles[i];
                kde.getSortedQuantiles(qseq);
                bc.add(kde.getMean());
                bc.add(kde.getVariance());
                bc.add(qseq);
                omp_set_lock(&(locks[i]));
                aggs[i].add(kde.getSum(), kde.getM2(), samples.size());
                omp_unset_lock(&(locks[i]));
            }
            bc.send();
        }
    }
    for (size_t i = 0; i < nColumns; ++i){
        omp_destroy_lock(&(locks[i]));
        const auto& agg = aggs[i];
        addStat(tableName, columnNames[i], agg.getSum(), agg.getM2(), agg.getCount());
    }
    sql = fmt::format("CREATE INDEX \"{}_{}\" ON \"{}\" ({})", PgManager::id, summaryTable, summaryTable, PgManager::id);
    res = PQexec(pg->conn.get(), sql.c_str());
    ck(pg->conn, res);
    PQclear(res);
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

void Stat::addStat(const string& tableName, const string& columnName, const long_double& sum, const long_double& m2, const long long& count){
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

void Stat::getStoMeanVars(const string& tableName, const string& columnName, const string& sqlId, vector<double>& means, vector<double>& vars){
    if (!isAnalyzed(tableName, columnName)){
        cerr << fmt::format("Column '{}' has not been analyzed in table '{}'\n", columnName, tableName);
        exit(1);
    }
    string sql = fmt::format("SELECT {}_mean,{}_variance FROM \"{}_summary\" WHERE {}", columnName, columnName, tableName, sqlId);
    auto res = PQexec(pg->conn.get(), sql.c_str());
    ck(pg->conn, res);
    for (int i = 0; i < PQntuples(res); ++i){
        means.emplace_back(atof(PQgetvalue(res, i, 0)));
        vars.emplace_back(atof(PQgetvalue(res, i, 1)));
    }
    PQclear(res);
}

void Stat::getDetAttrs(const string& tableName, const string& columnName, const string& sqlId, vector<double>& attrs){
    string sql, column, table;
    if (pg->getColumns(tableName)[columnName] == Column::array_type){
        column = columnName + "_mean";
        table = tableName + "_summary";
    } else{
        column = columnName;
        table = tableName;
    }
    sql = fmt::format("SELECT {} FROM \"{}\" WHERE {}", column, table, sqlId);
    auto res = PQexec(pg->conn.get(), sql.c_str());
    ck(pg->conn, res);
    for (int i = 0; i < PQntuples(res); ++i){
        attrs.emplace_back(atof(PQgetvalue(res, i, 0)));
    }
    PQclear(res);
}

pair<double, double> Stat::getDetMeanVar(const string& tableName, const string& columnName){
    if (!isAnalyzed(tableName, columnName)){
        cerr << fmt::format("Column '{}' has not been analyzed in table '{}'\n", columnName, tableName);
        exit(1);
    }
    string sql = fmt::format("SELECT sum/count,m2/count FROM \"{}\" WHERE table_name='{}' AND column_name='{}'", statTable, tableName, columnName);
    auto res = PQexec(pg->conn.get(), sql.c_str());
    ck(pg->conn, res);
    pair<double, double> result;
    result.first = atof(PQgetvalue(res, 0, 0));
    result.second = atof(PQgetvalue(res, 0, 1));
    PQclear(res);
    return result;
}

void Stat::getSamples(const string& tableName, const string& columnName, const long long& id, vector<double>& samples){
    string sql = fmt::format("SELECT {} FROM \"{}\" WHERE {}={}", columnName, tableName, PgManager::id, id);
    auto res = PQexec(pg->conn.get(), sql.c_str());
    ck(pg->conn, res);
    readArray(PQgetvalue(res, 0, 0), samples);
    PQclear(res);
}

void Stat::getSamples(const string& tableName, const string& columnName, const string& sqlId, const vector<size_t>& sampleIds, vector<vector<double>>& samples){
    vector<string> strSampleIds (sampleIds.size());
    for (size_t i = 0; i < sampleIds.size(); ++i) strSampleIds[i] = fmt::format("{}[{}]", columnName, sampleIds[i]);
    string sql = fmt::format("SELECT {} FROM \"{}\" WHERE {}", boost::join(strSampleIds, ","), tableName, sqlId);
    // auto res = PQexec(pg->conn.get(), sql.c_str());
    // ck(pg->conn, res);
    // readArray(PQgetvalue(res, 0, 0), samples);
    // PQclear(res);
}

void Stat::getQuantiles(const string& tableName, const string& columnName, const long long& id, vector<double>& quantiles){
    string sql = fmt::format("SELECT {}_quantiles FROM \"{}_summary\" WHERE {}={}", columnName, tableName, PgManager::id, id);
    auto res = PQexec(pg->conn.get(), sql.c_str());
    ck(pg->conn, res);
    readArray(PQgetvalue(res, 0, 0), quantiles);
    PQclear(res);
}