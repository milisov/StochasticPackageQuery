#include <boost/algorithm/string.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <fmt/core.h>
#include <libpq-fe.h>
#include <omp.h>
#include <iostream>

#include "stat.hpp"
#include "util/uconfig.hpp"
#include "util/udebug.hpp"
#include "util/unumeric.hpp"

using boost::algorithm::to_lower_copy;
using std::make_unique;
using std::cerr;

namespace ba = boost::accumulators;
typedef ba::accumulator_set<long double, ba::stats<ba::tag::sum, ba::tag::variance>> AccSet;

bool Stat::rebuildStat = to_lower_copy(Config::getInstance()->pt.get<string>("build.rebuild_stat")) == "true";
string Stat::statTable = Config::getInstance()->pt.get<string>("database.stat_table");

Stat::Stat(){
    pg = make_unique<PgManager>();
    if (!pg->existTable(statTable) || rebuildStat){
        rebuildStat = false;
        if (pg->existTable(statTable)) pg->dropTable(statTable);
        string sql = fmt::format("CREATE TABLE {} (\
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
    string sql = fmt::format("SELECT * FROM \"{}\" WHERE table_name='{}' AND column_name='{}'", statTable, tableName, columnName);
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
        for (size_t i = 0; i < summaryColumnNames.size(); ++i){
            if (!summaryColumns.count(summaryColumnNames[i])){
                pg->addColumn(summaryTable, summaryColumnNames[i], summaryColumnTypes[i]);
            }
        }
    }
}

void Stat::analyzeDeterministic(const string& tableName, const vector<string>& columnNames){
    auto columns = pg->getColumns(tableName);
    auto nPhysicalCores = Config::getInstance()->nPhysicalCores;
    auto intervals = divideInterval(1, pg->getTableSize(tableName), nPhysicalCores);
    string selectColumns = boost::join(columnNames, ",");
    vector<AccAggregator> aggs (columnNames.size());
    #pragma omp parallel num_threads(nPhysicalCores)
    {
        auto coreIndex = omp_get_thread_num();
        PgManager pg_;
        string sql = fmt::format("SELECT {} FROM \"{}\" WHERE {} BETWEEN {} AND {}", selectColumns, tableName, PgManager::id, intervals[coreIndex], intervals[coreIndex+1]-1);
        auto res = PQexec(pg_.conn.get(), sql.c_str());
        ck(pg_.conn, res);
        vector<AccSet> accs (columnNames.size());
        long long count = PQntuples(res);
        for (long long i = 0; i < count; ++i){
            for (size_t j = 0; j < columnNames.size(); ++j){
                if (columns[columnNames[j]] == Column::numeric_type){
                    accs[j](atof(PQgetvalue(res, i, j)));
                }
            }
        }
        PQclear(res);
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