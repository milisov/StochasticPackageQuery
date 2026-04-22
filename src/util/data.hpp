#include <iostream>
#include <unordered_map>
#include "util/udebug.hpp"
#include "spq/spq.hpp"
#include "core/stat.hpp"
#include "spq/bounder.hpp"
#include "core/checker.hpp"
#include "solver/starter.hpp"
#include "util/uio.hpp"
#include "util/udeclare.hpp"
#include "spq/cons.hpp"
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <set>
#include "gurobi_c++.h"
#include <cstring>
#include <cstdint>
#include <vector>
#pragma once

// Byte-swapping for Big-Endian to Little-Endian conversion
inline int16_t bswap16(int16_t v) { return __builtin_bswap16(v); }
inline int32_t bswap32(int32_t v) { return __builtin_bswap32(v); }
inline int64_t bswap64(int64_t v) { return __builtin_bswap64(v); }

inline int16_t readInt16(const char *&ptr)
{
    int16_t val;
    std::memcpy(&val, ptr, 2);
    ptr += 2;
    return bswap16(val);
}

inline int32_t readInt32(const char *&ptr)
{
    int32_t val;
    std::memcpy(&val, ptr, 4);
    ptr += 4;
    return bswap32(val);
}

inline int64_t readInt64(const char *&ptr)
{
    int64_t val;
    std::memcpy(&val, ptr, 8);
    ptr += 8;
    return bswap64(val);
}

inline float readFloat(const char *&ptr)
{
    uint32_t val;
    std::memcpy(&val, ptr, 4);
    ptr += 4;
    val = bswap32(val);
    float result;
    std::memcpy(&result, &val, 4);
    return result;
}

inline double readDouble(const char *&ptr)
{
    uint64_t val;
    std::memcpy(&val, ptr, 8);
    ptr += 8;
    val = bswap64(val);
    double result;
    std::memcpy(&result, &val, 8);
    return result;
}

// Universal Postgres Binary Array Parser
inline void parsePgBinaryArray(const char *&ptr, std::vector<double> &out, bool isFloat4)
{
    int32_t ndim = readInt32(ptr);
    int32_t hasnull = readInt32(ptr);
    int32_t elem_oid = readInt32(ptr);

    if (ndim == 0)
        return;

    int32_t dim = readInt32(ptr);
    int32_t lb = readInt32(ptr); // Skip lower bound

    out.reserve(dim); // Pre-allocate to prevent vector reallocation

    for (int i = 0; i < dim; ++i)
    {
        int32_t itemlen = readInt32(ptr);
        if (itemlen == -1)
        {
            out.push_back(0.0);
        }
        else
        {
            if (isFloat4)
                out.push_back(static_cast<double>(readFloat(ptr)));
            else
                out.push_back(readDouble(ptr));
        }
    }
}

// Routing strategies for binary parsing
enum class ParseType
{
    INT16,
    INT32,
    INT64,
    FLOAT32,
    FLOAT64,
    FLOAT32_ARRAY,
    FLOAT64_ARRAY,
    TEXT,
    IGNORE
};

inline ParseType getParseStrategy(const std::string& pgType) {
    // Standard scalar types
    if (pgType == "int2" || pgType == "smallint") return ParseType::INT16;
    if (pgType == "int4" || pgType == "integer") return ParseType::INT32;
    if (pgType == "int8" || pgType == "bigint") return ParseType::INT64;
    if (pgType == "float4" || pgType == "real") return ParseType::FLOAT32;
    if (pgType == "float8" || pgType == "double precision") return ParseType::FLOAT64;
    
    // Array types (udt_name always prefixes arrays with an underscore)
    if (pgType == "_float4") return ParseType::FLOAT32_ARRAY;
    if (pgType == "_float8") return ParseType::FLOAT64_ARRAY;
    
    return ParseType::IGNORE;
}

class Data
{
public:
    // Singleton access
    static Data &getInstance()
    {
        if (!instance)
        {
            throw std::runtime_error("Data::init() must be called before getInstance()");
        }
        return *instance;
    }

    static void init(std::shared_ptr<StochasticPackageQuery> spq)
    {
        if (!instance)
        {
            instance = std::unique_ptr<Data>(new Data(spq));
        }
    }

    static void reset()
    {
        instance.reset(); // Sets instance to nullptr
    }

    // Deterministic attributes: attrName -> vector<double>
    std::unordered_map<std::string, std::vector<double>> detAttrs;

    // Stochastic attributes: attrName -> vector<vector<double>>
    std::unordered_map<std::string, std::vector<std::vector<double>>> stochAttrs;

    // For expected profit values
    std::vector<double> stockExpectedProfit;
    std::vector<std::pair<int, double>> stockExpectedProfitSorted;

    // Size of tuples in the table
    std::shared_ptr<StochasticPackageQuery> spq;
    PgManager pg;
    int NTuples;
    std::string DB_optim;
    std::string DB_valid;
    int cntScenarios;
    std::set<std::string> columns;
    std::set<std::string> statColumns;

    // get the columns that are part of the query
    void fetchColumns(std::shared_ptr<StochasticPackageQuery> spq)
    {
        for (const auto &con : spq->cons)
        {
            std::shared_ptr<ProbConstraint> probCon;
            std::shared_ptr<AttrConstraint> attrCon;
            std::shared_ptr<BoundConstraint> boundCon;
            if (isStochastic(con, probCon, attrCon))
            {
                columns.insert(attrCon->attr);
            }
            else
            {
                boundCon = std::dynamic_pointer_cast<BoundConstraint>(con);
                attrCon = std::dynamic_pointer_cast<AttrConstraint>(con);
                if (boundCon && attrCon)
                {
                    columns.insert(attrCon->attr);
                }
            }
            std::shared_ptr<AttrObjective> attrObj;
            bool isDet = isDeterministic(spq->obj, attrObj);
            if (isDet && attrObj->objType == numeric_type)
            {
                columns.insert(attrObj->obj);
            }

            std::shared_ptr<AttrObjective> attrObj2;
            bool isDet2 = isDeterministic(spq->obj, attrObj2);
            if (isDet2 && attrObj2->objType == array_type)
            {
                columns.insert(attrObj2->obj);
                statColumns.insert(attrObj2->obj);
            }
            // initStockExpectedProfit(fmt::format("{}_{}", DB_optim, "summary"), "profit_quantiles");
        }
    }

    double computeAverage(std::vector<double> &values)
    {
        double sum = 0.0;
        for (int i = 0; i < values.size(); i++)
        {
            sum += values[i];
        }
        return sum / values.size();
    }

    // void fetchData()
    // {
    //     fetchColumns(spq);
    //     auto dataColumns = pg.getColumns(spq->tableName);
    //     std::vector<std::string> cols(columns.begin(), columns.end());
    //     std::string selectcols = fmt::format("{}", fmt::join(cols, ", "));

    //     std::string sql = fmt::format("SELECT {} FROM \"{}\" ORDER BY id", selectcols, spq->tableName);

    //     SingleRow sr(sql);
    //     int id = 0;
    //     while (sr.fetchRow())
    //     {
    //         for (int i = 0; i < cols.size(); i++)
    //         {
    //             std::string colName = cols[i];
    //             Column colType = dataColumns[colName];
    //             if (colType == numeric_type)
    //             {
    //                 detAttrs[colName].push_back(sr.getNumeric(i));
    //             }
    //             else
    //             {
    //                 std::vector<double> values;
    //                 sr.getArray(i, values);
    //                 stochAttrs[colName].push_back(values);
    //                 if (statColumns.count(colName) > 0)
    //                 {
    //                     double avg = computeAverage(values);
    //                     stockExpectedProfit.push_back(avg);
    //                     stockExpectedProfitSorted.emplace_back(id, avg);
    //                 }
    //             }
    //         }
    //         id++;
    //     }
    //     std::sort(stockExpectedProfitSorted.begin(), stockExpectedProfitSorted.end(),
    //               [](const auto &a, const auto &b)
    //               { return a.second > b.second; });

    //     //deb(stockExpectedProfitSorted, stockExpectedProfitSorted.size());
    // }

    void fetchData()
    {
        fetchColumns(spq);

        auto exactPgTypes = pg.getExactPgTypes(spq->tableName);
        auto dataColumns = pg.getColumns(spq->tableName);

        std::vector<std::string> cols(columns.begin(), columns.end());
        std::string selectcols = fmt::format("{}", fmt::join(cols, ", "));
        std::string sql = fmt::format("SELECT {} FROM \"{}\" ORDER BY id", selectcols, spq->tableName);

        // Pre-allocate memory to stop RAM thrashing
        if (NTuples > 0)
        {
            stockExpectedProfit.reserve(NTuples);
            stockExpectedProfitSorted.reserve(NTuples);
        }

        std::vector<std::vector<double> *> detDests(cols.size(), nullptr);
        std::vector<std::vector<std::vector<double>> *> stochDests(cols.size(), nullptr);
        std::vector<ParseType> parseStrategies(cols.size(), ParseType::IGNORE);
        std::vector<bool> isStat(cols.size(), false);

        // Map the parsing strategy per column
        for (size_t i = 0; i < cols.size(); i++)
        {
            const std::string &colName = cols[i];
            std::string exactType = exactPgTypes[colName];

            parseStrategies[i] = getParseStrategy(exactType);

            if (dataColumns[colName] == numeric_type)
            {
                detAttrs[colName].reserve(NTuples);
                detDests[i] = &detAttrs[colName];
            }
            else if (dataColumns[colName] == array_type)
            {
                stochAttrs[colName].reserve(NTuples);
                stochDests[i] = &stochAttrs[colName];
                if (statColumns.count(colName) > 0)
                    isStat[i] = true;
            }
        }

        cout<<"mapped columns, starting to fetch data"<<endl;

        // Pull the raw binary payload
        BulkFetch fetcher(pg.conn.get(), sql);
        std::vector<char> rawData = fetcher.fetchAll();
        cout<<"fetched"<<endl;
        if (rawData.empty())
            return;

        const char *ptr = rawData.data();
        const char *end = rawData.data() + rawData.size();

        // Verify Postgres Binary COPY Signature
        const char pg_sig[] = "PGCOPY\n\377\r\n\0";
        if (std::memcmp(ptr, pg_sig, 11) != 0)
        {
            throw std::runtime_error("Invalid binary COPY signature");
        }
        ptr += 11;

        // Skip flags and extension block
        ptr += 4;
        int32_t ext_len = readInt32(ptr);
        ptr += ext_len;

        // Parse the binary tuples
        int id = 0;
        while (ptr < end)
        {
            int16_t field_count = readInt16(ptr);
            if (field_count == -1)
                break; // EOF Marker

            for (int i = 0; i < field_count; i++)
            {
                int32_t field_len = readInt32(ptr);
                if (field_len == -1)
                    continue; // SQL NULL marker

                switch (parseStrategies[i])
                {
                case ParseType::INT16:
                    if (detDests[i])
                        detDests[i]->push_back(static_cast<double>(readInt16(ptr)));
                    else ptr += field_len;
                    break;
                case ParseType::INT32:
                    if (detDests[i])
                        detDests[i]->push_back(static_cast<double>(readInt32(ptr)));
                    else ptr += field_len;
                    break;
                case ParseType::INT64:
                    if (detDests[i])
                        detDests[i]->push_back(static_cast<double>(readInt64(ptr)));
                    else ptr += field_len;
                    break;
                case ParseType::FLOAT32:
                    if (detDests[i])
                        detDests[i]->push_back(static_cast<double>(readFloat(ptr)));
                    else ptr += field_len;
                    break;
                case ParseType::FLOAT64:
                    if (detDests[i])
                        detDests[i]->push_back(readDouble(ptr));
                    else ptr += field_len;
                    break;
                case ParseType::FLOAT32_ARRAY:
                case ParseType::FLOAT64_ARRAY:
                {
                    if (stochDests[i])
                    {
                        std::vector<double> values;
                        bool isFloat4 = (parseStrategies[i] == ParseType::FLOAT32_ARRAY);
                        parsePgBinaryArray(ptr, values, isFloat4);
                        stochDests[i]->push_back(std::move(values));

                        if (isStat[i])
                        {
                            double avg = computeAverage(stochDests[i]->back());
                            stockExpectedProfit.push_back(avg);
                            stockExpectedProfitSorted.emplace_back(id, avg);
                        }
                    }
                    else
                    {
                        ptr += field_len;
                    }
                    break;
                }
                case ParseType::TEXT:
                case ParseType::IGNORE:
                default:
                    ptr += field_len; // Advance pointer safely
                    break;
                }
            }
            id++;
        }

        std::sort(stockExpectedProfitSorted.begin(), stockExpectedProfitSorted.end(),
                  [](const auto &a, const auto &b)
                  { return a.second > b.second; });

        const auto& lastRowScenarios = stochAttrs["profit"].back();
        deb(stochAttrs["profit"].size());
        deb(lastRowScenarios);
    }

    // Compute expected profit values (needed for ExpSum constraints/obj)
    void initStockExpectedProfit(const std::string &table, const std::string &colName)
    {
        std::string sql = fmt::format(
            "SELECT id, AVG(val) AS avg_val "
            "FROM \"{}\", unnest({}) AS val "
            "GROUP BY id",
            table, colName);

        SingleRow sr(sql);

        while (sr.fetchRow())
        {
            double mean = sr.getNumeric(1);
            int id = sr.getBigInt(0) - 1;
            stockExpectedProfit[id] = mean;
        }

        stockExpectedProfitSorted.clear();
        for (int i = 0; i < (int)stockExpectedProfit.size(); i++)
        {
            stockExpectedProfitSorted.emplace_back(i, stockExpectedProfit[i]);
        }

        std::sort(stockExpectedProfitSorted.begin(), stockExpectedProfitSorted.end(),
                  [](const auto &a, const auto &b)
                  { return a.second > b.second; });
    }

private:
    // Private constructor
    Data(std::shared_ptr<StochasticPackageQuery> spq)
    {
        this->spq = spq;
        this->DB_optim = spq->tableName;
        this->DB_valid = fmt::format("{}_{}", DB_optim, "validate");
        this->NTuples = pg.getTableSize(spq->tableName);
        this->cntScenarios = pg.getColumnLength(spq->tableName, "profit");
        deb("I am creating data instance");
    }

    // static instance
    static std::unique_ptr<Data> instance;
};

// Definition of static member
inline std::unique_ptr<Data> Data::instance = nullptr;
