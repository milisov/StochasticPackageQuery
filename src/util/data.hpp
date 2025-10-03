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
#pragma once

class Data
{
public:
    // Singleton access
    static Data& getInstance() {
        if (!instance) {
            throw std::runtime_error("Data::init() must be called before getInstance()");
        }
        return *instance;
    }

    static void init(std::shared_ptr<StochasticPackageQuery> spq) {
        if (!instance) {
            instance = std::unique_ptr<Data>(new Data(spq));
        }
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

    void fetchData()
    {
        fetchColumns(spq);
        auto dataColumns = pg.getColumns(spq->tableName);
        std::vector<std::string> cols(columns.begin(), columns.end());
        std::string selectcols = fmt::format("{}", fmt::join(cols, ", "));

        std::string sql = fmt::format("SELECT {} FROM \"{}\" ORDER BY id", selectcols, spq->tableName);

        SingleRow sr(sql);
        int id = 0;
        while (sr.fetchRow())
        {
            for (int i = 0; i < cols.size(); i++)
            {
                std::string colName = cols[i];
                Column colType = dataColumns[colName];
                if (colType == numeric_type)
                {
                    detAttrs[colName].push_back(sr.getNumeric(i));
                }
                else
                {
                    std::vector<double> values;
                    sr.getArray(i, values);
                    stochAttrs[colName].push_back(values);
                    if (statColumns.count(colName) > 0)
                    {
                        double avg = computeAverage(values);
                        stockExpectedProfit.push_back(avg);
                        stockExpectedProfitSorted.emplace_back(id, avg);
                    }
                }
            }
            id++;
        }
        std::sort(stockExpectedProfitSorted.begin(), stockExpectedProfitSorted.end(),
                  [](const auto &a, const auto &b)
                  { return a.second > b.second; });

        //deb(stockExpectedProfitSorted, stockExpectedProfitSorted.size());
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
    }

    // static instance
    static std::unique_ptr<Data> instance;
};

// Definition of static member
inline std::unique_ptr<Data> Data::instance = nullptr;
