#include <iostream>
#include "util/data.hpp"
#include "solver/solversCore.hpp"
#include "solver/SummarySearch.hpp"
#include "spq/sdrformulator.hpp"
#include "spq/ssformulator.hpp"
#include <chrono>
#pragma once

struct FormHelper
{
    SDRFormulator formulatorDet;
    SDRFormulator formulatorCVaR;
    GRBModel modelDetLP;
    GRBModel modelCVaRLP;
};


class StochDualReducer : public Solver
{
private:
    double W0;
    double epsilonQ;

public:
    Data &data;
    double epsilon;
    std::string SPQPath;

    // Constructor
    StochDualReducer(int M = 1e4,
                    std::shared_ptr<StochasticPackageQuery> spq = nullptr,
                    double epsilon = 1e-5) : data(Data::getInstance())
    {
        this->M = M;
        this->spq = spq;
        this->DB_optim = spq->tableName;
        this->DB_valid = fmt::format("{}_{}", DB_optim, "validate");
        this->NTuples = pg.getTableSize(spq->tableName);
        this->cntScenarios = pg.getColumnLength(spq->tableName, "profit");
        this->probConstCnt = countProbConst(spq);
        this->epsilon = epsilon;
    }
    
    SolutionMetadata<int> stochDualReducer(std::shared_ptr<StochasticPackageQuery> spq,
                                           FormulateOptions &formOptions,
                                           std::map<std::string, Option> &curveFitOptions,                                                
                                           const std::chrono::steady_clock::time_point& start_time,
                                           double timeout_seconds);

    
    vector<double> stage1(shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions, GRBModel &modelDetLP, GRBModel &modelCVaRLP);

    vector<int> stage2(std::shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions, 
                       SDRFormulator &formulatorDet,
                       SDRFormulator &formulatorCVaR,
                       GRBModel &modelDetLP, 
                       GRBModel &modelCVaRLP,
                       vector<double> &solLP1);

    SolutionMetadata<int> stage3(shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions, map<string, Option> &curveFitOptions, int z, 
                                               const std::chrono::steady_clock::time_point& start_time,
                                               double timeout_seconds);
};
