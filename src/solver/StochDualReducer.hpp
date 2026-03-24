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
    StochDualReducer(std::shared_ptr<StochasticPackageQuery> spq = nullptr,
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

        auto &scenarios = data.stochAttrs["profit"];
        cout<<"In constructor SDR"<<endl;
        for(int j = 0; j < cntScenarios; j++)
        {
            cout<<scenarios[0][j]<<" ";
        }
        cout<<endl;
    }
    
    vector<double> solvelcvarRS(shared_ptr<StochasticPackageQuery>spq, FormulateOptions &formOptions);
};
