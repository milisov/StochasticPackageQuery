#include <iostream>
#include "solver/solversCore.hpp"
#include "spq/naiveformulator.hpp"
#include "util/udebug.hpp"
#include "spq/spq.hpp"
#include "core/stat.hpp"
#include "spq/bounder.hpp"
#include "core/checker.hpp"
#include "solver/starter.hpp"
#include "util/uio.hpp"
#include "spq/cons.hpp"
#include <fmt/core.h>
#include "gurobi_c++.h"
#pragma once
using namespace std;


class Solution
{
public:
    std::vector<int>x;
    bool isFeasible;
    double W_q;
    double Runtime;
    Solution(const std::vector<int>& xInit, double W_qinit, bool isFeasibleInit): x(xInit), W_q(W_qinit), isFeasible(isFeasibleInit) {}
};

class Naive : public Solver {
private:
    int M_hat;
    std::vector<GRBVar> yy;
    std::vector<GRBGenConstr> genCon;
    std::vector<GRBConstr> sumyCon;

public:
    Data& data;
    Naive(std::shared_ptr<StochasticPackageQuery> spq = nullptr) : data(Data::getInstance())
    {
        this->spq = spq;
        this->DB_optim = spq->tableName;
        this->DB_valid = fmt::format("{}_{}", DB_optim, "validate");
        this->NTuples = pg.getTableSize(spq->tableName);
        this->cntScenarios = pg.getColumnLength(spq->tableName, "profit");
        this->probConstCnt = countProbConst(spq);
        std::cout << "Success Constructor (Naive)" << std::endl;
    }

    template <typename T>
    SolutionMetadata<T> solveNaive(shared_ptr<StochasticPackageQuery> spq, FormulateOptions& formOptions)
    {
        deb("solveNaive");
        NaiveFormulator formulator(spq);
        GRBModel model = formulator.formulate(spq, formOptions);
        deb("here");
        vector<int>x;
        initializeVector(x,NTuples,0);
        SolveOptions options;
        options.reduced = formOptions.reduced;
        options.reducedIds = formOptions.reducedIds;
        options.computeActiveness = formOptions.computeActiveness;
        solve(model,x, options);
        validate(model, x, spq, options);
        SolutionMetadata<T> sol;
        if(isFeasible(r))
        {
            sol.x = x;
            sol.isFeasible = true; 
        }else
        {
            sol.x = x;
            sol.isFeasible = false;
        }
        return sol;
    }
};


