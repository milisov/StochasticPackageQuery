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

inline int computeIterationsToNTuples(int qTarget_initial, int NTuples, int qInc = 500)
{
    if (qTarget_initial >= NTuples)
        return 1; // At least 1 iteration
    double remaining = NTuples - qTarget_initial;
    return (int)std::ceil(std::log2(remaining / qInc + 1)) + 1;
}

inline int computeBinarySearchIterations(double T, double t, double c)
{
    deb(T, t, c);
    int d_max = (int)std::floor((1.25 - 1) * (T + t) / (t * (pow(1.25, c + 1) - 1)));
    d_max = std::max(d_max, 1);
    return d_max;
}

inline int computeBinarySearchIterationsForN(int cntScenarios, int NTuples, int N_start, SolveOptions &options)
{
    int d_max = static_cast<int>(options.hyperParams["cap"]);
    if (cntScenarios <= 10)
        return 1;
    if (cntScenarios >= 10000)
        return d_max;

    double log_t = (std::log2((double)cntScenarios) - std::log2(10.0)) / (std::log2((double)10000) - std::log2(10.0));

    return (int)std::floor(1 + (d_max - 1) * log_t);
}

inline double computeEpsilon(int cntScenarios, int NTuples, SolveOptions &options)
{
    int d = computeBinarySearchIterationsForN(cntScenarios, NTuples, static_cast<int>(options.hyperParams["reducedTuples"]), options);
    return std::sqrt(2.0 * M_PI) / std::pow(2.0, d);
}

class Solution
{
public:
    std::vector<int> x;
    bool isFeasible;
    double W_q;
    double Runtime;
    Solution(const std::vector<int> &xInit, double W_qinit, bool isFeasibleInit) : x(xInit), W_q(W_qinit), isFeasible(isFeasibleInit) {}
};

class Naive : public Solver
{
private:
    int M_hat;
    std::vector<GRBVar> yy;
    std::vector<GRBGenConstr> genCon;
    std::vector<GRBConstr> sumyCon;

public:
    Data &data;
    Naive(std::shared_ptr<StochasticPackageQuery> spq = nullptr, std::string validateTableNameBase = "") : data(Data::getInstance())
    {
        this->spq = spq;
        this->DB_optim = spq->tableName;
        this->DB_valid = validateTableNameBase;
        this->NTuples = pg.getTableSize(spq->tableName);
        this->cntScenarios = pg.getColumnLength(spq->tableName, "profit");
        this->probConstCnt = countProbConst(spq);
        spq->setTableName(DB_valid);
        this->checker = std::make_unique<SPQChecker>(this->spq);
        spq->setTableName(DB_optim);
        std::cout << "Success Constructor (Naive)" << std::endl;
    }

    // Formulate and solve a single iteration
    template <typename T>
    SolutionMetadata<T> naiveSolve(shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions, SolveOptions &solveOptions);

    // Formulate and solve the ILP with given formOptions and solveOptions
    template <typename T>
    SolutionMetadata<T> formulateAndSolve(shared_ptr<StochasticPackageQuery> spq,
                                          FormulateOptions &formOptions,
                                          SolveOptions &solveOptions);

    // Out-of-sample validation using checker
    template <typename T>
    void validateOutOfSample(SolutionMetadata<T> &sol,
                             FormulateOptions &formOptions,
                             SolveOptions &solveOptions);
};

template <typename T>
void Naive::validateOutOfSample(SolutionMetadata<T> &sol,
                                FormulateOptions &formOptions,
                                SolveOptions &solveOptions)
{
    if (sol.x.empty())
        return;

    gpro.stop("effectiveRuntime");

    vector<double> feasibilities;
    vector<double> surpluses;
    SolType res;
    for (int i = 0; i < NTuples; i++)
    {
        if (sol.x[i] > 0)
        {
            res[i + 1] = sol.x[i];
        }
    }
    bool deterFeasible, probFeasible;
    bool feas = checker->feasible(res, feasibilities, surpluses, deterFeasible, probFeasible);
    double validObj = checker->getObjective(res);

    for (int i = 0; i < surpluses.size(); i++)
    {
        if (surpluses[i] < 0)
        {
            probFeasible = false;
        }
    }

    // Order: runtime, optimFeas, optimSurplus, optimObj, deterFeasible, probFeasible, validFeas, validSurplus, validObj, qSz, Z
    double timeStamp;
    if (solveOptions.hyperParams["benchmark"] == 2)
    {
        timeStamp = gpro.getTime("effectiveRuntime") - gpro.getTime("fetchingToEndOfFirstStage");
    }
    else
    {
        timeStamp = gpro.getTime("effectiveRuntime");
    }

    sol.solutions.push_back(make_tuple(
        timeStamp,
        this->satisfiedScenarios,
        this->r,
        this->W_q,
        deterFeasible,
        probFeasible,
        feasibilities,
        surpluses,
        validObj,
        solveOptions.reducedIds.size(),
        formOptions.Z));

    gpro.clock("effectiveRuntime");
}

template <typename T>
SolutionMetadata<T> Naive::formulateAndSolve(shared_ptr<StochasticPackageQuery> spq,
                                             FormulateOptions &formOptions,
                                             SolveOptions &solveOptions)
{
    formOptions.kMADValues.clear();
    SolutionMetadata<int> bestSol;
    vector<double> low;
    vector<double> high;
    double epsilon = computeEpsilon(cntScenarios, NTuples, solveOptions);

    if (low.empty())
    {
        for (int i = 0; i < probConstCnt; i++)
        {
            low.push_back(0.0);
            high.push_back(std::sqrt(M_PI * 2));
            formOptions.kMADValues.push_back((low.back() + high.back()) / 2.0);
        }
    }

    bool breakLoop = false;
    while (true)
    {
        if (breakLoop)
        {
            break;
        }

        // Check time budget before formulating
        if (solveOptions.timeBudgeted)
        {
            gpro.stop("effectiveRuntime");
            double remainingTime = solveOptions.timeBudget - (gpro.getTime("effectiveRuntime") / 1000.0);
            gpro.clock("effectiveRuntime");
            if (remainingTime <= 0)
            {
                cout << "Time budget exceeded, exiting binary search loop" << endl;
                break;
            }
        }

        deb(low, high, formOptions.kMADValues);
       
        SolutionMetadata<int> sol;
        NaiveFormulator formulator(spq);
        GRBModel model = formulator.formulate(spq, formOptions);
        vector<int> x;
        initializeVector(x, NTuples, 0);
        {
            int terminateTuples = solveOptions.hyperParams.count("terminateTuples") ? static_cast<int>(solveOptions.hyperParams["terminateTuples"]) : NTuples;
            bool atTerminationSize = ((int)formOptions.reducedIds.size() >= NTuples) || ((int)formOptions.reducedIds.size() >= terminateTuples);
            solveOptions.enableRelaxation = atTerminationSize;
        }
        solve(model, x, solveOptions);
        if (!x.empty())
        {
            validate(model, x, spq, solveOptions);
            sol.x = x;
            sol.w = this->W_q;
            sol.bestRk = this->r;
            sol.isFeasible = isFeasible(this->r);
            validateOutOfSample(sol, formOptions, solveOptions);

            // Always add solutions to tracking list
            bestSol.solutions.insert(bestSol.solutions.end(), sol.solutions.begin(), sol.solutions.end());

            int sense = model.get(GRB_IntAttr_ModelSense);
            bool currentIsBetter = isCurrentBetterThanBest(this->r, this->W_q, bestSol, sense);
            if (currentIsBetter)
            {
                bestSol.x = sol.x;
                bestSol.w = sol.w;
                bestSol.bestRk = sol.bestRk;
                bestSol.isFeasible = sol.isFeasible;
                bestSol.qSz = formOptions.reducedIds.size();
                bestSol.qScenarios = formOptions.Z;
            }

            for (int i = 0; i < this->r.size(); i++)
            {
                if (this->r[i] < 0)
                {
                    low[i] = formOptions.kMADValues[i];
                }
                else
                {
                    high[i] = formOptions.kMADValues[i];
                }
            }

            if (solveOptions.foundSolutionUsingGurobiRelax)
            {
                breakLoop = true;
                break;
            }
        }
        else
        {
            for (int i = 0; i < probConstCnt; i++)
            {
                low[i] = formOptions.kMADValues[i];
            }
        }

        // Check convergence outside the if/else block
        for (int i = 0; i < probConstCnt; i++)
        {
            deb(i, high[i] - low[i], epsilon, formOptions.kMADValues[i]);
            if (high[i] - low[i] <= epsilon)
            {
                cout << "HERE I BREAK" << endl;
                breakLoop = true;
                break;
            }
        }

        for (int i = 0; i < probConstCnt; i++)
        {
            formOptions.kMADValues[i] = (low[i] + high[i]) / 2.0;
        }
    }

    return bestSol;
}

template <typename T>
SolutionMetadata<T> Naive::naiveSolve(shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions, SolveOptions &solveOptions)
{
    // deb(formOptions.reducedIds);
    NaiveFormulator formulator(spq);
    GRBModel model = formulator.formulate(spq, formOptions);
    SolutionMetadata<T> sol;
    vector<T> x;
    initializeVector(x, NTuples, T(0));
    solve(model, x, solveOptions);
    if (!x.empty())
    {
        validate(model, x, spq, solveOptions);
        sol.x = x;
        sol.w = this->W_q;
        sol.bestRk = this->r;
        sol.isFeasible = isFeasible(this->r);
        validateOutOfSample(sol, formOptions, solveOptions);
    }
    return sol;
}
