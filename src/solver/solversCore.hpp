#include <iostream>
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
#include "gurobi_c++.h"
#include <fmt/ranges.h>
#include <boost/algorithm/string/join.hpp>
#include <util/data.hpp>
#include <fstream>
#pragma once

struct SolveOptions
{
    int Z = 0;
    bool computeActiveness = false;
    bool reduced = false;
    vector<int> reducedIds;
};

// A small metadata helper
class BinarySearchMetadata
{
public:
    double low;
    double high;
    double alpha;
    BinarySearchMetadata(double lowInit, double highInit, double alphaInit)
        : low(lowInit), high(highInit), alpha(alphaInit) {}
};

template <typename T>
class SolutionMetadata
{
public:
    std::vector<T> x;
    double w;
    double objConsValue;
    double epsilon;
    bool isFeasible;
    bool isOptimal;
    int binarySearchSteps = -1;
    double bestEps;
    int Z = 0;
    int qSz = 0;
    double timeStage1 = 0.0;
    double timeStage2 = 0.0;
    double bestRk = -1e7;
    std::vector<pair<int, double>> bestPosActivenessRS;
    std::vector<pair<int, double>> bestNegActivenessRS;
    // Constructors
    SolutionMetadata() : w(0), epsilon(0), isFeasible(false), isOptimal(false) {}
    SolutionMetadata(const std::vector<T> &xInit, double wInit, double epsilonInit, bool isFeasibleInit)
        : x(xInit), w(wInit), epsilon(epsilonInit), isFeasible(isFeasibleInit), isOptimal(false) {}

    void setSolution(const std::vector<T> &xInit, double w, double epsilon, bool isFeasible, bool isOptimal, int Z)
    {
        x = xInit;
        this->w = w;
        this->epsilon = epsilon;
        this->isFeasible = isFeasible;
        this->isOptimal = isOptimal;
        this->Z = Z;
    }
};

class Solver
{
public:
    Data &data;
    double timeFetch = 0.0;
    double timeSolve = 0.0;
    Solver() : data(Data::getInstance()) {};
    PgManager pg;
    int M;
    std::shared_ptr<StochasticPackageQuery> spq;
    std::string DB_optim;
    std::string DB_valid;
    int NTuples;
    int cntScenarios;
    int probConstCnt;
    std::vector<double> r;
    // for each prob constraint, we carry a value per scenario
    std::vector<std::vector<pair<int, double>>> innerConstraints;
    std::vector<pair<int, double>> posActivenessRS;
    std::vector<pair<int, double>> negActivenessRS;
    double W_q;

    template <typename T>
    void validate(GRBModel &model, std::vector<T> &x,
                  std::shared_ptr<StochasticPackageQuery> spq,
                  SolveOptions &options);

    template <typename T>
    void solve(GRBModel &model, std::vector<T> &x, SolveOptions &options);

    template <typename T>
    double calculateObj(std::vector<T> &x,
                        std::vector<int> &selectIds,
                        std::shared_ptr<Objective> obj);

    template <typename T>
    double calculateCntObj(std::vector<T> &x,
                           std::vector<int> &selectIds,
                           std::shared_ptr<CountObjective> cntObj);

    template <typename T>
    double calculateSumObj(std::vector<T> &x,
                           std::vector<int> &selectIds,
                           std::shared_ptr<AttrObjective> attrObj);

    template <typename T>
    double calculateExpSumObj(std::vector<T> &x,
                              std::shared_ptr<AttrObjective> attrObj);

    bool isFeasible(std::vector<double> &r);

    double calculateRk(std::shared_ptr<ProbConstraint> probCon,
                       int satisfied,
                       int numScenarios,
                       std::shared_ptr<StochasticPackageQuery> spq);

    double calculateEpsilonQ(shared_ptr<StochasticPackageQuery> spq, double W_q, double W0);
    BinarySearchMetadata guessOptimalConservativenessBinarySearch(double low,
                                                                  double high,
                                                                  double rk);

    int countSatisfied(int numScenarios, 
                       std::vector<pair<int, double>> &innerConst, 
                       shared_ptr<ProbConstraint> probCon, 
                       shared_ptr<StochasticPackageQuery> spq);

    int countProbConst(std::shared_ptr<StochasticPackageQuery> spq);
};

inline void printVariableNames(GRBModel &model)
{
    try
    {
        // Get all the variables in the model
        int numVars = model.get(GRB_IntAttr_NumVars);
        GRBVar *vars = model.getVars();

        // Print the names of each variable
        for (int i = 0; i < numVars; i++)
        {
            GRBVar var = vars[i];
            std::string varName = var.get(GRB_StringAttr_VarName);
            std::cout << "Variable " << i << ": " << varName << std::endl;
        }

        delete[] vars;
    }
    catch (GRBException &e)
    {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
    catch (...)
    {
        std::cout << "Exception during printing variable names" << std::endl;
    }
}

template <typename T>
inline void findNonzero(std::vector<int> &selectIds, std::vector<T> &x)
{
    for (int i = 0; i < x.size(); i++)
    {
        if (x[i] != 0)
        {
            selectIds.push_back(i + 1);
        }
    }
}

template <typename T>
inline void initializeVector(std::vector<T> &v, int sz, T init)
{
    for (int i = 0; i < sz; i++)
    {
        v.push_back(init);
    }
}

// for each of the Stochastic Constraints we need to check how many are satisfied
// get the string sign from each prob constraint and perform the right operation
// the value of v is retrieved from the spq
inline int Solver::countSatisfied(int numScenarios, std::vector<pair<int, double>> &innerConst, shared_ptr<ProbConstraint> probCon, shared_ptr<StochasticPackageQuery> spq)
{
    int satisfied = 0;
    double epsilon = 1e-7;
    for (int i = 0; i < numScenarios; i++)
    {
        if (probCon->vsign == Inequality::gteq) // >=
        {
            if (innerConst[i].second >= spq->getValue(probCon->v) - epsilon)
            {
                posActivenessRS.emplace_back(i,innerConst[i].second - spq->getValue(probCon->v));
                satisfied++;
            }else
            {
                negActivenessRS.emplace_back(i,innerConst[i].second - spq->getValue(probCon->v));
            }
        }
        else
        {
            if (innerConst[i].second <= spq->getValue(probCon->v) + epsilon) // <=
            {
                satisfied++;
                posActivenessRS.emplace_back(i, spq->getValue(probCon->v) - innerConst[i].second);
            }else
            {
                negActivenessRS.emplace_back(i, spq->getValue(probCon->v) - innerConst[i].second);
            }
        }
    }
    return satisfied;
}

inline bool Solver::isFeasible(vector<double> &r)
{
    for (int k = 0; k < r.size(); k++)
    {
        if (r[k] < 0)
        {
            return false;
        }
    }
    return true;
}

inline double Solver::calculateRk(shared_ptr<ProbConstraint> probCon, int satisfied, int numScenarios, shared_ptr<StochasticPackageQuery> spq)
{
    double rk;
    // cout << "Satisfied: " << satisfied << endl;
    if (probCon->psign == Inequality::gteq)
    {
        rk = (double)satisfied / numScenarios - spq->getValue(probCon->p);
        // cout << "rk-gteq" << " " << rk << endl;
    }
    else
    {
        rk = (1 - (double)satisfied / numScenarios) - (1 - spq->getValue(probCon->p));
    }

    return rk;
}

template <typename T>
inline double Solver::calculateCntObj(std::vector<T> &x, std::vector<int> &selectIds, shared_ptr<CountObjective> cntObj)
{
    double sum = 0;
    for (int i = 0; i < x.size(); i++)
    {
        sum += x[i];
    }
    return sum;
}

template <typename T>
inline double Solver::calculateSumObj(std::vector<T> &x, std::vector<int> &selectIds, shared_ptr<AttrObjective> attrObj)
{
    double sum = 0;
    auto &detAttr = data.detAttrs[attrObj->obj];
    for (int i = 0; i < NTuples; i++)
    {
        sum += x[i] * detAttr[i];
    }
    return sum;
}

template <typename T>
inline double Solver::calculateExpSumObj(std::vector<T> &x, shared_ptr<AttrObjective> attrObj)
{
    double sum = 0;
    for (int i = 0; i < x.size(); i++)
    {
        sum += x[i] * data.stockExpectedProfit[i];
    }
    return sum;
}

template <typename T>
inline double Solver::calculateObj(std::vector<T> &x, std::vector<int> &selectIds, shared_ptr<Objective> obj)
{
    double w = 0;
    shared_ptr<CountObjective> cntObj = getCount(obj);
    if (cntObj)
    {
        w = calculateCntObj(x, selectIds, cntObj);
        return w;
    }

    shared_ptr<AttrObjective> attrObj;
    bool isDet = isDeterministic(obj, attrObj);
    if (isDet && attrObj->objType == numeric_type)
    {
        w = calculateSumObj(x, selectIds, attrObj);
        return w;
    }

    shared_ptr<AttrObjective> attrObj2;
    bool isDet2 = isDeterministic(obj, attrObj2);
    if (isDet2 && attrObj2->objType == array_type)
    {
        w = calculateExpSumObj(x, selectIds, attrObj2);
        return w;
    }
    return -1;
}

template <typename T>
inline void Solver::validate(GRBModel &model, std::vector<T> &x, shared_ptr<StochasticPackageQuery> spq, SolveOptions &options)
{
    // empty previous surplus vector
    // cout << "Validating!" << endl;
    // clear r from previous iteration -> holds all of the rk values for each iteration
    this->r.clear();
    // we need to clear the innerConstraints in order to update them with new values
    this->innerConstraints.clear();
    this->posActivenessRS.clear();
    this->negActivenessRS.clear();

    // because we get row by row we need a way to compute the innerConst for all scenarios in the same time
    int cons_num = spq->cons.size();
    // for all of the constrants, check if they are probConstraints
    // query the optim table
    // for each tuple update innerConst
    for (int i = 0; i < cons_num; i++)
    {
        shared_ptr<ProbConstraint> probCon;
        shared_ptr<AttrConstraint> attrCon;
        bool isStoch = isStochastic(spq->cons[i], probCon, attrCon);
        if (isStoch)
        {
            std::vector<pair<int, double>> innerConst;
            for (int i = 0; i < cntScenarios; i++)
            {
                innerConst.emplace_back(i, 0.0);
                //activenessRS.emplace_back(i, 0.0);
            }
            auto &scenarios = data.stochAttrs[attrCon->attr];
            if (options.reduced)
            {
                for (int i = 0; i < options.reducedIds.size(); i++)
                {
                    int id = options.reducedIds[i] - 1;
                    for (int j = 0; j < cntScenarios; j++)
                    {
                        innerConst[j].second += x[id] * scenarios[id][j];
                    }
                }
            }
            else
            {
                for (int i = 0; i < NTuples; i++)
                {
                    for (int j = 0; j < cntScenarios; j++)
                    {
                        innerConst[j].second += x[i] * scenarios[i][j];
                    }
                }
            }
            // push innerConst to innerConstraints
            this->innerConstraints.push_back(innerConst);
            int satisfied = countSatisfied(cntScenarios, innerConst, probCon, spq);
            double rk = calculateRk(probCon, satisfied, cntScenarios, spq);
            deb(rk);
            this->r.push_back(rk);
        }
    }
    this->W_q = model.get(GRB_DoubleAttr_ObjVal);
    deb(W_q);
}

// solve() -> model.optimize(); then translate the x[i] GRBVar vector into a x vector so you can manipulate it.
template <typename T>
inline void Solver::solve(GRBModel &model, std::vector<T> &x, SolveOptions &options)
{
    int numConstraints = model.get(GRB_IntAttr_NumConstrs);
    int numVars = model.get(GRB_IntAttr_NumVars);
    gpro.clock("time1");
    model.set(GRB_IntParam_Presolve, 0);
    model.optimize();
    gpro.stop("time1");
    int status = model.get(GRB_IntAttr_Status);

    if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE || model.get(GRB_IntAttr_Status) == GRB_INF_OR_UNBD)
    {
        cout << "Initial model is infeasible or unbounded" << std::endl;
        x.clear();
        return;
    }
    else if (status == GRB_TIME_LIMIT)
    {
        std::cout << "Time Limit Expried" << std::endl;
        x.clear();
        return;
    }
    else
    {
         std::cout << "Model optimization ended with status: "<< model.get(GRB_IntAttr_Status) << std::endl;
    }

    cout << "Solution Found" << endl;

    GRBVar *xx = model.getVars();
    if (options.reduced)
    {
        for (int i = 0; i < options.reducedIds.size(); i++)
        {
            int id = options.reducedIds[i] - 1;
            x[id] = static_cast<T>(xx[i].get(GRB_DoubleAttr_X));
        }
    }
    else
    {
        for (int i = 0; i < this->NTuples; i++)
        {
            x[i] = static_cast<T>(xx[i].get(GRB_DoubleAttr_X));
        }
    }
}

inline BinarySearchMetadata Solver::guessOptimalConservativenessBinarySearch(double low, double high, double rk)
{
    double alpha = low + (high - low) / 2;
    BinarySearchMetadata metadata(0.0, 0.0, 0.0);
    if (rk < 0)
    {
        // cout << "USE MORE CONSERVATIVE" << endl;
        //  the solution is infeasible -> use more conservative summary
        metadata.low = alpha;
        metadata.high = high;
        metadata.alpha = alpha;
    }
    else
    {
        // the solution is feasible but suboptimal or the system is infeasible -> use less conservative summary
        // cout << "USE LESS CONSERVATIVE" << endl;
        metadata.low = low;
        metadata.high = alpha;
        metadata.alpha = alpha;
    }

    return metadata;
}

inline int Solver::countProbConst(shared_ptr<StochasticPackageQuery> spq)
{
    int cnt = 0;
    int cons_num = spq->cons.size();

    for (int i = 0; i < cons_num; i++)
    {
        shared_ptr<ProbConstraint> probCon;
        shared_ptr<AttrConstraint> attrCon;
        bool isStoch = isStochastic(spq->cons[i], probCon, attrCon);
        if (isStoch)
        {
            cnt++;
        }
    }

    return cnt;
}

inline double Solver::calculateEpsilonQ(shared_ptr<StochasticPackageQuery> spq, double W_q, double W0)
{
    double epsilonQ;
    if (spq->obj->objSense == maximize)
    {
        epsilonQ = W0 / W_q - 1;
    }
    else
    {
        epsilonQ = (W_q / W0) - 1;
    }
    return epsilonQ;
}