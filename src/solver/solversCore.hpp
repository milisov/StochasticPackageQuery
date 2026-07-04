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
#include <tuple>
#include <set>
#pragma once

struct SolveOptions
{
    int Z = 0;
    double timeout_seconds = -1.0;
    bool timeBudgeted = false;
    double timeBudget = 0.0;
    bool includeObjectiveFunction = true;
    bool reduced = false;
    vector<int> reducedIds;
    bool enableRelaxation = false;
    bool enableRelaxationLP = false;
    bool foundSolutionUsingGurobiRelax = false;
    std::vector<std::set<int>> reducedScenariosPerConstraint;

    map<string, double> hyperParams;
    unsigned int randomSeed = 42;
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
    int qScenarios = 0;
    bool terminate = false;
    double timeStage1 = 0.0;
    double timeStage2 = 0.0;
    vector<double> bestRk;
    vector<std::tuple<double, vector<double>, vector<double>, double, bool, bool, vector<double>, vector<double>, double, int, int>> solutions; // (runtime, optimFeas, optimSurplus, optimObj, deterFeasible, probFeasible, validFeas, validSurplus, validObj, qSz, Z)
    // double bestRk = -1e7;
    vector<vector<pair<int, double>>> bestPosActiveness;
    vector<vector<pair<int, double>>> bestNegActiveness;
    vector<vector<pair<int, double>>> bestActiveness;
    vector<vector<pair<int, double>>> bestActivenessNoAbsValue;

    vector<double> meanPerVaR;
    vector<double> variancePerVaR;
    vector<double> minInnerConstPerVaR;
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
    Solver() : data(Data::getInstance()) { cout << "I am getting an instsance of Data" << endl; };
    PgManager pg;
    int M;
    std::shared_ptr<StochasticPackageQuery> spq;
    std::unique_ptr<SPQChecker> checker;
    std::string DB_optim;
    std::string DB_valid;
    int NTuples;
    int cntScenarios;
    int probConstCnt;
    std::vector<double> r; // surplus
    std::vector<double> satisfiedScenarios;
    // for each prob constraint, we carry a value per scenario
    vector<vector<pair<int, double>>> innerConstraints;
    vector<vector<pair<int, double>>> activeness;
    vector<vector<pair<int, double>>> activenessNoAbsValue;
    double W_q;
    vector<double> meanPerVaR;
    vector<double> variancePerVaR;
    vector<double> minInnerConstPerVaR;

    template <typename T>
    vector<std::tuple<int, double, double>> getZeroTuplesRanking(std::vector<T> &x, std::shared_ptr<StochasticPackageQuery> spq, SolveOptions &options);

    template <typename T>
    vector<pair<int, T>> getNonZeroTuplesRanking(std::vector<T> &x);

    template <typename T>
    void validate(GRBModel &model, std::vector<T> &x,
                  std::shared_ptr<StochasticPackageQuery> spq,
                  SolveOptions &options);

    template <typename T>
    void solve(GRBModel &model, std::vector<T> &x, SolveOptions &options);

    template <typename T>
    int getNonZeroCount(std::vector<T> &x);

    template <typename T>
    double calculateObj(vector<pair<int, T>> &package,
                        std::shared_ptr<Objective> obj);

    template <typename T>
    double calculateCntObj(vector<pair<int, T>> &package,
                           std::shared_ptr<CountObjective> cntObj);

    template <typename T>
    double calculateSumObj(vector<pair<int, T>> &package,
                           std::shared_ptr<AttrObjective> attrObj);

    template <typename T>
    double calculateExpSumObj(vector<pair<int, T>> &package,
                              std::shared_ptr<AttrObjective> attrObj);

    template <typename T>
    vector<pair<int, T>> getPackage(std::vector<T> &x);

    bool isFeasible(std::vector<double> &r);

    double calculateRk(std::shared_ptr<ProbConstraint> probCon,
                       int satisfied,
                       int numScenarios,
                       std::shared_ptr<StochasticPackageQuery> spq);

    double calculateEpsilonQ(shared_ptr<StochasticPackageQuery> spq, double W_q, double W0);
    BinarySearchMetadata guessOptimalConservativenessBinarySearch(double low,
                                                                  double high,
                                                                  double rk);

    int countSatisfied(SolveOptions &options,
                       int numScenarios,
                       std::vector<pair<int, double>> &innerConst,
                       shared_ptr<ProbConstraint> probCon,
                       shared_ptr<StochasticPackageQuery> spq);

    int countProbConst(std::shared_ptr<StochasticPackageQuery> spq);

    template <typename T>
    bool isCurrentBetterThanBest(vector<double> r, double W_q, SolutionMetadata<T> &bestSol, int sense);
};

template <typename T>
inline vector<pair<int, T>> Solver::getPackage(std::vector<T> &x)
{
    std::vector<pair<int, T>> package;
    for (int i = 0; i < x.size(); i++)
    {
        if (x[i] > 0)
        {
            package.push_back(make_pair(i, x[i]));
        }
    }
    return package;
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
inline int Solver::getNonZeroCount(std::vector<T> &x)
{
    int cnt = 0;
    for (int i = 0; i < x.size(); i++)
    {
        if (x[i] > 0)
        {
            cnt += 1;
        }
    }
    return cnt;
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
inline int Solver::countSatisfied(SolveOptions &options, int numScenarios, std::vector<pair<int, double>> &innerConst, shared_ptr<ProbConstraint> probCon, shared_ptr<StochasticPackageQuery> spq)
{
    // deb(innerConst);
    int satisfied = 0;
    double epsilon = Config::getInstance()->pt.get<double>("parameters.numeric_eps");
    vector<pair<int, double>> activenessForCons;
    vector<pair<int, double>> activenessForConsNoAbsValue;
    double var = 0.0;
    double minInnerConst = innerConst[0].second;
    for (int i = 0; i < numScenarios; i++)
    {
        double val = innerConst[i].second;
        if (probCon->vsign == Inequality::gteq) // >=
        {
            if (val >= spq->getValue(probCon->v) - epsilon)
            {
                satisfied++;
            }
            activenessForCons.emplace_back(i, abs(val - spq->getValue(probCon->v)));
            activenessForConsNoAbsValue.emplace_back(i, val - spq->getValue(probCon->v));
        }
        else
        {
            if (val <= spq->getValue(probCon->v) + epsilon) // <=
            {
                satisfied++;
            }
            activenessForCons.emplace_back(i, abs(val - spq->getValue(probCon->v)));
            activenessForConsNoAbsValue.emplace_back(i, val - spq->getValue(probCon->v));
        }
    }
    activeness.push_back(activenessForCons);
    activenessNoAbsValue.push_back(activenessForConsNoAbsValue);
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
inline double Solver::calculateCntObj(vector<pair<int, T>> &package, shared_ptr<CountObjective> cntObj)
{
    double sum = 0;
    for (auto &[id, val] : package)
        sum += static_cast<double>(val);
    return sum;
}

template <typename T>
inline double Solver::calculateSumObj(vector<pair<int, T>> &package, shared_ptr<AttrObjective> attrObj)
{
    double sum = 0;
    auto &detAttr = data.detAttrs[attrObj->obj];
    for (auto &[id, val] : package)
        sum += static_cast<double>(val) * detAttr[id];
    return sum;
}

template <typename T>
inline double Solver::calculateExpSumObj(vector<pair<int, T>> &package, shared_ptr<AttrObjective> attrObj)
{
    double sum = 0;
    for (auto &[id, val] : package)
        sum += static_cast<double>(val) * data.stockExpectedProfit[id];
    return sum;
}

template <typename T>
inline double Solver::calculateObj(vector<pair<int, T>> &package, shared_ptr<Objective> obj)
{
    double w = 0;
    shared_ptr<CountObjective> cntObj = getCount(obj);
    if (cntObj)
        return calculateCntObj(package, cntObj);

    shared_ptr<AttrObjective> attrObj;
    bool isDet = isDeterministic(obj, attrObj);
    if (isDet && attrObj->objType == numeric_type)
        return calculateSumObj(package, attrObj);

    shared_ptr<AttrObjective> attrObj2;
    bool isDet2 = isDeterministic(obj, attrObj2);
    if (isDet2 && attrObj2->objType == array_type)
        return calculateExpSumObj(package, attrObj2);

    return -1;
}

template <typename T>
inline void Solver::validate(GRBModel &model, std::vector<T> &x, shared_ptr<StochasticPackageQuery> spq, SolveOptions &options)
{
    vector<pair<int, T>> package = getPackage(x);
    // clear r from previous iteration -> holds all of the rk values for each iteration
    this->r.clear();
    this->satisfiedScenarios.clear();
    // we need to clear the innerConstraints in order to update them with new values
    this->innerConstraints.clear();
    this->activeness.clear();
    this->activenessNoAbsValue.clear();

    // because we get row by row we need a way to compute the innerConst for all scenarios in the same time
    int cons_num = spq->cons.size();
    // for all of the constrants, check if they are probConstraints
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
            }
            auto &scenarios = data.stochAttrs[attrCon->attr];

            for (auto &p : package)
            {
                int id = p.first;
                T val = p.second;
                for (int j = 0; j < cntScenarios; j++)
                {
                    innerConst[j].second += val * scenarios[id][j];
                }
            }
            // push innerConst to innerConstraints
            this->innerConstraints.push_back(innerConst);
            int satisfied = countSatisfied(options, cntScenarios, innerConst, probCon, spq);
            double rk = calculateRk(probCon, satisfied, cntScenarios, spq);
            double satisfiedScenarios = (double)satisfied / cntScenarios;
            deb(rk, satisfiedScenarios);
            this->r.push_back(rk);
            this->satisfiedScenarios.push_back(satisfiedScenarios);
        }
    }
    if (options.includeObjectiveFunction)
    {
        this->W_q = model.get(GRB_DoubleAttr_ObjVal);
    }
    else
    {
        this->W_q = getNonZeroCount(x);
        cout << "Non-zero count objective" << " " << W_q << endl;
    }
}

template <typename T>
inline vector<std::tuple<int, double, double>> Solver::getZeroTuplesRanking(std::vector<T> &x, std::shared_ptr<StochasticPackageQuery> spq, SolveOptions &options)
{
    vector<pair<int, T>> package = getPackage(x);

    // Find element with minimum value in the package
    double minValue = static_cast<double>(package[0].second);
    int minId = package[0].first;
    for (auto &p : package)
    {
        double val = static_cast<double>(p.second);
        if (val < minValue) 
        {
            minValue = val; 
            minId = p.first; 
        }
    }

    deb(minValue, minId);

    int cons_num = spq->cons.size();
    int conOrder = 0;
    vector<vector<pair<int, double>>> innerCons;
    vector<shared_ptr<ProbConstraint>> probCons;
    vector<shared_ptr<AttrConstraint>> attrCons;

    for (int i = 0; i < cons_num; i++)
    {
        shared_ptr<ProbConstraint> probCon;
        shared_ptr<AttrConstraint> attrCon;
        if (!isStochastic(spq->cons[i], probCon, attrCon)) continue;

        auto &reducedScenarios = options.reducedScenariosPerConstraint[conOrder];
        auto &scenarios = data.stochAttrs[attrCon->attr];

        vector<pair<int, double>> innerConst;
        for (int k : reducedScenarios)
        {
            innerConst.emplace_back(k, 0.0);
        }

        for (auto &p : package)
        {
            if (p.first == minId) continue;
            double val = static_cast<double>(p.second);
            int j = 0;
            for (int k : reducedScenarios)
            {
                innerConst[j++].second += val * scenarios[p.first][k];
            }
        }

        innerCons.push_back(innerConst);
        probCons.push_back(probCon);
        attrCons.push_back(attrCon);
        conOrder++;
    }
    int numProbCons = conOrder;

    // Base objective from package \ {minId}
    vector<pair<int, T>> basePackage;
    for (auto &p : package)
    {
        if (p.first != minId)
        {
            basePackage.push_back(p);
        }
    }
    double baseObj = calculateObj(basePackage, spq->obj);

    double epsilon = Config::getInstance()->pt.get<double>("parameters.numeric_eps");
    vector<std::tuple<int, double, double>> ranking;

    for (int i = 0; i < NTuples; i++)
    {
        if (x[i] != 0) continue;

        double worstRk = 1e7;
        for (int c = 0; c < numProbCons; c++)
        {
            auto &probCon = probCons[c];
            auto &scenarios = data.stochAttrs[attrCons[c]->attr];
            int Z = (int)innerCons[c].size();
            int satisfied = 0;
            for (auto &[k, baseInner] : innerCons[c])
            {
                double newInner = baseInner + minValue * scenarios[i][k];
                if (probCon->vsign == Inequality::gteq)
                {
                    satisfied += (newInner >= spq->getValue(probCon->v) - epsilon) ? 1 : 0;
                }
                else
                {
                    satisfied += (newInner <= spq->getValue(probCon->v) + epsilon) ? 1 : 0;
                }
            }
            double rk = calculateRk(probCon, satisfied, Z, spq);
            if (rk < worstRk)
            {
                worstRk = rk;
            }
        }

        // Objective delta: contribution of adding (i, minValue)
        vector<pair<int, T>> elem = {{i, static_cast<T>(minValue)}};
        double candidateObj = baseObj + calculateObj(elem, spq->obj);

        ranking.emplace_back(i+1, worstRk, candidateObj);
    }

    sort(ranking.begin(), ranking.end(), [](const auto &a, const auto &b)
    {
        if (std::get<1>(a) != std::get<1>(b))
        {
            return std::get<1>(a) > std::get<1>(b);
        }
        return std::get<2>(a) > std::get<2>(b);
    });

    return ranking;
}


template<typename T>
inline vector<pair<int, T>> Solver::getNonZeroTuplesRanking(vector<T> &x)
{
    vector<pair<int, T>> solutions;
    for(int i = 0; i < x.size(); i++)
    {
        if (x[i] != 0)
        {
            solutions.push_back(make_pair(i+1, x[i]));
        }
    }
    sort(solutions.begin(), solutions.end(), [](const pair<int, double> &a, const pair<int, double> &b)
    {
        return a.second > b.second;
    });
    return solutions;
}

template <typename T>
inline T extractGRBVar(GRBVar &var)
{
    double val = var.get(GRB_DoubleAttr_X);
    char vtype = var.get(GRB_CharAttr_VType);
    if (vtype == GRB_INTEGER)
    {
        return static_cast<T>(std::round(val));
    }
    else if (vtype == GRB_BINARY)
    {
        return static_cast<T>(std::round(std::clamp(val, 0.0, 1.0)));
    }
    else
    {
        return static_cast<T>(val);
    }
}

inline bool applyTimeBudget(GRBModel &model, SolveOptions &options)
{
    gpro.stop("effectiveRuntime");
    double remainingTime = options.timeBudget - (gpro.getTime("effectiveRuntime") / 1000);
    gpro.clock("effectiveRuntime");
    if (remainingTime <= 0)
    {
        cout << "Time Budget Exceeded" << endl;
        return false;
    }
    model.set(GRB_DoubleParam_TimeLimit, remainingTime);
    return true;
}

inline bool optimizeRelaxed(GRBModel &model)
{
    model.set(GRB_IntParam_OutputFlag, 0);
    model.set(GRB_IntParam_StartNodeLimit, 0);
    model.optimize();
    int status = model.get(GRB_IntAttr_Status);
    if (status == GRB_TIME_LIMIT)
    {
        if (model.get(GRB_IntAttr_SolCount) == 0)
        {
            cout << "Relaxed optimization: Time limit with no solution" << endl;
            return false;
        }
        cout << "Relaxed optimization: Time limit, using best incumbent" << endl;
    }
    else if (status != GRB_OPTIMAL && status != GRB_SUBOPTIMAL)
    {
        cout << "Relaxed optimization ended with status: " << status << endl;
        return false;
    }
    return true;
}

template <typename T>
inline void extractXSolution(GRBModel &model, std::vector<T> &x,
                              int numVars, bool reduced,
                              const vector<int> &reducedIds)
{
    std::fill(x.begin(), x.end(), T(0));
    if (reduced)
    {
        for (size_t i = 0; i < reducedIds.size(); i++)
        {
            int id = reducedIds[i] - 1;
            GRBVar var = model.getVarByName("xx[" + std::to_string(i) + "]");
            x[id] = extractGRBVar<T>(var);
        }
    }
    else
    {
        for (int i = 0; i < numVars; i++)
        {
            GRBVar var = model.getVarByName("xx[" + std::to_string(i) + "]");
            x[i] = extractGRBVar<T>(var);
        }
    }
}

template <typename T>
inline bool relaxVarianceConstraint(GRBModel &model, std::vector<T> &x,
                                    int numVars, bool reduced,
                                    const vector<int> &reducedIds,
                                    SolveOptions &options)
{
    try
    {
        GRBConstr *constrs = model.getConstrs();
        int numConstrs = model.get(GRB_IntAttr_NumConstrs);
        double *rhspen = new double[numConstrs];
        bool foundVarianceConstr = false;

        cout << "IN RELAX SOLVE" << numConstrs << endl;
        for (int i = 0; i < numConstrs; i++)
        {
            string name = constrs[i].get(GRB_StringAttr_ConstrName);
            if (name.find("mad_control") != string::npos ||
                name.find("variance_control") != string::npos ||
                name.find("prob_control") != string::npos ||
                name == "cvar")
            {
                rhspen[i] = 1.0;
                foundVarianceConstr = true;
                cout << "Relaxing constraint: " << name << endl;
            }
            else
                rhspen[i] = GRB_INFINITY;
        }

        double feasObj = model.feasRelax(0, true, 0, nullptr, nullptr, nullptr,
                                         numConstrs, constrs, rhspen);
        cout << "Variance relaxation amount: " << feasObj << endl;
        delete[] rhspen;
        delete[] constrs;

        if (options.timeBudgeted && !applyTimeBudget(model, options)) return false;
        if (!optimizeRelaxed(model)) return false;
        extractXSolution(model, x, numVars, reduced, reducedIds);
        return true;
    }
    catch (GRBException &e)
    {
        cout << "FeasRelax error: " << e.getErrorCode() << " - " << e.getMessage() << endl;
        return false;
    }
}

// solve() -> model.optimize(); then translate the x[i] GRBVar vector into a x vector so you can manipulate it.
template <typename T>
inline void Solver::solve(GRBModel &model, std::vector<T> &x, SolveOptions &options)
{
    deb(options.timeBudgeted, options.timeBudget, options.reduced);
    gpro.stop("effectiveRuntime");
    if (options.timeBudgeted)
    {
        double remainingTime = options.timeBudget - (gpro.getTime("effectiveRuntime") / 1000);
        // had a bug here --> if remainingTime is negative GRBExcept
        if (remainingTime <= 0)
        {
            std::cout << "Time Budget Exceeded before optimization" << std::endl;
            x.clear();
            return;
        }
        deb(remainingTime);
        model.set(GRB_DoubleParam_TimeLimit, remainingTime);
    }
    gpro.clock("effectiveRuntime");
    gpro.clock("gurobi");
    int numConstraints = model.get(GRB_IntAttr_NumConstrs);
    int numVars = model.get(GRB_IntAttr_NumVars);
    deb(numConstraints, numVars);
    model.set(GRB_IntParam_Presolve, 0);

    model.optimize();
    int status = model.get(GRB_IntAttr_Status);

    if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE || model.get(GRB_IntAttr_Status) == GRB_INF_OR_UNBD)
    {
        cout << "Initial model is infeasible or unbounded" << std::endl;

        if (options.enableRelaxation &&
            relaxVarianceConstraint(model, x, this->NTuples, options.reduced, options.reducedIds, options))
        {
            options.foundSolutionUsingGurobiRelax = true;
            cout << "Found solution with relaxed variance constraint" << std::endl;
            gpro.stop("gurobi");
            return;
        }else
        {
            x.clear();
            return;
        }
    }
    else if (status == GRB_TIME_LIMIT)
    {
        std::cout << "Time Limit Expired" << std::endl;
        int solCount = model.get(GRB_IntAttr_SolCount);
        if (solCount > 0)
        {
            std::cout << "Returning best incumbent solution found (" << solCount << " solutions found)" << std::endl;
            // Continue to extract the best solution below
        }
        else
        {
            std::cout << "No feasible solution found within time limit" << std::endl;
            x.clear();
            return;
        }
    }
    else
    {
        std::cout << "Model optimization ended with status: " << model.get(GRB_IntAttr_Status) << std::endl;
    }

    cout << "Solution Found" << endl;
    extractXSolution(model, x, this->NTuples, options.reduced, options.reducedIds);
    gpro.stop("gurobi");
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

inline double calculateE(vector<double> &x)
{
    double E = 0.0;
    for (int i = 0; i < x.size(); i++)
    {
        E += x[i];
    }

    return E;
}

inline void findUnion(vector<double> &solLP1, vector<double> &solLP2, vector<int> &reducedIds)
{
    for (int i = 0; i < solLP1.size(); i++)
    {
        if (solLP1[i] > 0 || solLP2[i] > 0)
        {
            reducedIds.push_back(i + 1);
        }
    }
}

template <typename T>
inline bool Solver::isCurrentBetterThanBest(vector<double> r, double W_q, SolutionMetadata<T> &bestSol, int sense)
{
    bool isFeas = isFeasible(r);
    if (isFeas && bestSol.isFeasible)
    {
        if (sense == GRB_MINIMIZE)
        {
            if (W_q < bestSol.w)
            {
                return true;
            }
        }
        else
        {
            if (W_q > bestSol.w)
            {
                return true;
            }
        }
    }
    else if (isFeas && !bestSol.isFeasible)
    {
        return true;
    }
    else if (!isFeas && !bestSol.isFeasible)
    {
        double worstRkSol = 1e7;
        double worstRkBestSol = 1e7;
        if (bestSol.bestRk.empty())
        {
            return true;
        }
        else
        {
            for (int i = 0; i < r.size(); i++)
            {
                if (r[i] < worstRkSol)
                {
                    worstRkSol = r[i];
                }
                if (bestSol.bestRk[i] < worstRkBestSol)
                {
                    worstRkBestSol = bestSol.bestRk[i];
                }
            }
            if (worstRkSol > worstRkBestSol)
            {
                return true;
            }
            else if (worstRkSol == worstRkBestSol)
            {
                if (sense == GRB_MINIMIZE)
                {
                    if (W_q < bestSol.w) // abusing maximization
                    {
                        return true;
                    }
                }
                else
                {
                    if (W_q > bestSol.w) // abusing maximization
                    {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}