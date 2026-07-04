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
#include "util/data.hpp"
#pragma once

void partition(int Z, std::vector<pair<double, double>> &realizationScore, std::vector<int> &shuffler, std::vector<std::vector<pair<double, double>>> &partitions);
int countProbConst(shared_ptr<StochasticPackageQuery> spq);
bool sortbysecASC(const pair<double, double> &a, const pair<double, double> &b);
bool sortbysecDESC(const pair<double, double> &a, const pair<double, double> &b);

// input: sorted partition, boolean that determines if we need min/max summary, alpha value to consider
// G(alpha) = |alpha*Пz| elements to consider in order to find min/max
template <typename T1, typename T2>
double calculateSummary(vector<double> &scenarios,std::vector<pair<T1, T2>> &innerCons, bool maxS, double alpha)
{
    //deb("CALCULATE SUMMARY", alpha, innerCons.size());
    int G_alpha = (int)ceil(alpha * (double)innerCons.size());
    int id = innerCons[0].first;
    double summary = scenarios[id];
    for (int i = 1; i < G_alpha; i++)
    {
        if (maxS)
        {
            id = innerCons[i].first;
            if (scenarios[id] > summary)
            {
                summary = scenarios[id];
            }
        }
        else
        {
            id = innerCons[i].first;
            if (scenarios[id] < summary)
            {
                summary = scenarios[id];
            }
        }
    }
    return summary;
}

//this function gets scenarios, for the row and the partition
template <typename T1, typename T2>
double calculateSummarySmooth(vector<double> &scenarios,std::vector<pair<T1, T2>> &innerCons, bool maxS, double alpha, double partitionMax)
{
    int topKplus1 = (int)ceil(alpha * (double)innerCons.size());
    int topK = topKplus1 - 1;

    double alpha0 = topK / (double) innerCons.size();
    double alpha1 = topKplus1 / (double) innerCons.size();

    // if(innerCons.size() < 5)
    // {
    //     deb(alpha, alpha0, alpha1, topKplus1, topK, innerCons);
    // }

    double summary1 = 0.0;
    double summary0;
    if(topK == 0)
    {
        summary0 = partitionMax + 2*abs(partitionMax);
    }else
    {
        int id = innerCons[0].first;
        summary0 = scenarios[id];
    }


    for (int i = 1; i < topK; i++)
    {
        if (maxS)
        {
            int id = innerCons[i].first;
            if (scenarios[id] > summary0)
            {
                summary0 = scenarios[id];
            }
        }
        else
        {
            int id = innerCons[i].first;
            if (scenarios[id] < summary0)
            {
                summary0 = scenarios[id];
            }
        }
    }
    int id = innerCons[topK].first; 
    double scenarioValue = scenarios[id];
    summary1 = summary0;
    if(maxS)
    {
        if(scenarioValue > summary1 || topK == 0)
        {
            summary1 = scenarioValue;
        }
    }else
    {
        if(scenarioValue < summary1 || topK == 0)
        {
            summary1 = scenarioValue;
        }
    }
    double step = (alpha - alpha0) / (alpha1 - alpha0);
    // if(innerCons.size() < 5)
    // {
    //     deb(summary0, summary1, step);
    // }
    double summary = summary0 + (summary1 - summary0) * step;
    return summary;
}




// function that initializes given vector
template <typename T>
void initializeVectorForm(std::vector<T> &v, int sz, T init)
{
    for (int i = 0; i < sz; i++)
    {
        v.push_back(init);
    }
}


enum class GrbVarType {
    Binary,       // GRB_BINARY
    Integer,      // GRB_INTEGER
    Continuous    // GRB_CONTINUOUS (real)
};


struct DecisionVarOptions{
    double lb;
    double ub; 
    double obj; 
    GrbVarType varType;
    string name; 
};

struct FormulateOptions{
    bool RS = false;
    bool activenessPartition = false;
    bool objCons = false; // if true, we formulate the objective constraint
    bool reduced = false;
    map<string, bool> ablate;
    bool reducedScenarios = false;
    bool partitionMostActive = false;
    bool partitionSpreadActiveness = false;
    bool includeObjectiveFunction = true;
    bool finalStage = false; // if true, in the summary search in stage 3, we use the given objective instead of minmax
    bool indicators = false;

    bool lcvarFormulation = false;
    bool clearLCVaRCache = false;
    bool changeVaRBound = false;
    vector<double> meanPerVaR;
    vector<double> variancePerVaR;
    vector<double> minInnerConstPerVaR;
    bool varianceControl = false;
    std::vector<std::set<int>> mostActiveScenariosPerConstraint;  // Fixed active scenarios per constraint (if set, used instead of computing)
    std::vector<std::set<int>> reducedScenariosPerConstraint;
    vector<double> kMADValues;


    vector<vector<pair<int, double>>> activeness;
    vector<vector<pair<int, double>>> activenessNoAbsValue;
    std::vector<int> reducedIds;
    std::vector<std::vector< pair<int, double>>> innerConstraints;
    std::vector<vector<vector<double>>> partitionMaximums;
    int M; 
    int Z;
    int Zinit;
    double low;
    double high;
    int qSz;
    double objValue;
    vector<double>alpha;
    int iteration = 0;
    DecisionVarOptions decisionVarOptions;
    std::vector<int> vbasis; //warmstart 
    std::vector<int> cbasis;
    double p;
    unsigned int randomSeed = 42;
};

void setDecisionVarOptions(DecisionVarOptions &options, double lb, double ub, double obj, GrbVarType GrbVarType);

/**
 * Abstract base class for all formulators
 */
class Formulator{
protected:
    GRBEnv env = GRBEnv(true);
    PgManager pg;
    string DB_optim;
    string DB_valid;
    int NTuples;
    int cntScenarios;
    static std::map<std::string, std::vector<double>> lcvarCoeffsCache;

    GRBVar addDecisionVar(GRBModel &model, DecisionVarOptions &options);
    std::vector<double> computeLCVaRCoeffs(std::shared_ptr<AttrConstraint> attrCon, double p, FormulateOptions &formOptions);
    
    public:
    double fetchRuntime = 0.0;
    shared_ptr<StochasticPackageQuery> spq;
    //needed for partition/summarize
    std::vector<int> shuffler;
    std::vector<vector<vector<pair<int, double>>>> partitions;
    Data& data; 
    Formulator();
    Formulator(shared_ptr<StochasticPackageQuery> spq);


    virtual GRBModel& getModel() {
        throw std::logic_error("Formulator has no model");
    }
    
    //Build and return a Gurobi model ready for optimization
    virtual GRBModel formulate(shared_ptr<StochasticPackageQuery> spq, FormulateOptions& FormOptions) = 0;
    
    //standard constraints that have same implementation for all solvers
    void formCountCons(GRBModel &model,shared_ptr<Constraint> cons, GRBVar *xx, FormulateOptions& options);
    void formSumCons(GRBModel &model,shared_ptr<Constraint> cons, GRBVar *xx,  FormulateOptions& options);
    void formExpCons(GRBModel &model,shared_ptr<Constraint> cons, GRBVar *xx, FormulateOptions& options);
    void formLCVaR(GRBModel &model,shared_ptr<Constraint> cons, GRBVar *xx, FormulateOptions& options);
    void formSumObj(GRBModel &model,shared_ptr<Objective> obj, GRBVar *xx, FormulateOptions& options);
    void formExpSumObj(GRBModel &model,shared_ptr<Objective> obj, GRBVar *xx, FormulateOptions& options);
    void formCntObj(GRBModel &model,shared_ptr<Objective> obj, GRBVar *xx, FormulateOptions& options);    
    void partition(int Z, FormulateOptions &formOptions, std::vector<pair<int, double>> &innerConstraints, std::vector<int> &shuffler, int conOrder, std::shared_ptr<AttrConstraint> attrCon);

    std::set<int> getReducedScenariosFromActiveness(vector<pair<int, double>> &activeness, int Z, double p);
    std::set<int> getMostActiveScenarios(vector<pair<int, double>> &activeness, int Z, double p);
    std::vector<int> getMostActiveScenariosPerPartition(int Z, std::vector<std::pair<int, double>> &sortedActiveness);
    void stage3Partition(FormulateOptions& options,double p, int conOrder, std::shared_ptr<AttrConstraint> attrCon);
    void activenessPartition(FormulateOptions& options, vector<pair<int, double>> &activeness, std::vector<pair<int, double>> &innerConstraints, double p);
    void createPartitions(shared_ptr<StochasticPackageQuery> spq, FormulateOptions& options);
    void reshuffleShuffler(std::vector<int> &shuffler);
    void populateShuffler(std::vector<int> &v);
    std::vector<std::vector<double>> summarize(FormulateOptions &formOptions,
                                                         std::shared_ptr<ProbConstraint> probCon,
                                                         std::shared_ptr<AttrConstraint> attrCon,
                                                         int conOrder);
};