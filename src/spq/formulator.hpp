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
    bool SAA = false;
    bool RS = false;
    bool objCons = false; // if true, we formulate the objective constraint
    bool reduced = false;
    bool reducedScenarios = false;
    vector< pair<int, double>> posActiveness;
    vector< pair<int, double>> negActiveness;
    bool computeActiveness = false;
    std::vector<int> reducedIds;
    std::vector<std::vector< pair<int, double>>> innerConstraints;
    std::map<std::string, double> cntoptions; 
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
    std::vector<int> cbasis; //warmstart
};

void setDecisionVarOptions(DecisionVarOptions &options, double lb, double ub, double obj, GrbVarType GrbVarType);

/**
 * Abstract base class for all formulators
 */
class Formulator{
protected:
GRBEnv env = GRBEnv();
PgManager pg;
string DB_optim;
string DB_valid;
int NTuples;
int cntScenarios;
//needed for partition/summarize

GRBVar addDecisionVar(GRBModel &model, DecisionVarOptions &options);

public:
    shared_ptr<StochasticPackageQuery> spq;
    std::vector<int> shuffler;
    std::vector<std::vector<pair<int, double>>> partitions;
    double fetchRuntime = 0.0;
    Data& data; 
    Formulator();
    Formulator(shared_ptr<StochasticPackageQuery> spq);

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
    void reshuffleShuffler(std::vector<int> &shuffler);
    void populateShuffler(std::vector<int> &v);
    void partition(int Z, std::vector<pair<int, double>> &innerConstraints, std::vector<int> &shuffler, std::vector<std::vector<pair<int, double>>> &partitions);
    std::vector<std::vector<double>> summarize(FormulateOptions &formOptions,
                                                         std::shared_ptr<ProbConstraint> probCon,
                                                         std::shared_ptr<AttrConstraint> attrCon,
                                                         int conOrder);
};