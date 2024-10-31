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
#include "CurveFit.hpp"

using namespace std;

class History
{
public:
    std::vector<int> xBest;
    double wBest;
    double epsilonBest;
    bool foundFeasible = false;

    //solution vector x, all alpha_k, indicator if it's feasible or not, objective value Wq
    std::vector<std::tuple<std::vector<int>,std::vector<double>,bool,double>>bestSolMetadata;
    //alpha_k, r_k pairs for each probCons K
    std::vector<std::vector<pair<double, double>>> curveFitMetadata;

    History(){};
};

class SolutionMetadata
{
public:
    std::vector<int>x;
    double w;
    double epsilon;
    bool isFeasible;
    double Runtime;
    int binarySearchSteps = -1;
    int Z = 0;
    SolutionMetadata(const std::vector<int>& xInit, double wInit, double epsilonInit, bool isFeasibleInit): x(xInit), w(wInit), epsilon(epsilonInit), isFeasible(isFeasibleInit) {}
    SolutionMetadata(){}
};


class BinarySearchMetadata
{
public: 
    double low;
    double high;
    double alpha;
    BinarySearchMetadata(double lowInit, double highInit, double alphaInit): low(lowInit), high(highInit), alpha(alphaInit) {}
};

class SummarySearch
{
private:
    std::vector<std::vector<pair<double, double>>> history;
    PgManager pg;
    string fitFunction = "atan";
    std::vector<std::vector<double>> innerConstraints;
    std::vector<GRBVar> yy;
    std::vector<GRBGenConstr> genCon;
    std::vector<GRBConstr> sumyCon;
    double W_q;
    double W0;
    double epsilonQ;

public:
    History H;
    GRBEnv env = GRBEnv();
    GRBModel modelILP = GRBModel(env);
    std::unique_ptr<GRBVar[]> xxILP;
    

    GRBModel modelLP1 = GRBModel(env);
    std::unique_ptr<GRBVar[]> xxLP1;

    GRBModel modelLP2 = GRBModel(env);
    std::unique_ptr<GRBVar[]> xxLP2;

    CurveFitter fitter;
    //std::vector<int> x;
    double epsilon;
    int cntScenarios;
    int M;
    int M_hat;
    int NTuples;
    string DB_optim;
    string DB_valid;
    shared_ptr<StochasticPackageQuery> spq;
    string SPQPath;
    std::vector<int> shuffler;
    int probConstCnt = 0;
    // r vector contains the surplus values for each of the constraints in iteration q
    std::vector<double> r;
    SummarySearch(int M = 1e4, int M_hat = 1e6, shared_ptr<StochasticPackageQuery> spq = nullptr, double epsilon = 1e-5);
    template<typename T>
    void validate(std::vector<T> &x, shared_ptr<StochasticPackageQuery> spq, int M_hat, string DB_optim);

    template<typename T>
    double calculateObj(std::vector<T> &x, std::vector<int> &selectIds, shared_ptr<Objective> obj, string DB_optim);

    template<typename T>
    double calculateCntObj(std::vector<T> &x, std::vector<int> &selectIds, shared_ptr<CountObjective> cntObj, string DB_optim);
    
    template<typename T>
    double calculateSumObj(std::vector<T> &x, std::vector<int> &selectIds, shared_ptr<AttrObjective> attrObj, string DB_optim);

    template<typename T>
    double calculateExpSumObj(std::vector<T> &x, std::vector<int> &selectIds, shared_ptr<AttrObjective> attrObj, string DB_optim);

    std::vector<std::vector<double>> summarize(std::vector<int> &x, int M, int Z, double alpha, string DB_optim, shared_ptr<ProbConstraint> probCon, shared_ptr<AttrConstraint> attrCon, int conOrder,
    vector<int>&reducedIds, bool reduced);
    void reshuffleShuffler(std::vector<int> &shuffler);
    void populateShuffler(vector<int> &v);

    void formulateSAA(GRBModel &model, std::vector<std::vector<std::vector<double>>> &summaries, int q, GRBVar *xx, vector<int>&reducedIds, bool reduced, map<string, Option>& cntoptions);
    void formulateLP(GRBModel &model, GRBVar *xx, bool reduced, map<string, Option>& cntoptions);
    void solveLP2(GRBModel &model, vector<double>&sol, GRBVar *xx, double ub,  vector<int>dummyVect);
    void formCountCons(GRBModel &model, shared_ptr<Constraint> cons, GRBVar *xx, vector<int>&reducedIds, bool reduced, map<string, Option>& options);
    void formSumCons(GRBModel &model, shared_ptr<Constraint> cons, GRBVar *xx, vector<int>&reducedIds, bool reduced);
    void formProbCons(GRBModel &model, shared_ptr<Constraint> cons, GRBVar *xx, std::vector<std::vector<std::vector<double>>> &summaries, 
    int &probConOrder, std::vector<GRBVar> &yy, std::vector<GRBGenConstr> &genCon, 
    std::vector<GRBConstr> &sumyCon, vector<int>&reducedIds, bool reduced);
    void formExpCons(GRBModel &model, shared_ptr<Constraint> cons, GRBVar *xx, vector<int>&reducedIds, bool reduced);
    void formSumObj(GRBModel &model, shared_ptr<Objective> obj, GRBVar *xx, vector<int>&reducedIds, bool reduced);
    void formExpSumObj(GRBModel &model, shared_ptr<Objective> obj, GRBVar *xx, vector<int>&reducedIds, bool reduced);
    void formCntObj(GRBModel &model, shared_ptr<Objective> obj, GRBVar *xx, vector<int>&reducedIds, bool reduced);

    void formCVaR(GRBModel &model, shared_ptr<Constraint> cons, GRBVar *xx);

    SolutionMetadata CSASolve(std::vector<int> &x, int M, int Z, vector<int>&reducedIds, bool reduced, map<string, Option>& cntoptions);
    SolutionMetadata CSASolveBinSearch(std::vector<int> &x, int M, int Z, vector<int>&reducedIds, bool reduced, map<string, Option>& cntoptions);
    BinarySearchMetadata guessOptimalConservativenessBinarySearch(double low, double high, double rk);


    template<typename T>
    void solve(GRBModel &model, std::vector<T> &x, GRBVar *xx, vector<int>&reducedIds, bool reduced);


    

    SolutionMetadata summarySearch(shared_ptr<StochasticPackageQuery> spq, int M_hat, int M, int m, int z, vector<int>&reducedIds, bool reduced, map<string, Option>& cntoptions, map<string, Option>& curveFitOptions);
    SolutionMetadata stochDualReducer(shared_ptr<StochasticPackageQuery> spq, int qSz, map<string, Option>& cntoptions, map<string, Option>& curveFitOptions);

    void guessOptimalConservativeness(std::vector<std::vector<pair<double, double>>> &history, std::vector<double> &alpha);
    void Best(shared_ptr<Objective> obj, std::vector<int> &x, History &H);
};