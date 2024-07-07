#include <iostream>
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

using namespace std;

class SummarySearch
{
private:
    // vector<vecor<pair<double,double>>> history
    PgManager pg;
    string fitFunction = "atan";
    vector<vector<double>> innerConstraints;
    vector<GRBVar> yy;
	vector<GRBGenConstr> genCon;
    vector<GRBConstr> sumyCon;
public:
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);
    std::unique_ptr<GRBVar[]> xx;
    int M;
    // M_hat
    int M_hat;
    int NTuples;
    string DB_optim;
    string DB_valid;
    shared_ptr<StochasticPackageQuery> spq;
    string SPQPath;
    vector<int> shuffler;
    vector<int> x;
    int probConstCnt = 0;
    // r vector contains the surplus values for each of the constraints in iteration q
    vector<double> r;
    SummarySearch(int M = 1e4, int M_hat = 1e6, shared_ptr<StochasticPackageQuery> spq = nullptr);
    void validate(vector<int> &x, shared_ptr<StochasticPackageQuery> spq, int M_hat, string DB_optim);
    vector<vector<double>> summarize(vector<int> &x, int Z, double alpha, string DB_optim, shared_ptr<ProbConstraint> probCon, shared_ptr<AttrConstraint> attrCon, int conOrder);
    void formulateSAA(vector<vector<vector<double>>> &summaries, int q);
    void formCountCons(shared_ptr<Constraint> cons, GRBVar *xx);
    void formSumCons(shared_ptr<Constraint> cons, GRBVar *xx);
    void formProbCons(shared_ptr<Constraint> cons, GRBVar *xx, vector<vector<vector<double>>> &summaries, int &probConOrder, vector<GRBVar>&yy, vector<GRBGenConstr>&genCon, vector<GRBConstr>&sumyCon);
    void formExpCons(shared_ptr<Constraint> cons, GRBVar *xx);
    void formSumObj(shared_ptr<Objective> obj, GRBVar *xx);
    void formExpSumObj(shared_ptr<Objective> obj, GRBVar *xx);
    void formCntObj(shared_ptr<Objective> obj, GRBVar *xx);
    void reshuffleShuffler(vector<int> &shuffler);
};