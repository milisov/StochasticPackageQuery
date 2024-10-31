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


class Solution
{
public:
    std::vector<int>x;
    bool isFeasible;
    double W_q;
    double Runtime;
    Solution(const std::vector<int>& xInit, double W_qinit, bool isFeasibleInit): x(xInit), W_q(W_qinit), isFeasible(isFeasibleInit) {}
};

class Naive
{
private:
    PgManager pg;
    double W_q;
    std::vector<std::vector<double>> innerConstraints;
    std::vector<GRBVar> yy;
    std::vector<GRBGenConstr> genCon;
    std::vector<GRBConstr> sumyCon;

public:
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);
    std::unique_ptr<GRBVar[]> xx;
    vector<double>r;
    int cntScenarios;
    int M;
    int M_hat;
    int NTuples;
    string DB_optim;
    string DB_valid;
    shared_ptr<StochasticPackageQuery> spq;
    int probConstCnt = 0;
    Naive(int M = 1e4, int M_hat = 1e6, shared_ptr<StochasticPackageQuery> spq = nullptr);
    void validate(std::vector<int> &x, shared_ptr<StochasticPackageQuery> spq, int M_hat, string DB_optim);
    void solve(std::vector<int> &x, GRBVar *xx);
    Solution solveNaive(shared_ptr<StochasticPackageQuery> spq, int m, int M, int M_hat);

    void formulateSAA();
    void formCountCons(shared_ptr<Constraint> cons, GRBVar *xx);
    void formSumCons(shared_ptr<Constraint> cons, GRBVar *xx);
    void formProbCons(shared_ptr<Constraint> cons, GRBVar *xx);
    void formExpCons(shared_ptr<Constraint> cons, GRBVar *xx);
    void formSumObj(shared_ptr<Objective> obj, GRBVar *xx);
    void formExpSumObj(shared_ptr<Objective> obj, GRBVar *xx);
    void formCntObj(shared_ptr<Objective> obj, GRBVar *xx);

    double calculateObj(std::vector<int> &x, std::vector<int> &selectIds, shared_ptr<Objective> obj, string DB_optim);
    double calculateCntObj(std::vector<int> &x, std::vector<int> &selectIds, shared_ptr<CountObjective> cntObj, string DB_optim);
    double calculateSumObj(std::vector<int> &x, std::vector<int> &selectIds, shared_ptr<AttrObjective> attrObj, string DB_optim);
    double calculateExpSumObj(std::vector<int> &x, std::vector<int> &selectIds, shared_ptr<AttrObjective> attrObj, string DB_optim);
};