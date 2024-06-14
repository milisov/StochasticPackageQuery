#ifndef TAYLOR_HPP
#define TAYLOR_HPP

#include <memory>
#include <map>
#include <unordered_map>
#include <set>
#include <string>
#include <vector>
#include <utility>

#include "util/udebug.hpp"
#include "util/udeclare.hpp"
#include "util/unumeric.hpp"
#include "spq/spq.hpp"
#include "core/stat.hpp"
#include "core/optim.hpp"

using std::shared_ptr;
using std::unique_ptr;
using std::map;
using std::unordered_map;
using std::string;
using std::pair;
using std::set;
using std::vector;

class Taylor{
private:
    // Persistent members
    unique_ptr<RMSprop> optim;
    static const double adjustCoef;
    unique_ptr<UniqueIndexer> idx;
    unique_ptr<Stat> stat;
    shared_ptr<StochasticPackageQuery> spq;
    ObjectiveSense objSense;
    int nMaxIters, nCores;
    bool isSoftDetCon;
    map<string, Option> options;
    SolIndType sol, bestSol;
    // double sqn, gamma, preVio; 
    // Persistent members_Deterministic
    vector<vector<double>> detCons;
    vector<double> obj;
    // // Persistent members_Stochastic
    // vector<double> nsix;
    // vector<vector<double>> stoMeans, stoVars;
    // map<string, pair<double, double>> adjustments;
    // // System-update
    // INIT(pro);
    double bestObj;
    // double maxSolutionSize, minVio, bestObj;
    // unordered_map<size_t, vector<SolIndType>> hashedSols;
    // // System update_Gurobi
    // vector<int> vStart, cStart;
    // vector<double> pStart, dStart, preViolations;
    // Stochastic update
    vector<vector<double>> stoXs;
    // Deterministic update
    double objValue;
    vector<double> detXs;
private:
    // void doAdjustment();
    // void undoAdjustment();
    void solve(SolIndType& nextSol);
    void update(const SolIndType& step);
    // void update(const double& vio, const double& objValue, const SolIndType& sol);
public:
    TaylorStatus status;
    Taylor(shared_ptr<StochasticPackageQuery> spq, const vector<long long>& ids={}, const SolIndType& initSol={}, const map<string, Option>& options={});
    void solve();
    // SolType getSol(const SolIndType& sol={}) const;
};

#endif