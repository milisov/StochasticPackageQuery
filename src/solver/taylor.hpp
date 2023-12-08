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
    static const double adjustCoef;
    unique_ptr<UniqueIndexer> idx;
    unique_ptr<Stat> stat;
    shared_ptr<StochasticPackageQuery> spq;
    ObjectiveSense objSense;
    int nMaxIters, nCores;
    bool isSoftDetCon, isDependentVar;
    double sqn, gamma; 
    map<string, Option> options;
    // Persistent members_Deterministic
    vector<vector<double>> detCons;
    vector<double> detNorms, obj;
    // Persistent members_Stochastic
    vector<double> nsix;
    vector<vector<double>> stoMeans, stoVars;
    map<string, pair<double, double>> adjustments;
    // Persistent members_Stochastic_IsDependentVar
    // System-update
    INIT(pro);
    int iter;
    double maxSolutionSize, minVio, bestObj;
    SolIndType sol, bestSol;
    unordered_map<size_t, vector<SolIndType>> hashedSols;
    // System update_Gurobi
    vector<int> vStart, cStart;
    vector<double> pStart, dStart, preViolations;
    // Stochastic update
    vector<vector<double>> stoXs;
    // Stochastic update_IsDependentVar
    vector<double> varXs;
    vector<vector<double>> zeroPds;
    vector<set<size_t>> jvInds;
    // Deterministic update
    double objValue;
    vector<double> detXs;
private:
    void doAdjustment();
    void undoAdjustment();
    void solve(SolIndType& nextSol);
    void update(const SolIndType& step);
    void update(const double& vio, const double& objValue, const SolIndType& sol);
public:
    TaylorStatus status;
    Taylor(shared_ptr<StochasticPackageQuery> spq, const vector<long long>& ids={}, const map<string, Option>& options={});
    void solve();
    SolType getSol(const SolIndType& sol={}) const;
};

#endif