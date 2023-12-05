#ifndef TAYLOR_HPP
#define TAYLOR_HPP

#include <memory>
#include <map>
#include <unordered_map>
#include <string>
#include <vector>

#include "util/udebug.hpp"
#include "util/udeclare.hpp"
#include "spq/spq.hpp"
#include "core/stat.hpp"

using std::shared_ptr;
using std::unique_ptr;
using std::map;
using std::unordered_map;
using std::string;
using std::vector;

class Taylor{
private:
    ObjectiveSense objSense;
    int nMaxIters, iter;
    bool isSoftDetCon, isDependencyVar;
    double sqn, maxSolutionSize, objValue, minVio, bestObj;
    unique_ptr<Stat> stat;
    shared_ptr<StochasticPackageQuery> spq;
    vector<long long> ids;
    vector<int> vStart, cStart;
    vector<vector<double>> detCons, stoMeans, stoXs;
    vector<double> detNorms, obj, detXs, pStart, dStart, preViolations;
    map<string, Option> options;
    SolIndType sol, bestSol;
    unordered_map<size_t, vector<SolIndType>> hashedSols;
    INIT(pro);
private:
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