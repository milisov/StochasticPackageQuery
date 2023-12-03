#ifndef TAYLOR_HPP
#define TAYLOR_HPP

#include <memory>
#include <map>
#include <string>
#include <vector>

#include "util/udeclare.hpp"
#include "spq/spq.hpp"
#include "core/stat.hpp"

using std::shared_ptr;
using std::unique_ptr;
using std::map;
using std::string;
using std::vector;

class Taylor{
private:
    double sqn, objValue, objNorm, minVio, bestObj;
    unique_ptr<Stat> stat;
    shared_ptr<StochasticPackageQuery> spq;
    vector<long long> ids;
    vector<vector<double>> detCons, stoXs, stoMeans;
    vector<double> detNorms, detXs, obj;
    map<string, Option> options;
    SolIndType sol, bestSol;
private:
    void solve(SolIndType& nextSol);
    void update(const SolIndType& step);
    void update(const double& vio, const double& objValue, const SolIndType& sol);
public:
    Taylor(shared_ptr<StochasticPackageQuery> spq, const vector<long long>& ids={}, const map<string, Option>& options={});
    void solve();
    SolType getSol() const;
};

#endif