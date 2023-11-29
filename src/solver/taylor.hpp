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
    double sqn;
    unique_ptr<Stat> stat;
    shared_ptr<StochasticPackageQuery> spq;
    vector<long long> ids;
    vector<vector<double>> detCons, stoXs;
    vector<double> detNorms, detXs, obj;
    map<string, Option> options;
    map<int, double> sol;
public:
    void solve(map<int, double>& nextSol) const;
public:
    Taylor(shared_ptr<StochasticPackageQuery> spq, const vector<long long>& ids={}, const map<string, Option>& options={});
};

#endif