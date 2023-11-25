#ifndef BOUNDER_HPP
#define BOUNDER_HPP

#include <memory>
#include <vector>
#include <map>

#include "stat.hpp"
#include "spq.hpp"
#include "core/kde.hpp"
#include "util/udeclare.hpp"

using std::shared_ptr;
using std::unique_ptr;
using std::vector;

class Bounder{
private:
    static const double hardEps, hardLimit;
    shared_ptr<StochasticPackageQuery> spq;
    unique_ptr<Stat> stat;
    size_t N;
    double E;
    vector<unique_ptr<KDE>> stoKDEs, detKDEs;
    vector<double> steps;
    vector<bool> reverses;
    map<double, map<string, double>> varTables;
private:
    void setBound(const double& h, const Bound& bound, double value);
public:
    Bounder(shared_ptr<StochasticPackageQuery> spq, const size_t& N, const double& E);
    void generate(const vector<double>& hards);
    void set(const double& h);
};

#endif