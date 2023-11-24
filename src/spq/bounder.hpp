#ifndef BOUNDER_HPP
#define BOUNDER_HPP

#include <memory>
#include <vector>

#include "spq/spq.hpp"
#include "core/kde.hpp"

using std::shared_ptr;
using std::unique_ptr;
using std::vector;

class Bounder{
private:
    static const double hardEps, hardLimit;
    shared_ptr<StochasticPackageQuery> spq;
    size_t N;
    double E;
    vector<unique_ptr<KDE>> kdes;
    vector<double> steps;
    vector<bool> reverses;
public:
    Bounder(shared_ptr<StochasticPackageQuery> spq, const size_t& N, const double& E);
};

#endif