#ifndef STARTER_HPP
#define STARTER_HPP

#include <vector>
#include <map>
#include <string>

#include "spq/spq.hpp"
#include "util/udeclare.hpp"
#include "core/stat.hpp"

using std::shared_ptr;
using std::unique_ptr;
using std::string;
using std::vector;

class Starter{
private:
    unique_ptr<Stat> stat;
    shared_ptr<StochasticPackageQuery> spq;
    map<string, Option> options;
public:
    SolIndType sol;
    Starter(shared_ptr<StochasticPackageQuery> spq, const vector<long long>& ids={}, const map<string, Option>& options={});
    
};

#endif