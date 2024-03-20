#include "starter.hpp"

using std::make_unique;

map<string, Option> starterOptions = {
};

Starter::Starter(shared_ptr<StochasticPackageQuery> spq, const vector<long long>& ids, const map<string, Option>& options): spq(spq), options(options){
    stat = make_unique<Stat>();
    for (const auto& p : starterOptions){
        if (!options.count(p.first)) this->options[p.first] = p.second;
    }
    
}