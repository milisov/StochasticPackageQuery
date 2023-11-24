#include <pcg_random.hpp>
#include <random>
#include <cmath>
#include <unordered_map>
#include <utility>
#include <iostream>
#include <fmt/core.h>
#include <boost/math/special_functions/erf.hpp>

#include "bounder.hpp"
#include "stat.hpp"
#include "util/udebug.hpp"
#include "util/uconfig.hpp"
#include "util/uio.hpp"

using std::unordered_map;
using std::pair;
using std::dynamic_pointer_cast;
using std::cerr;

const double Bounder::hardEps = Config::getInstance()->pt.get<double>("parameters.hard_eps");
const double Bounder::hardLimit = -log(Bounder::hardEps/(1-Bounder::hardEps));

Bounder::Bounder(shared_ptr<StochasticPackageQuery> spq, const size_t& N, const double& E): spq(spq), N(N), E(E){
    PgManager pg;
    Stat stat;
    int stoSize = spq->countStochastic();
    if (stoSize){
        pcg32 gen (Config::getInstance()->seed());
        long long n = pg.getTableSize(spq->tableName);
        size_t iE;
        double whole;
        double fE = modf(E, &whole);
        iE = static_cast<size_t>(whole)+1;
        vector<vector<long long>> choices (N, vector<long long>(iE));
        std::uniform_int_distribution<long long> dist (0LL, n-1);
        unordered_map<long long, pair<double, double>> storages;
        storages.reserve(N*iE);
        for (size_t i = 0; i < N; ++i){
            for (size_t j = 0; j < iE; ++j){
                choices[i][j] = dist(gen);
                storages.try_emplace(choices[i][j]);
            }
        }
        for (auto con : spq->cons){
            shared_ptr<ProbConstraint> probCon = dynamic_pointer_cast<ProbConstraint>(con);
            if (probCon){
                shared_ptr<AttrConstraint> attrCon = dynamic_pointer_cast<AttrConstraint>(con);
                if (attrCon){
                    stat.getStochasticMeanVar(spq->tableName, attrCon->attr, storages);
                } else{
                    cerr << fmt::format("Stochastic constraint '{}' needs to have an attribute\n", static_cast<string>(*con));
                    exit(1);
                }
                double c = 0;
                shared_ptr<VarConstraint> varCon = dynamic_pointer_cast<VarConstraint>(con);
                if (varCon){
                    double p = spq->getValue(varCon->p);
                    if (varCon->vsign == Ineq::gteq) p = 1-p;
                    c = sqrt(2)*boost::math::erf_inv(2*p-1);
                    bool reverse = false;
                    if (varCon->vsign == varCon->psign) reverse = true;
                    reverses.push_back(reverse);
                }
                double step = 0;
                vector<double> G (N);
                for (size_t i = 0; i < N; ++i){
                    double mean = 0, var = 0;
                    for (size_t j = 0; j < iE-1; ++j){
                        const auto& pd = storages[choices[i][j]];
                        mean += pd.first;
                        var += pd.second;
                    }
                    const auto& pd = storages[choices[i].back()];
                    mean += pd.first*fE;
                    var += pd.second*fE;
                    G[i] = mean + sqrt(var)*c;
                    step += var;
                }
                kdes.push_back(std::make_unique<KDE>(G, true));
                step = sqrt(step/N)*abs(c);
                steps.push_back(step);
            }
        }
    }
}