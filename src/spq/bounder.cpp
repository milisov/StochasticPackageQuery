#include <pcg_random.hpp>
#include <random>
#include <cmath>
#include <unordered_map>
#include <utility>
#include <iostream>
#include <fmt/core.h>
#include <boost/variant.hpp>
#include <boost/math/special_functions/erf.hpp>

#include "bounder.hpp"
#include "util/udebug.hpp"
#include "util/uconfig.hpp"
#include "util/uio.hpp"
#include "util/unumeric.hpp"

using std::unordered_map;
using std::pair;
using std::dynamic_pointer_cast;
using std::cerr;

const double Bounder::hardEps = Config::getInstance()->pt.get<double>("parameters.hard_eps");
const double Bounder::hardLimit = -log(Bounder::hardEps/(1-Bounder::hardEps));

Bounder::Bounder(shared_ptr<StochasticPackageQuery> spq, const size_t& N, const double& E): spq(spq), N(N), E(E){
    PgManager pg;
    stat = std::make_unique<Stat>();
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
        unordered_map<long long, double> detrages;
        storages.reserve(N*iE);
        for (size_t i = 0; i < N; ++i){
            for (size_t j = 0; j < iE; ++j){
                choices[i][j] = dist(gen);
                storages.try_emplace(choices[i][j]);
                detrages.try_emplace(choices[i][j]);
            }
        }
        for (auto con : spq->cons){
            shared_ptr<BoundConstraint> boundCon = dynamic_pointer_cast<BoundConstraint>(con);
            if (boundCon){
                shared_ptr<AttrConstraint> attrCon = dynamic_pointer_cast<AttrConstraint>(con);
                if (attrCon){
                    stat->getDetAttrs(spq->tableName, attrCon->attr, detrages);
                    vector<double> G (N);
                    for (size_t i = 0; i < N; ++i){
                        double attr = 0;
                        for (size_t j = 0; j < iE-1; ++j) attr += detrages[choices[i][j]];
                        attr += detrages[choices[i].back()]*fE;
                        G[i] = attr;
                    }
                    detKDEs.push_back(std::make_unique<KDE>(G, true));
                }
            }
            shared_ptr<ProbConstraint> probCon = dynamic_pointer_cast<ProbConstraint>(con);
            if (probCon){
                shared_ptr<AttrConstraint> attrCon = dynamic_pointer_cast<AttrConstraint>(con);
                if (attrCon){
                    stat->getStoMeanVar(spq->tableName, attrCon->attr, storages);
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
                stoKDEs.push_back(std::make_unique<KDE>(G, true));
                step = sqrt(step/N)*abs(c);
                steps.push_back(step);
            }
        }
    }
}

void Bounder::setBound(const double& h, const Bound& bound, double value){
    if (bound.which() == 0){
        string var = boost::get<string>(bound);
        varTables[h][var] = value;
    }
}

void Bounder::generate(const vector<double>& hards){
    auto n = hards.size();
    vector<double> probs (n);
    for (size_t i = 0; i < n; ++i){
        varTables.try_emplace(hards[i]);
        probs[i] = sigmoid(-hards[i]);
    }
    int stoIndex = 0;
    int detIndex = 0;
    for (auto con : spq->cons){
		shared_ptr<BoundConstraint> boundCon = dynamic_pointer_cast<BoundConstraint>(con);
		if (boundCon){
		    shared_ptr<CountConstraint> countCon = dynamic_pointer_cast<CountConstraint>(con);
            if (countCon){
                for (const auto& h : hards){
                    setBound(h, countCon->lb, E);
                    setBound(h, countCon->ub, E);
                }
            }
            shared_ptr<AttrConstraint> attrCon = dynamic_pointer_cast<AttrConstraint>(con);
            if (attrCon){
                auto& kde = detKDEs[detIndex];
                if (boundCon->lb.which() == 0 && boundCon->ub.which() == 0){
                    vector<double> quantiles (2*n);
                    for (size_t i = 0; i < n; ++i){
                        quantiles[2*i] = (1-probs[i])/2;
                        quantiles[2*i+1] = (1+probs[i])/2;
                    }
                    kde->getQuantiles(quantiles);
                    for (size_t i = 0; i < n; ++i){
                        setBound(hards[i], boundCon->lb, quantiles[2*i]);
                        setBound(hards[i], boundCon->ub, quantiles[2*i+1]);
                    }
                } else if (boundCon->ub.which() == 0){
                    vector<double> quantiles = probs;
                    kde->getQuantiles(quantiles);
                    for (size_t i = 0; i < n; ++i) setBound(hards[i], boundCon->ub, quantiles[i]);
                } else if (boundCon->lb.which() == 0){
                    vector<double> quantiles (n);
                    for (size_t i = 0; i < n; ++i) quantiles[i] = 1-probs[i];
                    for (size_t i = 0; i < n; ++i) setBound(hards[i], boundCon->lb, quantiles[i]);
                }
                detIndex ++;
            }
        }
        shared_ptr<ProbConstraint> probCon = dynamic_pointer_cast<ProbConstraint>(con);
        if (probCon && probCon->v.which() == 0){
            string var = boost::get<string>(probCon->v);
            auto& kde = stoKDEs[stoIndex];
            vector<double> quantiles = probs;
            if (reverses[stoIndex]){
                for (size_t i = 0; i < n; ++i) quantiles[i] = 1-quantiles[i];
            }
            kde->getQuantiles(quantiles);
            for (size_t i = 0; i < n; ++i){
                double h = hards[i];
                double easyValue = kde->getMax() + steps[stoIndex]*(h-hardLimit);
                double hardValue = kde->getMin() + steps[stoIndex]*(hardLimit+h);
                if (reverses[stoIndex]) std::swap(easyValue, hardValue);
                if (h < -hardLimit) varTables[h][var] = easyValue;
                else if (h > hardLimit) varTables[h][var] = hardValue;
                else varTables[h][var] = quantiles[i];
            }
            stoIndex ++;
        }
    }
}

void Bounder::set(const double& h){
    if (!varTables.count(h)) generate({h});
    for (const auto& p : varTables[h]){
        spq->setVariable(p.first, p.second);
    }
}