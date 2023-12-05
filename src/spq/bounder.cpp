#include <pcg_random.hpp>
#include <random>
#include <cmath>
#include <utility>
#include <iostream>
#include <fmt/core.h>
#include <omp.h>
#include <boost/variant.hpp>
#include <boost/math/special_functions/erf.hpp>

#include "bounder.hpp"
#include "util/udebug.hpp"
#include "util/uconfig.hpp"
#include "util/uio.hpp"
#include "util/unumeric.hpp"

using std::pair;
using std::cerr;

const double Bounder::hardLimit = -log(NUMERIC_EPS/(1-NUMERIC_EPS));
const double countRelaxationRatio = 0.5;

Bounder::Bounder(shared_ptr<StochasticPackageQuery> spq, const size_t& N, const double& E): spq(spq), N(N), E(E){
    PgManager pg;
    stat = std::make_unique<Stat>();
    long long n = pg.getTableSize(spq->tableName);
    size_t iE;
    double whole;
    double fE = modf(E, &whole);
    iE = static_cast<size_t>(whole)+1;
    pcg32 seedGen (Config::getInstance()->seed());
    auto nCores = Config::getInstance()->nPhysicalCores*2;
    vector<unsigned int> seeds (nCores);
    for (size_t i = 0; i < nCores; ++i) seeds[i] = seedGen();
    auto intervals = divideInterval(0, N-1, nCores);
    vector<string> joinIds (nCores);

    #pragma omp parallel num_threads(nCores)
    {
        auto coreIndex = omp_get_thread_num();
        pcg32 gen (seeds[coreIndex]);
        std::uniform_int_distribution<long long> dist (0LL, n-1);
        vector<string> strIds; strIds.reserve((intervals[coreIndex+1]-intervals[coreIndex])*iE);
        for (size_t i = intervals[coreIndex]; i < intervals[coreIndex+1]; ++i){
            for (size_t j = 0; j < iE; ++j){
                strIds.push_back(to_string(dist(gen)));
            }
        }
        joinIds[coreIndex] = boost::join(strIds, ",");
    }

    for (auto con : spq->cons){
        shared_ptr<BoundConstraint> boundCon;
        shared_ptr<ProbConstraint> probCon;
        shared_ptr<AttrConstraint> attrCon;
        if (isDeterministic(con, boundCon, attrCon) && attrCon){
            vector<double> G (N);
            #pragma omp parallel num_threads(nCores)
            {
                auto coreIndex = omp_get_thread_num();
                Stat stat;
                vector<double> attrs; attrs.reserve((intervals[coreIndex+1]-intervals[coreIndex])*iE);
                stat.getDetAttrs(spq->tableName, attrCon->attr, joinIds[coreIndex], attrs);
                for (size_t i = intervals[coreIndex]; i < intervals[coreIndex+1]; ++i){
                    double attr = 0;
                    size_t ind = (i-intervals[coreIndex])*iE;
                    for (size_t j = 0; j < iE-1; ++j) attr += attrs[ind+j];
                    attr += attrs[ind+iE-1]*fE;
                    G[i] = attr;
                }
            }
            detKDEs.push_back(std::make_unique<KDE>(G, true));
        }
        if (isStochastic(con, probCon, attrCon) && attrCon){
            double c = 0;
            if (getVar(con)){
                double p = spq->getValue(probCon->p);
                if (probCon->vsign == Inequality::gteq) p = 1-p;
                c = sqrt(2)*boost::math::erf_inv(2*p-1);
                bool reverse = false;
                if (probCon->vsign == probCon->psign) reverse = true;
                reverses.push_back(reverse);
            }
            double step = 0;
            vector<double> G (N);
            #pragma omp parallel num_threads(nCores)
            {
                auto coreIndex = omp_get_thread_num();
                Stat stat;
                size_t sz = (intervals[coreIndex+1]-intervals[coreIndex])*iE;
                vector<double> means; means.reserve(sz);
                vector<double> vars; vars.reserve(sz);
                stat.getStoMeanVars(spq->tableName, attrCon->attr, joinIds[coreIndex], means, vars);
                double step_ = 0;
                for (size_t i = intervals[coreIndex]; i < intervals[coreIndex+1]; ++i){
                    double mean = 0, var = 0;
                    size_t ind = (i-intervals[coreIndex])*iE;
                    for (size_t j = 0; j < iE-1; ++j){
                        mean += means[ind+j];
                        var += vars[ind+j];
                    }
                    mean += means[ind+iE-1]*fE;
                    var += vars[ind+iE-1]*fE;
                    G[i] = mean + sqrt(var)*c;
                    step_ += var;
                }
                #pragma omp atomic
                step += step_;
            }
            stoKDEs.push_back(std::make_unique<KDE>(G, true));
            step = sqrt(step/N)*abs(c);
            steps.push_back(step);
        }
    }
}

void Bounder::setBound(const double& h, const Bound& bound, const double& value){
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
		shared_ptr<BoundConstraint> boundCon;
        shared_ptr<ProbConstraint> probCon;
        shared_ptr<AttrConstraint> attrCon;
        if (isDeterministic(con, boundCon, attrCon)){
            if (getCount(con)){
                for (const auto& h : hards){
                    if (boundCon->lb.which() == 0) setBound(h, boundCon->lb, E*(1-countRelaxationRatio));
                    if (boundCon->ub.which() == 0) setBound(h, boundCon->ub, E*(1+countRelaxationRatio));
                }
            }
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
        if (isStochastic(con, probCon, attrCon) && attrCon && probCon->v.which() == 0){
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