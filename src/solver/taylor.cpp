#include <iostream>
#include <fmt/core.h>
#include <omp.h>
#include <numeric>
#include <algorithm>
#include <gurobi_c.h>

#include "taylor.hpp"
#include "core/kde.hpp"
#include "core/checker.hpp"
#include "util/uconfig.hpp"

using std::cerr;
using std::dynamic_pointer_cast;
using std::min;
using std::to_string;
using std::make_unique;
using std::make_shared;
using std::fill;
using std::iota;
using std::copy;

map<string, Option> taylorOptions = {
    {"number_of_cores", static_cast<int>(Config::getInstance()->nPhysicalCores)}, 
    {"soft_deterministic_constraint", false}, 
    {"time_limit", 60.0},
    {"max_number_of_iterations", 50}
};

const double Taylor::adjustCoef = Config::getInstance()->pt.get<double>("parameters.sample_adjustment_coefficient");

Taylor::Taylor(shared_ptr<StochasticPackageQuery> spq, const vector<long long>& ids, const SolIndType& initSol, const map<string, Option>& options): spq(spq), options(options){
    stat = make_unique<Stat>();
    for (const auto& p : taylorOptions){
        if (!options.count(p.first)) this->options[p.first] = p.second;
    }
    isSoftDetCon = boost::get<bool>(this->options.at("soft_deterministic_constraint"));
    nMaxIters = boost::get<int>(options.at("max_number_of_iterations"));
    nCores = boost::get<int>(this->options["number_of_cores"]);
    status = TaylorStatus::not_yet_found;

    auto tableSize = stat->pg->getTableSize(spq->tableName);
    for (const auto& id : ids){
        if (id > tableSize || id < 1){
            cerr << fmt::format("Cannot process id '{}' in table '{}'\n", id, spq->tableName);
            exit(1);
        }
    }
    idx = make_unique<UniqueIndexer>(nCores, tableSize, ids);
    size_t n = idx->size();
    shared_ptr<BoundConstraint> boundCon;
    shared_ptr<ProbConstraint> probCon;
    shared_ptr<AttrConstraint> attrCon;
    for (const auto& con : spq->cons){
        if (isStochastic(con, probCon, attrCon) && attrCon){
            size_t N = stat->pg->getColumnLength(spq->tableName, attrCon->attr);
            stoXs.emplace_back(N, 0);
        }
        if (isDeterministic(con, boundCon, attrCon)){
            detXs.push_back(0);
            if (attrCon){
                vector<double> attrs (n);
                #pragma omp parallel num_threads(nCores)
                {
                    auto coreIndex = omp_get_thread_num();
                    const auto& interval = idx->getInterval(coreIndex);
                    size_t n_ = interval.second - interval.first;
                    vector<double> attrs_; attrs_.reserve(n_);
                    Stat stat_;
                    stat_.getDetAttrs(spq->tableName, attrCon->attr, idx->getSql(coreIndex), attrs_);
                    copy(attrs_.begin(), attrs_.end(), attrs.begin()+interval.first);
                }
                detCons.emplace_back(attrs);
            }
            if (getCount(con)){
                detCons.emplace_back(n, 1.0);
            }
        }
    }

    objValue = 0;
    obj.resize(n);
    fill(obj.begin(), obj.end(), 0);
    objSense = ObjectiveSense::minimize;
    if (spq->obj){
        objSense = spq->obj->objSense;
        shared_ptr<AttrObjective> attrObj;
        if (isDeterministic(spq->obj, attrObj)){
            if (attrObj){
                #pragma omp parallel num_threads(nCores)
                {
                    auto coreIndex = omp_get_thread_num();
                    const auto& interval = idx->getInterval(coreIndex);
                    size_t n_ = interval.second - interval.first;
                    vector<double> attrs_; attrs_.reserve(n_);
                    Stat stat_;
                    stat_.getDetAttrs(spq->tableName, attrObj->obj, idx->getSql(coreIndex), attrs_);
                    copy(attrs_.begin(), attrs_.end(), obj.begin()+interval.first);
                }
            }
            if (getCount(spq->obj)) fill(obj.begin(), obj.end(), 1.0);
        }
        if (objSense == ObjectiveSense::maximize) bestObj = NEG_INF;
        else bestObj = POS_INF;
    }
    update(initSol);
}

void Taylor::update(const SolIndType& step){
    if (!step.size()) return;
    auto n = idx->size();
    auto m = spq->cons.size();
    size_t detInd = 0;
    size_t stoInd = 0;
    shared_ptr<ProbConstraint> probCon;
    shared_ptr<AttrConstraint> attrCon;

    for (size_t j = 0; j < m; ++j){
        const auto& con = spq->cons[j];
        if (isStochastic(con, probCon, attrCon) && attrCon){
            size_t N = stat->pg->getColumnLength(spq->tableName, attrCon->attr);
            vector<double> samples; samples.reserve(N);
            for (const auto& p : step){
                stat->getSamples(spq->tableName, attrCon->attr, idx->at(p.first), samples);
                for (size_t k = 0; k < N; ++k){
                    stoXs[stoInd][k] += samples[k]*p.second;
                }
            }
            stoInd ++;
        }
        if (isDeterministic(con)){
            for (const auto& p : step) detXs[detInd] += detCons[detInd][p.first]*p.second;
            detInd ++;
        }
    }
    // vector<size_t> inds; inds.reserve(step.size());
    // for (const auto& p : step) inds.push_back(p.first);
    // vector<size_t> stoInds (m);
    // size_t stoInd = 0, detInd = 0;
    // for (size_t j = 0; j < m; ++j){
    //     stoInds[j] = stoInd;
    //     const auto& con = spq->cons[j];
    //     if (isStochastic(con)) stoInd ++;
    //     if (isDeterministic(con)){
    //         for (const auto& p : step) detXs[detInd] += detCons[detInd][p.first]*p.second;
    //         detInd ++;
    //     }
    // }
    // for (const auto& p : step) objValue += obj[p.first]*p.second;
    // #pragma omp parallel num_threads(nCores)
    // {
    //     Stat stat_;
    //     shared_ptr<ProbConstraint> probCon;
    //     shared_ptr<AttrConstraint> attrCon;
    //     #pragma omp for schedule(dynamic) collapse(2)
    //     for (size_t j = 0; j < m; ++j){
    //         for (size_t i = 0; i < inds.size(); ++i){
    //             const auto& con = spq->cons[j];
    //             if (isStochastic(con, probCon, attrCon) && attrCon){
    //                 size_t N = stat_.pg->getColumnLength(spq->tableName, attrCon->attr);
    //                 long long id = idx->at(inds[i]);
    //                 vector<double> samples; samples.reserve(N);
    //                 stat_.getSamples(spq->tableName, attrCon->attr, id, samples);
    //                 auto m = step.at(inds[i]);
    //                 for (size_t k = 0; k < N; ++k){
    //                     #pragma omp atomic
    //                     stoXs[stoInds[j]][k] += samples[k]*m;
    //                 }
    //             }
    //         }
    //     }
    // }
}

void Taylor::solve(SolIndType& nextSol){
    size_t nInds = idx->size();
    size_t numvars = nInds + spq->getCvarCount() + spq->isStochasticObjective();
    size_t m = spq->cons.size() + spq->isStochasticObjective();
    size_t nzN = sol.size();
    vector<int> cbeg, cind;
    vector<char> sense; sense.reserve(2*m);
    vector<double> cval, rhs;
    cbeg.reserve(2*m+1); cind.reserve(2*m*numvars); cval.reserve(2*m*numvars); rhs.reserve(2*m);
    cbeg.push_back(0);
    deb(m);
}

void Taylor::solve(){
    SolIndType nextSol; solve(nextSol);
}

// SolType Taylor::getSol(const SolIndType& sol) const{
//     SolType res;
//     auto n = idx->size();
//     if (sol.size()){
//         for (const auto& p : sol){
//             if (p.first < n && p.first >= 0) res[idx->at(p.first)] = p.second;
//         }
//     } else{
//         for (const auto& p : bestSol){
//             if (p.first < n && p.first >= 0) res[idx->at(p.first)] = p.second;
//         }
//     }
//     return res;
// }