#include <Highs.h>
#include <iostream>
#include <fmt/core.h>
#include <omp.h>
#include <numeric>
#include <algorithm>

#include "taylor.hpp"
#include "core/kde.hpp"
#include "core/optim.hpp"
#include "core/checker.hpp"
#include "util/udebug.hpp"
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

#define DYNAMIC_CHUNK 10

map<string, Option> defaultOptions = {
    {"number_of_cores", static_cast<int>(Config::getInstance()->nPhysicalCores)}, 
    {"soft_deterministic_constraint", false}, 
    {"dependency_var", false}, 
    {"time_limit", 60.0},
    {"max_number_of_iterations", 50}
};

Taylor::Taylor(shared_ptr<StochasticPackageQuery> spq, const vector<long long>& ids, const map<string, Option>& options): spq(spq), ids(ids), options(options){
    stat = make_unique<Stat>();
    for (const auto& p : defaultOptions){
        if (!options.count(p.first)) this->options[p.first] = p.second;
    }
    auto tableSize = stat->pg->getTableSize(spq->tableName);
    for (const auto& id : ids){
        if (id > tableSize){
            cerr << fmt::format("Cannot process id '{}' in table '{}'\n", id, spq->tableName);
            exit(1);
        }
    }
    vector<string> strIds (ids.size());
    for (size_t i = 0; i < ids.size(); ++ i) strIds[i] = to_string(this->ids[i]);
    string joinId = boost::join(strIds, ",");
    if (!ids.size()){
        this->ids.resize(tableSize);
        iota(this->ids.begin(), this->ids.end(), 1LL);
    }
    size_t n = this->ids.size();
    sqn = sqrt(n);
    int nCores = boost::get<int>(this->options["number_of_cores"]);
    shared_ptr<BoundConstraint> boundCon;
    shared_ptr<ProbConstraint> probCon;
    shared_ptr<AttrConstraint> attrCon;
    for (const auto& con : spq->cons){
        if (isStochastic(con, probCon, attrCon) && attrCon){
            size_t N = stat->pg->getColumnLength(spq->tableName, attrCon->attr);
            stoXs.emplace_back(N, 0);
            vector<double> means;
            if (!ids.size()){
                means.resize(n);
                auto intervals = divideInterval(1, n, nCores);
                #pragma omp parallel num_threads(nCores)
                {
                    auto coreIndex = omp_get_thread_num();
                    string betweenId = fmt::format("BETWEEN {} AND {}", intervals[coreIndex], intervals[coreIndex+1]-1);
                    size_t n_ = intervals[coreIndex+1]-intervals[coreIndex];
                    vector<double> means_; means_.reserve(n_);
                    vector<double> _; _.reserve(n_);
                    Stat stat_;
                    stat_.getStoMeanVars(spq->tableName, attrCon->attr, betweenId, means_, _);
                    copy(means_.begin(), means_.end(), means.begin()+intervals[coreIndex]-1);
                }
            } else {
                vector<double> _; _.reserve(n);
                means.reserve(n);
                stat->getStoMeanVars(spq->tableName, attrCon->attr, joinId, means, _);
            }
            stoMeans.emplace_back(means);
        }
        if (isDeterministic(con, boundCon, attrCon)){
            detXs.push_back(0);
            if (attrCon){
                vector<double> attrs;
                if (!ids.size()){
                    attrs.resize(n);
                    auto intervals = divideInterval(1, n, nCores);
                    #pragma omp parallel num_threads(nCores)
                    {
                        auto coreIndex = omp_get_thread_num();
                        string betweenId = fmt::format("BETWEEN {} AND {}", intervals[coreIndex], intervals[coreIndex+1]-1);
                        size_t n_ = intervals[coreIndex+1]-intervals[coreIndex];
                        vector<double> attrs_; attrs_.reserve(n_);
                        Stat stat_;
                        stat_.getDetAttrs(spq->tableName, attrCon->attr, betweenId, attrs_);
                        copy(attrs_.begin(), attrs_.end(), attrs.begin()+intervals[coreIndex]-1);
                    }
                } else {
                    attrs.reserve(n);
                    stat->getDetAttrs(spq->tableName, attrCon->attr, joinId, attrs);
                }
                detCons.emplace_back(attrs);
                detNorms.push_back(norm(attrs));
            }
            if (getCount(con)){
                detCons.emplace_back(n, 1.0);
                detNorms.push_back(sqn);
            }
        }
    }
    objValue = 0;
    minVio = POS_INF;
    maxSolutionSize = 0;
    bestObj = 0;
    obj.resize(n);
    fill(obj.begin(), obj.end(), 0);
    if (spq->obj){
        shared_ptr<AttrObjective> attrObj;
        if (isDeterministic(spq->obj, attrObj)){
            if (attrObj){
                vector<double> attrs;
                if (!ids.size()){
                    attrs.resize(n);
                    auto intervals = divideInterval(1, n, nCores);
                    #pragma omp parallel num_threads(nCores)
                    {
                        auto coreIndex = omp_get_thread_num();
                        string betweenId = fmt::format("BETWEEN {} AND {}", intervals[coreIndex], intervals[coreIndex+1]-1);
                        size_t n_ = intervals[coreIndex+1]-intervals[coreIndex];
                        vector<double> attrs_; attrs_.reserve(n_);
                        Stat stat_;
                        stat_.getDetAttrs(spq->tableName, attrObj->obj, betweenId, attrs_);
                        copy(attrs_.begin(), attrs_.end(), attrs.begin()+intervals[coreIndex]-1);
                    }
                } else {
                    attrs.reserve(n);
                    stat->getDetAttrs(spq->tableName, attrObj->obj, joinId, attrs);
                }
                copy(attrs.begin(), attrs.end(), obj.begin());
            }
            if (getCount(spq->obj)) fill(obj.begin(), obj.end(), 1.0);
        }
        if (spq->obj->objSense == ObjectiveSense::maximize) bestObj = NEG_INF;
        else bestObj = POS_INF;
    }
    objNorm = norm(obj);

    isSoftDetCon = boost::get<bool>(this->options.at("soft_deterministic_constraint"));
    isDependencyVar = boost::get<bool>(this->options.at("dependency_var"));
    nIters = boost::get<int>(options.at("max_number_of_iterations"));
}

void Taylor::solve(SolIndType& nextSol){
    int nCores = boost::get<int>(this->options.at("number_of_cores"));
    size_t nInds = ids.size();
    size_t n = nInds + spq->nCvar + spq->isStoObj;
    vector<int> inds (nInds);
    iota(inds.begin(), inds.end(), 0);
    size_t m = spq->cons.size() + spq->isStoObj;
    size_t nzN = sol.size();
    HighsModel model;
    model.lp_.sense_ = ObjSense::kMinimize;
    model.lp_.offset_ = 0;
    model.lp_.a_matrix_.format_ = MatrixFormat::kRowwise;
    model.lp_.num_col_ = n;
    model.lp_.col_lower_ = vector<double>(n, 0);
    model.lp_.col_upper_ = vector<double>(n, POS_INF);
    if (spq->repeat != StochasticPackageQuery::NO_REPEAT) fill(model.lp_.col_upper_.begin(), model.lp_.col_upper_.begin()+nInds, spq->repeat+1);
    if (maxSolutionSize > 0){
        double stepSize = (maxSolutionSize - sum(sol)) / nIters;
        deb(stepSize, maxSolutionSize, sum(sol), nIters);
        assert(stepSize > 0);
        for (size_t i = 0; i < nInds; ++i){
            if (sol.count(i)){
                model.lp_.col_lower_[i] = max(0.0, sol.at(i)-stepSize);
                model.lp_.col_upper_[i] = sol.at(i)+stepSize;
            } else model.lp_.col_upper_[i] = stepSize;
        }
        if (spq->repeat != StochasticPackageQuery::NO_REPEAT){
            for (size_t i = 0; i < nInds; ++i) model.lp_.col_upper_[i] = min(model.lp_.col_upper_[i], spq->repeat+1.0);
        }
    }
    vector<int> start_, index_;
    vector<double> violations (m, 0), detMuls (m, 0), Fxvs, value_;
    start_.reserve(m); index_.reserve(m*model.lp_.num_col_); value_.reserve(m*model.lp_.num_col_);
    start_.push_back(0);
    vector<shared_ptr<KDE>> kdes;
    vector<double> softObj (n, 0);
    if (objNorm > 0){
        if (!sameSense(spq->obj->objSense, model.lp_.sense_)){
            for (size_t i = 0; i < n; ++i) softObj[i] = -obj[i]/objNorm;
        } else for (size_t i = 0; i < n; ++i) softObj[i] = obj[i]/objNorm;
    }
    size_t detInd = 0, stoInd = 0;
    shared_ptr<ProbConstraint> probCon;
    shared_ptr<BoundConstraint> boundCon;
    for (const auto& con : spq->cons){
        if (isStochastic(con, probCon)){
            auto kde = make_shared<KDE>(stoXs[stoInd], true);
            Fxvs.push_back(kde->getQuickCdf(spq->getValue(probCon->v)));
            if (getVar(con)){
                double vio = Fxvs.back();
                double p = spq->getValue(probCon->p);
                if (probCon->vsign == Inequality::lteq) vio -= p;
                else if (probCon->vsign == Inequality::gteq) vio -= (1-p);
                if (probCon->vsign != probCon->psign) vio = -vio;
                double maxVio = p;
                if (probCon->psign == Inequality::lteq) maxVio = 1-p;
                violations[detInd+stoInd] = max(0.0, vio) / maxVio;
            }
            kdes.push_back(kde);
            stoInd ++; 
        }
        if (isDeterministic(con, boundCon)){
            if (isSoftDetCon){
                shared_ptr<BoundConstraint> boundCon = dynamic_pointer_cast<BoundConstraint>(con);
                if (boundCon){
                    double res = detXs[detInd];
                    auto lb = spq->getValue(boundCon->lb);
                    auto ind = detInd+stoInd;
                    if (lb != NEG_INF && lb > res){
                        violations[ind] = detViolate((lb-res)/detNorms[detInd]*sqn);
                        detMuls[ind] = -static_cast<double>(model.lp_.sense_);
                    }
                    auto ub = spq->getValue(boundCon->ub);
                    if (ub != POS_INF && ub < res){
                        violations[ind] = detViolate((res-ub)/detNorms[detInd]*sqn);
                        detMuls[ind] = static_cast<double>(model.lp_.sense_);
                    }
                }
            }
            detInd ++;
        }
    }
    deb(violations);
    double vio = std::accumulate(violations.begin(), violations.end(), 0.0);
    update(vio, objValue, sol);
    double logitSum = logLogit(violations);
    for (size_t i = 0; i < violations.size(); ++ i){
        if (violations[i] > 0) violations[i] = exp(logit(violations[i])-logitSum);
    }
    size_t rowInd = 0, nzInd = 0;
    detInd = 0; stoInd = 0;
    shared_ptr<AttrConstraint> attrCon;
    for (const auto& con : spq->cons){
        double multiplier = violations[detInd+stoInd]*sqn;
        if (isStochastic(con, probCon, attrCon) && attrCon){
            double v = spq->getValue(probCon->v);
            double p = spq->getValue(probCon->p);
            vector<double> stoCon (nInds);
            double stoBound;
            Inequality stoIneq = Inequality::lteq;
            if (getVar(con)){
                size_t N = stat->pg->getColumnLength(spq->tableName, attrCon->attr);
                if (isDependencyVar){
                    
                } else{
                    const auto& kde = kdes[stoInd];
                    double fxv = max(kde->getPdf(v), MACHINE_EPS);
                    vector<double> dF (nInds);
                    #pragma omp parallel num_threads(nCores)
                    {
                        Stat stat_;
                        #pragma omp for schedule(dynamic, DYNAMIC_CHUNK)
                        for (size_t i = 0; i < nInds; ++i){
                            if (!sol.count(i)){
                                dF[i] = -fxv*stoMeans[stoInd][i];
                            } else{
                                vector<double> samples, quantiles, XiBar (N); 
                                samples.reserve(N); quantiles.reserve(2*N+1);
                                stat_.getSamples(spq->tableName, attrCon->attr, ids[i], samples);
                                stat_.getQuantiles(spq->tableName, attrCon->attr, ids[i], quantiles);
                                for (size_t j = 0; j < N; ++j) XiBar[j] = stoXs[stoInd][j] - sol.at(i)*samples[j];
                                KDE kdeXiBar (XiBar, true);
                                dF[i] = kdeXiBar.convolve(quantiles, v, sol.at(i));
                            }
                            stoCon[i] = dF[i];
                        }
                    }
                    stoBound = p;
                    if (probCon->vsign == Inequality::gteq) stoBound = 1-p;
                    double dot = 0;
                    for (const auto& p : sol) dot += p.second*dF[p.first];
                    stoBound += dot-Fxvs[stoInd];
                    if (probCon->vsign != probCon->psign) stoIneq = Inequality::gteq;
                }
            }
            if (multiplier > 0){
                if (!sameSense(stoIneq, model.lp_.sense_)) multiplier *= -1;
                double stoNorm = norm(stoCon);
                for (size_t i = 0; i < nInds; ++i) softObj[i] += multiplier*stoCon[i]/stoNorm;
            } else{
                auto preNzInd = nzInd;
                nzInd += n;
                start_.push_back(nzInd);
                index_.resize(nzInd); value_.resize(nzInd);
                copy(inds.begin(), inds.end(), index_.begin()+preNzInd);
                copy(stoCon.begin(), stoCon.end(), value_.begin()+preNzInd);
                if (stoIneq == Inequality::lteq){
                    model.lp_.row_lower_.push_back(NEG_INF);
                    model.lp_.row_upper_.push_back(stoBound);
                } else if (stoIneq == Inequality::gteq){
                    model.lp_.row_lower_.push_back(stoBound);
                    model.lp_.row_upper_.push_back(POS_INF);
                }
                rowInd ++;
            }
            stoInd ++;
        }
        if (isDeterministic(con, boundCon)){
            if (multiplier > 0){
                multiplier *= detMuls[detInd+stoInd];
                for (size_t i = 0; i < nInds; ++i) softObj[i] += multiplier*detCons[detInd][i]/detNorms[detInd];
            } else{
                auto preNzInd = nzInd;
                nzInd += n;
                start_.push_back(nzInd);
                index_.resize(nzInd); value_.resize(nzInd);
                copy(inds.begin(), inds.end(), index_.begin()+preNzInd);
                copy(detCons[detInd].begin(), detCons[detInd].end(), value_.begin()+preNzInd);
                model.lp_.row_lower_.push_back(spq->getValue(boundCon->lb));
                model.lp_.row_upper_.push_back(spq->getValue(boundCon->ub));
                rowInd ++;
            }
            detInd ++;
        }
    }
    model.lp_.num_row_ = rowInd;
    model.lp_.col_cost_ = softObj;
    model.lp_.a_matrix_.start_ = start_;
    model.lp_.a_matrix_.index_ = index_;
    model.lp_.a_matrix_.value_ = value_;

    Highs highs;
	// highs.setOptionValue("output_flag", false);
	// highs.setOptionValue("log_to_console", false);
	highs.setOptionValue("presolve", "off");
	highs.setOptionValue("random_seed", abs(static_cast<int>(Config::getInstance()->seed())));
	highs.setOptionValue("simplex_strategy", 0);
	highs.setOptionValue("infinite_bound", POS_INF);
    highs.setOptionValue("small_matrix_value", MACHINE_EPS);
    auto status = highs.passModel(model);
    assert(status==HighsStatus::kOk || status==HighsStatus::kWarning);
    const HighsLp& lp = highs.getLp();
    status = highs.run();
    assert(status==HighsStatus::kOk);
    const auto& solution = highs.getSolution();
    const bool hasValues = highs.getInfo().primal_solution_status;
    if (!hasValues){
        cerr << "No solutions found\n";
        exit(1);
    } else{
        for (size_t i=0; i < n; ++i) {
            if (solution.col_value[i] > 0){
                nextSol[i] = solution.col_value[i];
            }
        }
        double solutionSize = 0;
        for (const auto& p : nextSol){
            if (p.first < ids.size()) solutionSize += p.second;
        }
        maxSolutionSize = max(maxSolutionSize, solutionSize);
    }
}

void Taylor::update(const SolIndType& step){
    assert(step.size());
    int nCores = boost::get<int>(this->options.at("number_of_cores"));
    vector<size_t> inds; inds.reserve(step.size());
    for (const auto& p : step) inds.push_back(p.first);
    vector<size_t> stoInds (spq->cons.size());
    size_t stoInd = 0, detInd = 0;
    for (size_t j = 0; j < spq->cons.size(); ++j){
        stoInds[j] = stoInd;
        const auto& con = spq->cons[j];
        if (isStochastic(con)) stoInd ++;
        if (isDeterministic(con)){
            for (const auto& p : step) detXs[detInd] += detCons[detInd][p.first]*p.second;
            detInd ++;
        }
    }
    for (const auto& p : step) objValue += obj[p.first]*p.second;
    #pragma omp parallel num_threads(nCores)
    {
        Stat stat_;
        shared_ptr<ProbConstraint> probCon;
        shared_ptr<AttrConstraint> attrCon;
        #pragma omp for schedule(dynamic) collapse(2) nowait
        for (size_t j = 0; j < spq->cons.size(); ++j){
            for (size_t i = 0; i < inds.size(); ++i){
                const auto& con = spq->cons[j];
                if (isStochastic(con, probCon, attrCon) && attrCon){
                    size_t N = stat_.pg->getColumnLength(spq->tableName, attrCon->attr);
                    long long id = ids[inds[i]];
                    vector<double> samples; samples.reserve(N);
                    stat_.getSamples(spq->tableName, attrCon->attr, id, samples);
                    auto m = step.at(inds[i]);
                    for (size_t k = 0; k < N; ++k){
                        #pragma omp atomic
                        stoXs[stoInds[j]][k] += samples[k]*m;
                    }
                }
            }
        }
    }
    deb(sol, objValue);
    sol = add(1.0, sol, 1.0, step);
}

void Taylor::update(const double& vio, const double& objValue, const SolIndType& sol){
    if (minVio > vio){
        minVio = vio;
        bestObj = objValue;
        bestSol = sol;
    } else if (minVio == vio && spq->obj){
        if (spq->obj->objSense == ObjectiveSense::maximize && bestObj < objValue){
            bestObj = objValue;
            bestSol = sol;
        } else if (spq->obj->objSense == ObjectiveSense::minimize && bestObj > objValue){
            bestObj = objValue;
            bestSol = sol;
        }
    }
}

void Taylor::solve(){
    RMSprop optim;
    SPQChecker chk (spq);
    while (nIters > 0){
        SolIndType nextSol; solve(nextSol);
        deb(nextSol);
        update(optim.towards(nextSol));
        chk.display(getSol(sol));
        nIters --;
    }
    chk.display({{6205,0.953924},{9854,74.0461}});
}

SolType Taylor::getSol(const SolIndType& sol) const{
    SolType res;
    if (sol.size()){
        for (const auto& p : sol){
            if (p.first < ids.size() && p.first >= 0) res[ids[p.first]] = p.second;
        }
    } else{
        for (const auto& p : bestSol){
            if (p.first < ids.size() && p.first >= 0) res[ids[p.first]] = p.second;
        }
    }
    return res;
}