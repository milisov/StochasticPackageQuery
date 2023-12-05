#include <iostream>
#include <fmt/core.h>
#include <omp.h>
#include <numeric>
#include <algorithm>
#include <gurobi_c.h>

#include "taylor.hpp"
#include "core/kde.hpp"
#include "core/optim.hpp"
#include "core/checker.hpp"
#include "util/uconfig.hpp"
#include "util/unumeric.hpp"

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
    objSense = ObjectiveSense::minimize;
    if (spq->obj){
        objSense = spq->obj->objSense;
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
        if (objSense == ObjectiveSense::maximize) bestObj = NEG_INF;
        else bestObj = POS_INF;
    }

    isSoftDetCon = boost::get<bool>(this->options.at("soft_deterministic_constraint"));
    isDependencyVar = boost::get<bool>(this->options.at("dependency_var"));
    nMaxIters = boost::get<int>(options.at("max_number_of_iterations"));
    iter = 0;
    status = TaylorStatus::not_yet_found;
}

void Taylor::solve(SolIndType& nextSol){
    iter ++;
    int nCores = boost::get<int>(this->options.at("number_of_cores"));
    size_t nInds = ids.size();
    vector<int> inds (nInds);
    iota(inds.begin(), inds.end(), 0);
    size_t n = nInds + spq->nCvar + spq->isStoObj;
    size_t m = spq->cons.size() + spq->isStoObj;
    size_t nzN = sol.size();
    vector<double> lb = vector<double>(n, 0);
    vector<double> ub = vector<double>(n, POS_INF);
    // if (spq->repeat != StochasticPackageQuery::NO_REPEAT) fill(ub.begin(), ub.begin()+nInds, spq->repeat+1);
    if (maxSolutionSize > 0){
        double stepSize = maxSolutionSize / nMaxIters;
        for (size_t i = 0; i < nInds; ++i){
            if (sol.count(i)){
                lb[i] = max(0.0, sol.at(i)-stepSize);
                ub[i] = sol.at(i)+stepSize;
            } else ub[i] = stepSize;
        }
        if (spq->repeat != StochasticPackageQuery::NO_REPEAT){
            for (size_t i = 0; i < nInds; ++i) ub[i] = min(ub[i], spq->repeat+1.0);
        }
    } else fill(ub.begin(), ub.begin()+nInds, 1.0);
    vector<int> cbeg, cind;
    vector<double> violations (m, 0), detMuls (m, 0), Fxvs, cval, rhs;
    vector<char> sense; sense.reserve(2*m);
    cbeg.reserve(2*m+1); cind.reserve(2*m*n); cval.reserve(2*m*n); rhs.reserve(2*m);
    cbeg.push_back(0);
    vector<shared_ptr<KDE>> kdes;
    vector<double> softObj = obj;
    size_t detInd = 0, stoInd = 0;
    shared_ptr<ProbConstraint> probCon;
    shared_ptr<BoundConstraint> boundCon;
    for (const auto& con : spq->cons){
        if (isStochastic(con, probCon)){
            auto kde = make_shared<KDE>(stoXs[stoInd], true);
            Fxvs.push_back(kde->getQuickCdf(spq->getValue(probCon->v)));
            // deb(stoXs[stoInd], spq->getValue(probCon->v));
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
                        detMuls[ind] = 1 - 2*to_index(objSense);
                    }
                    auto ub = spq->getValue(boundCon->ub);
                    if (ub != POS_INF && ub < res){
                        violations[ind] = detViolate((res-ub)/detNorms[detInd]*sqn);
                        detMuls[ind] = 2*to_index(objSense) - 1;
                    }
                }
            }
            detInd ++;
        }
    }
    deb(iter, violations, objValue);
    double vio = std::accumulate(violations.begin(), violations.end(), 0.0);
    if (vio > 0) fill(softObj.begin(), softObj.end(), 0.0);
    update(vio, objValue, sol);
    double logitSum = logLogit(violations);
    for (size_t i = 0; i < violations.size(); ++ i){
        if (violations[i] > 0) violations[i] = exp(logit(violations[i])-logitSum);
    }
    size_t rowInd = 0, nzInd = 0;
    detInd = 0; stoInd = 0;
    shared_ptr<AttrConstraint> attrCon;
    for (const auto& con : spq->cons){
        double multiplier = violations[detInd+stoInd];
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
                if (!sameSense(stoIneq, objSense)) multiplier *= -1;
                double stoNorm = norm(stoCon);
                for (size_t i = 0; i < nInds; ++i) softObj[i] += multiplier*stoCon[i]/stoNorm;
            } else{
                auto preNzInd = nzInd;
                nzInd += nInds + (getCvar(con) != nullptr);
                cbeg.push_back(nzInd);
                cind.resize(nzInd); cval.resize(nzInd);
                copy(inds.begin(), inds.end(), cind.begin()+preNzInd);
                copy(stoCon.begin(), stoCon.end(), cval.begin()+preNzInd);
                rhs.push_back(stoBound);
                if (stoIneq == Inequality::lteq) sense.push_back(GRB_LESS_EQUAL);
                else if (stoIneq == Inequality::gteq) sense.push_back(GRB_GREATER_EQUAL);
                rowInd ++;
            }
            stoInd ++;
        }
        if (isDeterministic(con, boundCon)){
            if (multiplier > 0){
                multiplier *= detMuls[detInd+stoInd];
                for (size_t i = 0; i < nInds; ++i) softObj[i] += multiplier*detCons[detInd][i]/detNorms[detInd];
            } else{
                auto lbValue = spq->getValue(boundCon->lb);
                if (lbValue != NEG_INF){
                    auto preNzInd = nzInd;
                    nzInd += nInds;
                    cbeg.push_back(nzInd);
                    cind.resize(nzInd); cval.resize(nzInd);
                    copy(inds.begin(), inds.end(), cind.begin()+preNzInd);
                    copy(detCons[detInd].begin(), detCons[detInd].end(), cval.begin()+preNzInd);
                    rhs.push_back(lbValue);
                    sense.push_back(GRB_GREATER_EQUAL);
                    rowInd ++;
                }
                auto ubValue = spq->getValue(boundCon->ub);
                if (ubValue != POS_INF){
                    auto preNzInd = nzInd;
                    nzInd += nInds;
                    cbeg.push_back(nzInd);
                    cind.resize(nzInd); cval.resize(nzInd);
                    copy(inds.begin(), inds.end(), cind.begin()+preNzInd);
                    copy(detCons[detInd].begin(), detCons[detInd].end(), cval.begin()+preNzInd);
                    rhs.push_back(ubValue);
                    sense.push_back(GRB_LESS_EQUAL);
                    rowInd ++;
                }
            }
            detInd ++;
        }
    }

	GRBenv* env = NULL;
	GRBmodel* model = NULL;
	ckg(GRBemptyenv(&env), env);
	// ckg(GRBsetintparam(env, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_CONSERVATIVE), env);
	ckg(GRBsetintparam(env, GRB_INT_PAR_SIFTING, 2), env);
    ckg(GRBsetintparam(env, GRB_INT_PAR_LPWARMSTART, 1), env);
	// ckg(GRBsetintparam(env, GRB_INT_PAR_METHOD, GRB_METHOD_DUAL), env);
    ckg(GRBsetintparam(env, GRB_INT_PAR_METHOD, GRB_METHOD_CONCURRENT), env);
	ckg(GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0), env);
	ckg(GRBstartenv(env), env);
    vector<char> vtype (n, GRB_CONTINUOUS);
	ckgb(GRBnewmodel(env, &model, NULL, n, softObj.data(), lb.data(), ub.data(), vtype.data(), NULL), env, model);
    auto modelSense = objSense == ObjectiveSense::maximize ? GRB_MAXIMIZE : GRB_MINIMIZE;
	ckgb(GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, modelSense), env, model);
	ckgb(GRBaddconstrs(model, rowInd, cind.size(), cbeg.data(), cind.data(), cval.data(), sense.data(), rhs.data(), NULL), env, model);
    
	ckgb(GRBupdatemodel(model), env, model);
    if (preViolations.size()){
        bool isWarmStart = true;
        for (size_t i = 0; i < preViolations.size(); ++i){
            if (preViolations[i] > 0 && violations[i] == 0) isWarmStart = false;
            if (preViolations[i] == 0 && violations[i] > 0) isWarmStart = false;
            if (!isWarmStart) break;
        }
        if (isWarmStart){
            ckgb(GRBsetintattrarray(model, GRB_INT_ATTR_VBASIS, 0, n, vStart.data()), env, model);
            ckgb(GRBsetintattrarray(model, GRB_INT_ATTR_CBASIS, 0, rowInd, cStart.data()), env, model);
            ckgb(GRBsetdblattrarray(model, GRB_DBL_ATTR_PSTART, 0, n, pStart.data()), env, model);
            ckgb(GRBsetdblattrarray(model, GRB_DBL_ATTR_DSTART, 0, rowInd, dStart.data()), env, model);
        }
    }
    preViolations = violations;
    CLK(pro, "a");
	ckgb(GRBoptimize(model), env, model);
    STP(pro, "a");
	int status;
	ckgb(GRBgetintattr(model, GRB_INT_ATTR_STATUS, &status), env, model);
	if (status == GRB_OPTIMAL){
        pStart.resize(n);
    	ckgb(GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, n, pStart.data()), env, model);
        for (size_t i = 0; i < n; ++i) {
            if (pStart[i] > 0){
                nextSol[i] = pStart[i];
            }
        }
        double solutionSize = 0;
        for (const auto& p : nextSol){
            if (p.first < nInds) solutionSize += p.second;
        }
        maxSolutionSize = max(maxSolutionSize, solutionSize);

        vStart.resize(n); cStart.resize(rowInd);
        ckgb(GRBgetintattrarray(model, GRB_INT_ATTR_VBASIS, 0, n, vStart.data()), env, model);
        ckgb(GRBgetintattrarray(model, GRB_INT_ATTR_CBASIS, 0, rowInd, cStart.data()), env, model);
        dStart.resize(rowInd);
        ckgb(GRBgetdblattrarray(model, GRB_DBL_ATTR_PI, 0, rowInd, dStart.data()), env, model);
        
        // vector<size_t> basics;
        // getBasicVariables(model, n, basics);
        // deb(nextSol, basics);
	} else{
        cerr << getGurobiStatus(status) << '\n';
        status = minVio == 0 ?  TaylorStatus::found : TaylorStatus::not_found;
    }
    if (model) GRBfreemodel(model);
    if (env) GRBfreeenv(env);
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
    // deb(sol, step, objValue);
    // sol = add(1.0, sol, 1.0, step);
    for (const auto& p : step){
        double val = p.second;
        if (sol.count(p.first)) val += sol.at(p.first);
        if (val > MACHINE_EPS) sol[p.first] = val;
        else sol.erase(p.first);
    }
    size_t hashed = hashSol(sol);
    if (hashedSols.count(hashed)){
        bool isCycled = false;
        for (const auto& hashedSol : hashedSols.at(hashed)){
            if (isEqual(sol, hashedSol)){
                isCycled = true;
                break;
            }
        }
        if (!isCycled) hashedSols[hashed].emplace_back(sol);
        else status = minVio > 0 ? TaylorStatus::cycled : TaylorStatus::found;
    } else hashedSols[hashed] = {sol};
}

void Taylor::update(const double& vio, const double& objValue, const SolIndType& sol){
    if (minVio > vio){
        minVio = vio;
        bestObj = objValue;
        bestSol = sol;
    } else if (minVio == vio && spq->obj){
        if (objSense == ObjectiveSense::maximize && bestObj < objValue){
            bestObj = objValue;
            bestSol = sol;
        } else if (objSense == ObjectiveSense::minimize && bestObj > objValue){
            bestObj = objValue;
            bestSol = sol;
        }
    }
}

void Taylor::solve(){
    unique_ptr<Optim> optim = make_unique<Direct>();
    for (int i = 0; i < nMaxIters; ++i){
        CLK(pro, "b");
        SolIndType nextSol; solve(nextSol);
        update(optim->towards(nextSol));
        STP(pro, "b");
        PRINT(pro);
        if (status != TaylorStatus::not_yet_found) break;
    }
    if (status == TaylorStatus::not_yet_found) status = minVio == 0 ? TaylorStatus::found : TaylorStatus::not_found;
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