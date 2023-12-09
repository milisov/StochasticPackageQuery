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

using std::cerr;
using std::dynamic_pointer_cast;
using std::min;
using std::to_string;
using std::make_unique;
using std::make_shared;
using std::fill;
using std::iota;
using std::copy;

map<string, Option> defaultOptions = {
    {"number_of_cores", static_cast<int>(Config::getInstance()->nPhysicalCores)}, 
    {"soft_deterministic_constraint", false}, 
    {"dependency_var", false}, 
    {"time_limit", 60.0},
    {"max_number_of_iterations", 50}
};

const double Taylor::adjustCoef = Config::getInstance()->pt.get<double>("parameters.sample_adjustment_coefficient");
const double budget = 1e6;

double getKernelRange(const vector<double>& info={}){
    return info.at(0)*sqrt(info.at(1));
}

double getKernelIntegral(const double& x, const vector<double>& info={}){
    return -0.125*x*x;
}

Taylor::Taylor(shared_ptr<StochasticPackageQuery> spq, const vector<long long>& ids, const map<string, Option>& options): spq(spq), options(options){
    stat = make_unique<Stat>();
    for (const auto& p : defaultOptions){
        if (!options.count(p.first)) this->options[p.first] = p.second;
    }
    isSoftDetCon = boost::get<bool>(this->options.at("soft_deterministic_constraint"));
    isDependentVar = boost::get<bool>(this->options.at("dependency_var"));
    nMaxIters = boost::get<int>(options.at("max_number_of_iterations"));
    nCores = boost::get<int>(this->options["number_of_cores"]);
    iter = 0;
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
    sqn = sqrt(n);
    shared_ptr<BoundConstraint> boundCon;
    shared_ptr<ProbConstraint> probCon;
    shared_ptr<AttrConstraint> attrCon;
    for (const auto& con : spq->cons){
        if (isStochastic(con, probCon, attrCon) && attrCon){
            size_t N = stat->pg->getColumnLength(spq->tableName, attrCon->attr);
            stoXs.emplace_back(N, 0);
            varXs.push_back(0.0);
            zeroPds.emplace_back(n, 0);
            nsix.push_back(pow(static_cast<double>(N), -1.0/6.0));
            vector<double> means (n), vars (n);
            #pragma omp parallel num_threads(nCores)
            {
                auto coreIndex = omp_get_thread_num();
                const auto& interval = idx->getInterval(coreIndex);
                size_t n_ = interval.second - interval.first;
                vector<double> means_; means_.reserve(n_);
                vector<double> vars_; vars_.reserve(n_);
                Stat stat_;
                stat_.getStoMeanVars(spq->tableName, attrCon->attr, idx->getSql(coreIndex), means_, vars_);
                copy(means_.begin(), means_.end(), means.begin()+interval.first);
                copy(vars_.begin(), vars_.end(), vars.begin()+interval.first);
            }
            stoMeans.emplace_back(means);
            stoVars.emplace_back(vars);
            double p = spq->getValue(probCon->p);
            if (getVar(con)){
                double rhsp = probCon->vsign == Inequality::lteq ? p : 1-p;
                double denom = 2*(adjustCoef*adjustCoef+N);
                double tmp = 2*rhsp-1;
                double delta = adjustCoef*sqrt(denom/2-N*tmp*tmp)/denom;
                double base = 0.5 + N*tmp/denom;
                double res = probCon->vsign == probCon->psign ? base-delta : base+delta;
                if (probCon->vsign == Inequality::gteq) res = 1-res;
                adjustments[boost::get<string>(probCon->p)] = {p, res};
            }
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
                detNorms.push_back(norm(attrs));
            }
            if (getCount(con)){
                detCons.emplace_back(n, 1.0);
                detNorms.push_back(sqn);
            }
        }
    }
    jvInds.resize(varXs.size());
    objValue = 0;
    minVio = POS_INF;
    maxSolutionSize = 0;
    gamma = 1.0;
    bestObj = 0;
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
}

void Taylor::solve(SolIndType& nextSol){
    CLOCK("pre");
    iter ++;
    size_t nInds = idx->size();
    size_t numvars = nInds + spq->nCvar + spq->isStoObj;
    size_t m = spq->cons.size() + spq->isStoObj;
    size_t nzN = sol.size();
    vector<double> lb = vector<double>(numvars, 0);
    vector<double> ub = vector<double>(numvars, POS_INF);
    for (size_t i = 0; i < nInds; ++i){
        if (sol.count(i)){
            lb[i] = max(0.0, sol.at(i)-gamma);
            ub[i] = sol.at(i)+gamma;
        } else ub[i] = gamma;
    }
    if (spq->repeat != StochasticPackageQuery::NO_REPEAT){
        for (size_t i = 0; i < nInds; ++i) ub[i] = min(ub[i], spq->repeat+1.0);
    }
    vector<int> cbeg, cind;
    vector<double> violations (m, 0), detMuls (m, 0), Fxvs, cval, rhs;
    vector<char> sense; sense.reserve(2*m);
    cbeg.reserve(2*m+1); cind.reserve(2*m*numvars); cval.reserve(2*m*numvars); rhs.reserve(2*m);
    cbeg.push_back(0);
    vector<shared_ptr<KDE>> kdes;
    vector<double> softObj = obj;
    size_t detInd = 0, stoInd = 0;
    shared_ptr<ProbConstraint> probCon;
    shared_ptr<BoundConstraint> boundCon;
    STOP("pre");
    CLOCK("vio");
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
    STOP("vio");
    CLOCK("posvio");
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
                CLOCK("sto");
                const auto& kde = kdes[stoInd];
                double fxv = kde->getPdf(v);
                vector<double> dF (nInds, 0);
                #pragma omp parallel num_threads(nCores)
                {
                    Stat stat_;
                    if (isDependentVar){
                        #pragma omp for schedule(dynamic)
                        for (size_t i = 0; i < nInds; ++i){
                            if (!sol.count(i)){
                                if (isEqual(fxv, 0)){
                                    dF[i] = -stoMeans[stoInd][i]/max(sqrt(stoVars[stoInd][i]), 1.0);
                                } else dF[i] = zeroPds[stoInd][i];
                            } else{
                                vector<double> samples, XiBar (N); 
                                samples.reserve(N);
                                stat_.getSamples(spq->tableName, attrCon->attr, idx->at(i), samples);
                                AccSet acc;
                                for (size_t j = 0; j < N; ++j){
                                    XiBar[j] = stoXs[stoInd][j] - sol.at(i)*samples[j];
                                    acc(XiBar[j]);
                                }
                                double r1 = getKernelRange({nsix[stoInd], static_cast<double>(ba::variance(acc))});
                                double r2 = getKernelRange({nsix[stoInd], stoVars[stoInd][i]});
                                for (size_t j = 0; j < N; ++j){
                                    double left = max(samples[j]-r2, (v-XiBar[j]-r1)/sol.at(i));
                                    double right = min(samples[j]+r2, (v-XiBar[j]+r1)/sol.at(i));
                                    if (left < right) dF[i] += (getKernelIntegral(right)-getKernelIntegral(left));
                                }
                            }
                            stoCon[i] = dF[i];
                        }
                    } else{
                        #pragma omp for schedule(dynamic)
                        for (size_t i = 0; i < nInds; ++i){
                            if (!sol.count(i)){
                                if (isEqual(fxv, 0)){
                                    dF[i] = -stoMeans[stoInd][i]/max(sqrt(stoVars[stoInd][i]), 1.0);
                                } else dF[i] = -fxv*stoMeans[stoInd][i];
                            } else{
                                vector<double> samples, quantiles, XiBar (N); 
                                samples.reserve(N); quantiles.reserve(2*N+1);
                                stat_.getSamples(spq->tableName, attrCon->attr, idx->at(i), samples);
                                stat_.getQuantiles(spq->tableName, attrCon->attr, idx->at(i), quantiles);
                                for (size_t j = 0; j < N; ++j) XiBar[j] = stoXs[stoInd][j] - sol.at(i)*samples[j];
                                KDE kdeXiBar (XiBar, true);
                                dF[i] = kdeXiBar.convolve(quantiles, v, sol.at(i));
                            }
                            stoCon[i] = dF[i];
                        }
                    }
                }
                stoBound = p;
                if (probCon->vsign == Inequality::gteq) stoBound = 1-p;
                double dot = 0;
                for (const auto& p : sol) dot += p.second*dF[p.first];
                stoBound += dot-Fxvs[stoInd];
                if (probCon->vsign != probCon->psign) stoIneq = Inequality::gteq;
                STOP("sto");
            }
            CLOCK("stocon");
            if (multiplier > 0){
                if (!sameSense(stoIneq, objSense)) multiplier *= -1;
                double stoNorm = norm(stoCon);
                for (size_t i = 0; i < nInds; ++i) softObj[i] += multiplier*stoCon[i]/stoNorm;
            } else{
                auto preNzInd = nzInd;
                nzInd += nInds + (getCvar(con) != nullptr);
                cbeg.push_back(nzInd);
                cind.resize(nzInd); cval.resize(nzInd);
                copy(idx->inds.begin(), idx->inds.end(), cind.begin()+preNzInd);
                copy(stoCon.begin(), stoCon.end(), cval.begin()+preNzInd);
                rhs.push_back(stoBound);
                if (stoIneq == Inequality::lteq) sense.push_back(GRB_LESS_EQUAL);
                else if (stoIneq == Inequality::gteq) sense.push_back(GRB_GREATER_EQUAL);
                rowInd ++;
            }
            stoInd ++;
            STOP("stocon");
        }
        if (isDeterministic(con, boundCon)){
            CLOCK("det");
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
                    copy(idx->inds.begin(), idx->inds.end(), cind.begin()+preNzInd);
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
                    copy(idx->inds.begin(), idx->inds.end(), cind.begin()+preNzInd);
                    copy(detCons[detInd].begin(), detCons[detInd].end(), cval.begin()+preNzInd);
                    rhs.push_back(ubValue);
                    sense.push_back(GRB_LESS_EQUAL);
                    rowInd ++;
                }
            }
            detInd ++;
            STOP("det");
        }
    }
    STOP("posvio");
    CLOCK("pregu");
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
    vector<char> vtype (numvars, GRB_CONTINUOUS);
	ckgb(GRBnewmodel(env, &model, NULL, numvars, softObj.data(), lb.data(), ub.data(), vtype.data(), NULL), env, model);
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
            ckgb(GRBsetintattrarray(model, GRB_INT_ATTR_VBASIS, 0, numvars, vStart.data()), env, model);
            ckgb(GRBsetintattrarray(model, GRB_INT_ATTR_CBASIS, 0, rowInd, cStart.data()), env, model);
            ckgb(GRBsetdblattrarray(model, GRB_DBL_ATTR_PSTART, 0, numvars, pStart.data()), env, model);
            ckgb(GRBsetdblattrarray(model, GRB_DBL_ATTR_DSTART, 0, rowInd, dStart.data()), env, model);
        }
    }
    preViolations = violations;
    STOP("pregu");
    CLOCK("gu");
	ckgb(GRBoptimize(model), env, model);
    STOP("gu");
    CLOCK("posgu");
	int status;
	ckgb(GRBgetintattr(model, GRB_INT_ATTR_STATUS, &status), env, model);
	if (status == GRB_OPTIMAL){
        pStart.resize(numvars);
    	ckgb(GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, numvars, pStart.data()), env, model);
        for (size_t i = 0; i < numvars; ++i) {
            if (pStart[i] > 0){
                nextSol[i] = pStart[i];
            }
        }
        double solutionSize = 0;
        for (const auto& p : nextSol){
            if (p.first < nInds) solutionSize += p.second;
        }
        maxSolutionSize = max(maxSolutionSize, solutionSize);
        gamma *= 1-1/maxSolutionSize;

        vStart.resize(numvars); cStart.resize(rowInd);
        ckgb(GRBgetintattrarray(model, GRB_INT_ATTR_VBASIS, 0, numvars, vStart.data()), env, model);
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
    STOP("posgu");
}

void Taylor::update(const SolIndType& step){
    CLOCK("c1");
    assert(step.size());
    auto n = idx->size();
    int m = spq->cons.size();
    vector<size_t> inds; inds.reserve(step.size());
    for (const auto& p : step) inds.push_back(p.first);
    vector<size_t> stoInds (m);
    size_t stoInd = 0, detInd = 0;
    for (size_t j = 0; j < m; ++j){
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
        #pragma omp for schedule(dynamic) collapse(2)
        for (size_t j = 0; j < m; ++j){
            for (size_t i = 0; i < inds.size(); ++i){
                const auto& con = spq->cons[j];
                if (isStochastic(con, probCon, attrCon) && attrCon){
                    size_t N = stat_.pg->getColumnLength(spq->tableName, attrCon->attr);
                    long long id = idx->at(inds[i]);
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
    STOP("c1");
    if (isDependentVar){
        CLOCK("c2");
        stoInd = 0;
        double nMaxSampleChanges = budget / (n*(spq->nStochastic-spq->nCvar));
        shared_ptr<ProbConstraint> probCon;
        shared_ptr<AttrConstraint> attrCon;
        vector<vector<double>> mulChanges (m);
        vector<vector<size_t>> indexChanges (m);
        for (const auto& con : spq->cons){
            if (isStochastic(con, probCon, attrCon) && attrCon){
                if (getVar(con)){
                    double v = spq->getValue(probCon->v);
                    AccSet acc;
                    for (auto x : stoXs[stoInd]) acc(x);
                    varXs[stoInd] = ba::variance(acc);
                    auto N = stoXs[stoInd].size();
                    double range = getKernelRange({nsix[stoInd], varXs[stoInd]});
                    vector<pair<double, size_t>> coefs (N);
                    for (size_t j = 0; j < N; ++j) coefs[j] = {-fabs(v-stoXs[stoInd][j])/range, j};
                    std::sort(coefs.begin(), coefs.end());
                    int cnt = 0;
                    for (const auto& pr : coefs){
                        double coef = -pr.first;
                        size_t ind = pr.second;
                        if (coef <= 1){
                            if (!jvInds[stoInd].count(ind) && (cnt++ < nMaxSampleChanges)){
                                indexChanges[stoInd].push_back(ind);
                                mulChanges[stoInd].push_back(1.0);
                                jvInds[stoInd].insert(ind);
                            }
                        } else{
                            if (jvInds[stoInd].count(ind)){
                                indexChanges[stoInd].push_back(ind);
                                mulChanges[stoInd].push_back(-1.0);
                                jvInds[stoInd].erase(ind);
                            }
                        }
                    }
                    deb(indexChanges[stoInd], mulChanges[stoInd]);
                }
                stoInd ++;
            }
        }
        STOP("c2");
        CLOCK("c3");
        #pragma omp parallel num_threads(nCores)
        {
            auto coreIndex = omp_get_thread_num();
            Stat stat_;
            shared_ptr<ProbConstraint> probCon;
            shared_ptr<AttrConstraint> attrCon;
            size_t stoInd = 0;
            const auto& interval = idx->getInterval(coreIndex);
            size_t n_ = interval.second - interval.first;
            for (const auto& con : spq->cons){
                if (isStochastic(con, probCon, attrCon) && attrCon){
                    auto sampleSize = indexChanges[stoInd].size();
                    if (getVar(con) && sampleSize){
                        auto N = stoXs[stoInd].size();
                        vector<vector<double>> samples; samples.reserve(n_);
                        stat_.getSamples(spq->tableName, attrCon->attr, idx->getSql(coreIndex), indexChanges[stoInd], samples);
                        for (size_t i = interval.first; i < interval.second; ++i){
                            size_t ind = i-interval.first;
                            double a = -0.5*getKernelRange({nsix[stoInd], stoVars[stoInd][i]});
                            for (size_t j = 0; j < sampleSize; ++j) zeroPds[stoInd][i] += mulChanges[stoInd][j]*samples[ind][j]*a;
                        }
                    }
                    stoInd ++;
                }
            }
        }
        STOP("c3");
    }
    CLOCK("c4");
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
    STOP("c4");
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

void Taylor::doAdjustment(){
    for (const auto& pr : adjustments) spq->setVariable(pr.first, pr.second.second);
}

void Taylor::undoAdjustment(){
    for (const auto& pr : adjustments) spq->setVariable(pr.first, pr.second.first);
}

void Taylor::solve(){
    unique_ptr<Optim> optim = make_unique<Direct>();
    doAdjustment();
    for (int i = 0; i < nMaxIters; ++i){
        CLOCK("b");
        SolIndType nextSol; solve(nextSol);
        deb(nextSol);
        STOP("b");
        CLOCK("c");
        update(optim->towards(nextSol));
        STOP("c");
        PRINT(pro);
        if (status != TaylorStatus::not_yet_found) break;
    }
    undoAdjustment();
    if (status == TaylorStatus::not_yet_found) status = minVio == 0 ? TaylorStatus::found : TaylorStatus::not_found;
}

SolType Taylor::getSol(const SolIndType& sol) const{
    SolType res;
    auto n = idx->size();
    if (sol.size()){
        for (const auto& p : sol){
            if (p.first < n && p.first >= 0) res[idx->at(p.first)] = p.second;
        }
    } else{
        for (const auto& p : bestSol){
            if (p.first < n && p.first >= 0) res[idx->at(p.first)] = p.second;
        }
    }
    return res;
}