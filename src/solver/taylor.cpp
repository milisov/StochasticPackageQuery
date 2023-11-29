#include <Highs.h>
#include <iostream>
#include <fmt/core.h>
#include <omp.h>
#include <numeric>
#include <algorithm>

#include "taylor.hpp"
#include "core/kde.hpp"
#include "util/udebug.hpp"
#include "util/uconfig.hpp"

using std::cerr;
using std::dynamic_pointer_cast;
using std::min;
using std::to_string;
using std::make_unique;
using std::make_shared;
using std::fill;
using std::copy;

map<string, Option> defaultOptions = {
    {"number_of_cores", static_cast<int>(Config::getInstance()->nPhysicalCores)}, 
    {"soft_deterministic_constraint", false}, 
    {"dependency_var", false}, 
    {"time_limit", 60.0}
};

Taylor::Taylor(shared_ptr<StochasticPackageQuery> spq, const vector<long long>& ids, const map<string, Option>& options): spq(spq), ids(ids), options(options){
    stat = make_unique<Stat>();
    for (const auto& p : defaultOptions){
        if (!options.count(p.first)) this->options[p.first] = p.second;
    }
    // if (!options.count("number_of_scenarios")){
    //     for (const auto& con : spq->cons){
    //         shared_ptr<ProbConstraint> probCon = dynamic_pointer_cast<ProbConstraint>(con);
	// 	    if (probCon){
    //             shared_ptr<AttrConstraint> attrCon = dynamic_pointer_cast<AttrConstraint>(con);
	// 	        if (attrCon){
    //                 if (!this->options.count("number_of_scenarios")){
    //                     this->options["number_of_scenarios"] = pg->getColumnLength(spq->tableName, attrCon->attr);
    //                 } else this->options["number_of_scenarios"] = min(boost::get<int>(this->options["number_of_scenarios"]), pg->getColumnLength(spq->tableName, attrCon->attr));
    //             }
    //         }
    //     }
    // }
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
        std::iota(this->ids.begin(), this->ids.end(), 1LL);
    }
    size_t n = this->ids.size();
    sqn = sqrt(n);
    int nCores = boost::get<int>(this->options["number_of_cores"]);
    for (const auto& con : spq->cons){
        shared_ptr<ProbConstraint> probCon = dynamic_pointer_cast<ProbConstraint>(con);
        if (probCon){
            shared_ptr<AttrConstraint> attrCon = dynamic_pointer_cast<AttrConstraint>(con);
            if (attrCon){
                size_t N = stat->pg->getColumnLength(spq->tableName, attrCon->attr);
                stoXs.emplace_back(N, 0);
            }
        } else{
            shared_ptr<AttrConstraint> attrCon = dynamic_pointer_cast<AttrConstraint>(con);
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
                        vector<double> attrs_; attrs.reserve(n_);
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
            } else{
                shared_ptr<CountConstraint> countCon = dynamic_pointer_cast<CountConstraint>(con);
                if (countCon){
                    detCons.emplace_back(n, 1.0);
                    detNorms.push_back(sqn);
                }
            }
        }
    }
    obj.resize(n);
    fill(obj.begin(), obj.end(), 0);
    if (spq->obj){
        shared_ptr<AttrObjective> attrObj = dynamic_pointer_cast<AttrObjective>(spq->obj);
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
                    vector<double> attrs_; attrs.reserve(n_);
                    Stat stat_;
                    stat_.getDetAttrs(spq->tableName, attrObj->obj, betweenId, attrs_);
                    copy(attrs_.begin(), attrs_.end(), attrs.begin()+intervals[coreIndex]-1);
                }
            } else {
                attrs.reserve(n);
                stat->getDetAttrs(spq->tableName, attrObj->obj, joinId, attrs);
            }
            copy(attrs.begin(), attrs.end(), obj.begin());
        } else{
            if (dynamic_pointer_cast<CountObjective>(spq->obj)){
                fill(obj.begin(), obj.end(), 1.0);
            }
        }
    }
    normalize(obj);
    for (size_t i = 0; i < obj.size(); ++i) obj[i] = obj[i]*sqn;
}

void Taylor::solve(map<int, double>& nextSols) const{
    size_t n = ids.size();
    size_t m = spq->cons.size();
    size_t nzN = sol.size();
    bool isSoftDetCon = boost::get<bool>(options.at("soft_deterministic_constraint"));
    HighsModel model;
    model.lp_.sense_ = ObjSense::kMinimize;
    model.lp_.a_matrix_.format_ = MatrixFormat::kRowwise;
    model.lp_.num_col_ = n;
    model.lp_.col_lower_ = vector<double>(n, 0);
    model.lp_.col_upper_ = vector<double>(n, POS_INF);
    if (spq->repeat != StochasticPackageQuery::NO_REPEAT) model.lp_.col_upper_ = vector<double>(n, spq->repeat+1.0);
    vector<int> start_, index_;
    vector<double> violations (m, 0), Fxvs, detMuls, value_;
    start_.reserve(m+1); index_.reserve((m+1)*model.lp_.num_col_); value_.reserve((m+1)*model.lp_.num_col_);
    start_.push_back(0);
    vector<vector<double>> stoCons (spq->countStochastic());
    vector<shared_ptr<KDE>> kdes;
    vector<double> softObj = obj;
    if ((model.lp_.sense_ == ObjSense::kMinimize && spq->obj->objSense == ObjectiveSense::maximize)
    || (model.lp_.sense_ == ObjSense::kMaximize && spq->obj->objSense == ObjectiveSense::minimize)){
        for (size_t i = 0; i < obj.size(); ++i) softObj[i] = -obj[i];
    }
    size_t detInd = 0, stoInd = 0;
    for (const auto& con : spq->cons){
        shared_ptr<ProbConstraint> probCon = dynamic_pointer_cast<ProbConstraint>(con);
        if (probCon){
            auto kde = make_shared<KDE>(stoXs[stoInd], true);
            if (dynamic_pointer_cast<VarConstraint>(con)){
                Fxvs.push_back(kde->getQuickCdf(spq->getValue(probCon->v)));
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
        } else{
            shared_ptr<BoundConstraint> boundCon = dynamic_pointer_cast<BoundConstraint>(con);
            if (boundCon){
                double res = detXs[detInd];
                auto lb = spq->getValue(boundCon->lb);
                if (lb != NEG_INF && lb > res){
                    violations[detInd+stoInd] = detViolate((lb-res)/detNorms[detInd]*sqn);
                    detMuls.push_back(-static_cast<double>(model.lp_.sense_));
                }
                auto ub = spq->getValue(boundCon->ub);
                if (ub != POS_INF && ub < res){
                    violations[detInd+stoInd] = detViolate((res-ub)/detNorms[detInd]*sqn);
                    detMuls.push_back(static_cast<double>(model.lp_.sense_));
                }
            }
            detInd ++;
        }
    }
    double logitSum = logLogit(violations);
    for (size_t i = 0; i < violations.size(); ++ i){
        if (violations[i] > 0) violations[i] = exp(logit(violations[i])-logitSum);
    }
    size_t varInd = 0, rowInd = 0, nzInd = 0;
    detInd = 0; stoInd = 0;
    for (const auto& con : spq->cons){
        double multiplier = violations[detInd+stoInd]*sqn;
        shared_ptr<ProbConstraint> probCon = dynamic_pointer_cast<ProbConstraint>(con);
        if (probCon){
            if (dynamic_pointer_cast<VarConstraint>(con)){

                varInd ++;
            }
            stoInd ++;
        } else{
            if (multiplier > 0 && isSoftDetCon){
                multiplier *= detMuls[detInd];
                for (size_t i = 0; i < n; ++i) softObj[i] += multiplier*detCons[detInd][i]/detNorms[detInd];
            } else{
                nzInd += n;
                start_.push_back(nzInd);
                for (size_t i = 0; i < n; ++i){
                    index_.push_back(i);
                    value_.push_back(detCons[detInd][i]);
                }
                shared_ptr<BoundConstraint> boundCon = dynamic_pointer_cast<BoundConstraint>(con);
                if (boundCon){
                    model.lp_.row_lower_.push_back(spq->getValue(boundCon->lb));
                    model.lp_.row_upper_.push_back(spq->getValue(boundCon->ub));
                    rowInd ++;
                }
            }
            detInd ++;
        }
    }
    model.lp_.num_row_ = rowInd;
    model.lp_.col_cost_ = softObj;
    model.lp_.a_matrix_.start_ = start_;
    model.lp_.a_matrix_.index_ = index_;
    model.lp_.a_matrix_.value_ = value_;
    deb(start_, index_.size(), value_.size());
    deb(violations);
}