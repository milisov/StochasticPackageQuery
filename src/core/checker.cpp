#include <iostream>
#include <boost/algorithm/string/join.hpp>

#include "util/udebug.hpp"
#include "core/kde.hpp"
#include "checker.hpp"

using std::make_unique;
using std::dynamic_pointer_cast;
using std::cout;

Checker::Checker(){
    stat = make_unique<Stat>();
}

SPQChecker::SPQChecker(shared_ptr<StochasticPackageQuery> spq): Checker(), spq(spq){
}

double SPQChecker::getObjective(const SolType& sol) const{
    double res = 0;
    if (spq->obj){
        shared_ptr<AttrObjective> attrObj = dynamic_pointer_cast<AttrObjective>(spq->obj);
        if (attrObj){
            auto n = sol.size();
            vector<double> attrs; attrs.reserve(n);
            vector<string> strIds; strIds.reserve(n);
            for (const auto& p : sol) strIds.push_back(to_string(p.first));
            stat->getDetAttrs(spq->tableName, attrObj->obj, boost::join(strIds, ","), attrs);
            size_t i = 0;
            for (const auto& p : sol) res += p.second*attrs[i++];
        } else{
            if (dynamic_pointer_cast<CountObjective>(spq->obj)){
                for (const auto& p : sol) res += p.second;
            }
        }
    }
    return res;
}

double SPQChecker::getConIndicator(const SolType& sol, shared_ptr<Constraint> con) const{
    double res = 0;
    shared_ptr<ProbConstraint> probCon = dynamic_pointer_cast<ProbConstraint>(con);
    if (probCon){
        shared_ptr<AttrConstraint> attrCon = dynamic_pointer_cast<AttrConstraint>(con);
        if (attrCon){
            size_t N = stat->pg->getColumnLength(spq->tableName, attrCon->attr);
            vector<double> X (N, 0);
            for (const auto& p : sol){
                vector<double> samples; samples.reserve(N);
                stat->getSamples(spq->tableName, attrCon->attr, p.first, samples);
                for (size_t i = 0; i < N; ++i) X[i] += samples[i]*p.second;
            }
            KDE kde (X, true);
            double v = spq->getValue(probCon->v);
            if (dynamic_pointer_cast<VarConstraint>(con)){
                res = kde.getQuickCdf(v);
                if (probCon->vsign == Inequality::gteq) res = 1-res;
            }
        }
    } else{
		shared_ptr<AttrConstraint> attrCon = dynamic_pointer_cast<AttrConstraint>(con);
		if (attrCon){
            auto n = sol.size();
            vector<double> attrs; attrs.reserve(n);
            vector<string> strIds; strIds.reserve(n);
            for (const auto& p : sol) strIds.push_back(to_string(p.first));
            stat->getDetAttrs(spq->tableName, attrCon->attr, boost::join(strIds, ","), attrs);
            size_t i = 0;
            for (const auto& p : sol) res += p.second*attrs[i++];
        } else{
            if (dynamic_pointer_cast<CountConstraint>(con)){
                for (const auto& p : sol) res += p.second;
            }
        }
    }
    return res;
}

bool SPQChecker::feasible(const SolType& sol) const{
    for (const auto& p : sol) if (p.second < 0) return false;
    if (spq->repeat != StochasticPackageQuery::NO_REPEAT){
        for (const auto& p : sol) if (p.second > spq->repeat+1) return false;
    }
    auto tableSize = stat->pg->getTableSize(spq->tableName);
    for (const auto& p : sol) if (p.first < 1 || p.first > tableSize) return false;
    for (const auto& con : spq->cons){
        double indicator = getConIndicator(sol, con);
    }
    return true;
}

void SPQChecker::display(const SolType& sol) const{
    deb(sol);
    for (const auto& p : sol) if (p.second < 0){
        cout << fmt::format("sol[{}]={}{}{}<0\n", p.first, RED, p.second, RESET);
    }
    if (spq->repeat != StochasticPackageQuery::NO_REPEAT){
        for (const auto& p : sol) if (p.second > spq->repeat+1){
            cout << fmt::format("sol[{}]={}{}{}>{}\n", p.first, RED, p.second, RESET, spq->repeat+1);
        }
    }
    auto tableSize = stat->pg->getTableSize(spq->tableName);
    for (const auto& p : sol) if (p.first < 1 || p.first > tableSize){
        cout << fmt::format("sol has index {}{}{} outside the range [{},{}] of table '{}'\n", RED, p.first, RESET, 1, tableSize, spq->tableName);
    }
	cout << fmt::format("SELECT PACKAGE({}) FROM {}", spq->strAttrList(), spq->tableName);
	if (spq->cons.size()) cout << " SUCH THAT\n";
	vector<string> strCons;
	for (const auto& con : spq->cons){
        double indicator = getConIndicator(sol, con);
		if (con) strCons.push_back("\t" + spq->substitute(con, {indicator}));
	}
	cout << boost::join(strCons, " AND\n") << '\n';
    if (spq->obj) cout << spq->obj->toStr({getObjective(sol)}) + '\n';
}