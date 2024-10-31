#include <iostream>
#include <boost/algorithm/string/join.hpp>

#include "util/unumeric.hpp"
#include "util/udebug.hpp"
#include "core/mapop.hpp"
#include "core/kde.hpp"
#include "checker.hpp"

using std::make_unique;
using std::cout;

Checker::Checker(){
    stat = make_unique<Stat>();
}

SPQChecker::SPQChecker(shared_ptr<StochasticPackageQuery> spq): Checker(), spq(spq){
    validateTableName = spq->tableName + "_validate";
}

double SPQChecker::getObjective(const SolType& sol) const{
    double res = 0;
    if (spq->obj){
        shared_ptr<AttrObjective> attrObj;
        if (isDeterministic(spq->obj, attrObj)){
            if (attrObj){
                auto n = sol.size();
                vector<double> attrs; attrs.reserve(n);
                vector<string> strIds; strIds.reserve(n);
                for (const auto& p : sol) strIds.push_back(to_string(p.first));
                stat->getDetAttrs(spq->tableName, attrObj->obj, fmt::format("{} IN ({}) ORDER BY {}", PgManager::id, boost::join(strIds, ","), PgManager::id), attrs);
                //stat->getDetAttrs(validateTableName, attrObj->obj, fmt::format("{} IN ({}) ORDER BY {}", PgManager::id, boost::join(strIds, ","), PgManager::id), attrs);
                size_t i = 0;
                for (const auto& p : sol) res += p.second*attrs[i++];
            }
            if (getCount(spq->obj)) res += sum(sol);
        }
    }
    return res;
}

double SPQChecker::getConIndicator(const SolType& sol, shared_ptr<Constraint> con) const{
    double res = 0;
    shared_ptr<ProbConstraint> probCon;
    shared_ptr<BoundConstraint> boundCon;
    shared_ptr<AttrConstraint> attrCon;
    if (isStochastic(con, probCon, attrCon) && attrCon){
        size_t N = stat->pg->getColumnLength(validateTableName, attrCon->attr);
        vector<double> X (N, 0);
        for (const auto& p : sol){
            vector<double> samples; samples.reserve(N);
            stat->getSamples(validateTableName, attrCon->attr, p.first, samples);
            for (size_t i = 0; i < N; ++i) X[i] += samples[i]*p.second;
        }
        KDE kde (X, true);
        double v = spq->getValue(probCon->v);
        if (getVar(con)){
            res = kde.getQuickCdf(v);
            if (probCon->vsign == Inequality::gteq) res = 1-res;
        }
    }
    if (isDeterministic(con, boundCon, attrCon)){
        if (attrCon){
            auto n = sol.size();
            vector<double> attrs; attrs.reserve(n);
            vector<string> strIds; strIds.reserve(n);
            for (const auto& p : sol) strIds.push_back(to_string(p.first));
            stat->getDetAttrs(spq->tableName, attrCon->attr, fmt::format("{} IN ({}) ORDER BY {}", PgManager::id, boost::join(strIds, ","), PgManager::id), attrs);
            size_t i = 0;
            for (const auto& p : sol) res += p.second*attrs[i++];
        }
        if (getCount(con)) res += sum(sol);
    }
    return res;
}

bool SPQChecker::feasible(const SolType& sol) const{
    for (const auto& p : sol) if (isLess(p.second, 0)) return false;
    if (spq->repeat != StochasticPackageQuery::NO_REPEAT){
        for (const auto& p : sol) if (isGreater(p.second, spq->repeat+1)) return false;
    }
    auto tableSize = stat->pg->getTableSize(spq->tableName);
    for (const auto& p : sol) if (p.first < 1 || p.first > tableSize) return false;
    for (const auto& con : spq->cons){
        if (con->isViolate({getConIndicator(sol, con)})) return false;
    }
    return true;
}

void SPQChecker::display(const SolType& sol) const{
    // deb(sol);
    for (const auto& p : sol) if (isLess(p.second, 0)){
        cout << fmt::format("sol[{}]={}{}{}<0\n", p.first, RED, p.second, RESET);
    }
    if (spq->repeat != StochasticPackageQuery::NO_REPEAT){
        for (const auto& p : sol) if (isGreater(p.second, spq->repeat+1)){
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
		if (con) strCons.push_back("\t" + spq->substitute(con, {getConIndicator(sol, con)}));
	}
	cout << boost::join(strCons, " AND\n") << '\n';
    if (spq->obj) cout << spq->obj->toStr({getObjective(sol)}) + '\n';
}