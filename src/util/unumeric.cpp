#include <fmt/core.h>
#include <algorithm>
#include <map>
#include <boost/math/special_functions/logsumexp.hpp>
#include <boost/functional/hash.hpp>
#include <boost/algorithm/string.hpp>

#include "unumeric.hpp"
#include "udebug.hpp"
#include "uio.hpp"

using std::iota;
using std::to_string;
using std::map;

size_t hashSol(const SolIndType& sol){
    size_t seed = 0;
    for (const auto& p : sol){
        boost::hash_combine(seed, p.first);
        boost::hash_combine(seed, static_cast<long long>(p.second/NUMERIC_EPS));
    }
    return seed;
}

bool isEqual(const SolIndType& sol1, const SolIndType& sol2){
    for (const auto& p : sol1){
        if (!sol2.count(p.first)) return false;
        if (abs(p.second-sol2.at(p.first)) > MACHINE_EPS) return false;
    }
    for (const auto& p : sol2){
        if (!sol1.count(p.first)) return false;
    }
    return true;
}

bool isEqual(const double& a, const double& b){
    return fabs(a-b) < NUMERIC_EPS;
}

bool isLess(const double& a, const double& b){
    return a<b && !isEqual(a, b);
}

bool isGreater(const double& a, const double& b){
    return a>b && !isEqual(a, b);
}

bool isLessEqual(const double& a, const double& b){
    return a<b || isEqual(a, b);
}

bool isGreaterEqual(const double& a, const double& b){
    return a>b || isEqual(a, b);
}

string strf(const long_double& value) {
    return value.str();
}

vector<long long> divideInterval(const long long& start, const long long& end, const int& div){
    vector<long long> result;
    long long intervalLength = end - start + 1;
    long long size = intervalLength / div;
    long long remainder = intervalLength % div;
    long long current = start;
    for (int i = 0; i < div; ++i) {
        result.push_back(current);
        current += size + (i < remainder ? 1 : 0);
    }
    result.push_back(end + 1);
    return result;
}

UniqueIndexer::UniqueIndexer(const int& nCores, const long long& size, const vector<long long>& ids){
    if (ids.size()){
        auto tmp = ids;
        std::sort(tmp.begin(), tmp.end());
        sortedIds.reserve(tmp.size());
        for (size_t i = 0; i < tmp.size()-1; ++i){
            if (tmp[i] != tmp[i+1]) sortedIds.push_back(tmp[i]);
        }
        sortedIds.push_back(tmp.back());
    } else{
        sortedIds.resize(size);
        iota(sortedIds.begin(), sortedIds.end(), 1LL);
    }
    intervals = divideInterval(0, sortedIds.size()-1, nCores);
    sqlIds.resize(nCores);
    for (int i = 0; i < nCores; ++i){
        if (ids.size()){
            vector<string> strIds (intervals[i+1]-intervals[i]);
            for (long long j = intervals[i]+1; j <= intervals[i+1]; ++j) strIds[j] = to_string(j);
            sqlIds[i] = fmt::format("{} IN ({}) ORDER BY {}", PgManager::id, boost::join(strIds, ","), PgManager::id);
        } else{
            sqlIds[i] = fmt::format("{} BETWEEN {} AND {} ORDER BY {}", PgManager::id, intervals[i]+1, intervals[i+1], PgManager::id);
        }
    }
    inds.resize(sortedIds.size());
    iota(inds.begin(), inds.end(), 0);
}

string UniqueIndexer::getSql(const int& coreIndex) const{
    return sqlIds.at(coreIndex);
}

pair<long long, long long> UniqueIndexer::getInterval(const int& coreIndex) const{
    return {intervals[coreIndex], intervals[coreIndex+1]};
}

long long UniqueIndexer::at(const size_t& ind) const{
    return sortedIds.at(ind);
}

size_t UniqueIndexer::size() const{
    return sortedIds.size();
}

Indexer::Indexer(const vector<long long>& ids){
    map<long long, size_t> sortedIds;
    for (const auto& id : ids) sortedIds.emplace(id, 0);
    sz = sortedIds.size();
    size_t i = 0;
    vector<string> strIds (sortedIds.size());
    for (auto& pr : sortedIds){
        strIds[i] = to_string(pr.first);
        pr.second = i++;
    }
    indexOfIds.resize(ids.size());
    for (size_t i = 0; i < ids.size(); ++i) indexOfIds[i] = sortedIds.at(ids.at(i));
    sqlId = fmt::format("{} IN ({}) ORDER BY {}", PgManager::id, boost::join(strIds, ","), PgManager::id);
}

size_t Indexer::crushedIndex(const size_t& idIndex) const{
    return indexOfIds.at(idIndex);
}

size_t Indexer::crushedSize() const{
    return sz;
}

double sigmoid(const double& x, const double& k){
    return exp(-boost::math::logaddexp(0.0, -k*x));
}

double normalize(vector<double>& vec){
    double norm = 0;
    for (const auto& v : vec) norm += v*v;
    norm = sqrt(norm);
    if (norm > 0) std::transform(vec.begin(), vec.end(), vec.begin(), [norm](const double& val) { return val/norm; });
    return norm;
}

double norm(const vector<double>& vec){
    double norm = 0;
    for (const auto& v : vec) norm += v*v;
    return sqrt(norm);
}

double detViolate(const double& posDif){
    return 1-1.0/(posDif+1);
}

double logit(const double& x){
    if (x == 1.0) return -LOG_MACHINE_EPS;
    if (x == 0.0) return LOG_MACHINE_EPS;
    return log(x)-log(1-x);
}

double logLogit(const vector<double>& vios){
    vector<double> posVios;
    for (const auto& v: vios) if (v > 0) posVios.push_back(logit(v));
    if (!posVios.size()) return 0;
    return boost::math::logsumexp(posVios);
}

AccAggregator::AccAggregator(){
    count = 0;
    sum = 0;
    M2 = 0;
}

void AccAggregator::add(const long double& sum, const long double& M2, const long long& count){
    if (this->count){
        long long nextCount = this->count + count;
        long_double delta = this->sum/this->count - sum/count;
        this->M2 += M2 + delta*delta*this->count/nextCount*count;
        this->sum += sum;
        this->count = nextCount;
    } else{
        this->count = count;
        this->sum = sum;
        this->M2 = M2;
    }
}

long long AccAggregator::getCount() const{
    return count;
}

long_double AccAggregator::getSum() const{
    return sum;
}

long_double AccAggregator::getM2() const{
    return M2;
}

double AccAggregator::getMean() const{
    if (count) return (sum / count).convert_to<double>();
    return 0;
}

double AccAggregator::getStd() const{
    return sqrt(getVariance());
}

double AccAggregator::getVariance() const{
    if (count) return (M2 / count).convert_to<double>();
    return 0;
}

AccAggregator::operator string() const{
    return fmt::format("Count:{} Mean:{} Std:{}", getCount(), getMean(), getStd());
}

ostream& operator<<(ostream& os, const AccAggregator& agg){
    os << static_cast<string>(agg);
    return os;
}

