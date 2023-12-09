#ifndef UNUMERIC_HPP
#define UNUMERIC_HPP

#include <vector>
#include <ostream>
#include <string>
#include <utility>
#include <cmath>
#include <algorithm>

#include <boost/multiprecision/gmp.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include "udeclare.hpp"
#include "uconfig.hpp"

using std::string;
using std::vector;
using std::ostream;
using std::pair;
using std::min;

namespace ba = boost::accumulators;
typedef ba::accumulator_set<long double, ba::stats<ba::tag::sum, ba::tag::variance>> AccSet;
using long_double = boost::multiprecision::number<boost::multiprecision::backends::gmp_float<1<<10>>;

const double MACHINE_EPS = 1e-12;
const double LOG_MACHINE_EPS = log(MACHINE_EPS);
const double NUMERIC_EPS = Config::getInstance()->pt.get<double>("parameters.numeric_eps");

size_t hashSol(const SolIndType& sol);
bool isEqual(const SolIndType& sol1, const SolIndType& sol2);

bool isEqual(const double& a, const double& b);
bool isLess(const double& a, const double& b);
bool isGreater(const double& a, const double& b);
bool isLessEqual(const double& a, const double& b);
bool isGreaterEqual(const double& a, const double& b);

template<typename T>
vector<pair<T, size_t>> topAbsK(const vector<T>& vec, const size_t& k){
    vector<pair<T, size_t>> tmpVec (vec.size());
    auto sz = min(vec.size(), k);
    for (size_t i = 0; i < vec.size(); ++i) tmpVec[i] = {abs(vec[i]), i};
    std::partial_sort(tmpVec.begin(), tmpVec.begin() + sz, tmpVec.end(), std::greater<pair<T, size_t>>());
    return vector<pair<T, size_t>>(tmpVec.begin(), tmpVec.begin() + sz);
}

template<typename T>
vector<pair<T, size_t>> topK(const vector<T>& vec, const size_t& k){
    vector<pair<T, size_t>> tmpVec (vec.size());
    auto sz = min(vec.size(), k);
    for (size_t i = 0; i < vec.size(); ++i) tmpVec[i] = {vec[i], i};
    std::partial_sort(tmpVec.begin(), tmpVec.begin() + sz, tmpVec.end(), std::greater<pair<T, size_t>>());
    return vector<pair<T, size_t>>(tmpVec.begin(), tmpVec.begin() + sz);
}

template<typename T>
vector<pair<T, size_t>> lowK(const vector<T>& vec, const size_t& k){
    vector<pair<T, size_t>> tmpVec (vec.size());
    auto sz = min(vec.size(), k);
    for (size_t i = 0; i < vec.size(); ++i) tmpVec[i] = {vec[i], i};
    std::partial_sort(tmpVec.begin(), tmpVec.begin() + sz, tmpVec.end());
    return vector<pair<T, size_t>>(tmpVec.begin(), tmpVec.begin() + sz);
}

template<typename T>
T sortedMedian(const vector<T>& arr, size_t st, size_t fn){
    auto len = fn - st + 1;
    if (len%2) return arr[st+len/2];
    len /= 2;
    return (arr[st+len] + arr[st+len-1]) / 2;
}

template<typename T>
T median(const vector<T>& arr){
    auto tmp = arr;
    std::sort(tmp.begin(), tmp.end());
    return sortedMedian(tmp, 0, arr.size()-1);
}

template<typename T>
T sortedIQR(const vector<T>& arr){
    size_t half = arr.size()/2;
    if (arr.size()%2) return sortedMedian(arr, half+1, arr.size()-1) - sortedMedian(arr, 0, half-1);
    return sortedMedian(arr, half, arr.size()-1) - sortedMedian(arr, 0, half-1);
}

string strf(const long_double& value);

/**
 * @brief Divide the interval [start, end] inclusively into div intervals
 * 
 */
vector<long long> divideInterval(const long long& start, const long long& end, const int& div);

class UniqueIndexer{
private:
    vector<string> sqlIds;
    vector<long long> sortedIds, intervals;
public:
    vector<int> inds;
    UniqueIndexer(const int& nCores, const long long& size, const vector<long long>& ids);
    string getSql(const int& coreIndex) const;
    pair<long long, long long> getInterval(const int& coreIndex) const;
    long long at(const size_t& ind) const;
    size_t size() const;
};

class Indexer{
private:
    vector<size_t> indexOfIds;
    size_t sz;
public:
    string sqlId;
    Indexer(const vector<long long>& ids);
    size_t crushedIndex(const size_t& idIndex) const;
    size_t crushedSize() const;
};

/**
 * @brief Compute sigmoid function with multiplicity k
 * 
 * @param x 
 * @param k 
 * @return double 
 */
double sigmoid(const double& x, const double& k=1.0);

double normalize(vector<double>& vec);
double norm(const vector<double>& vec);
double detViolate(const double& posDif);
double logit(const double& x);
double logLogit(const vector<double>& vios);



class AccAggregator{
private:
    long_double sum, M2;
    long long count;
public:
    AccAggregator();
    void add(const long double& sum, const long double& M2, const long long& count);
    long long getCount() const;
    long_double getSum() const;
    long_double getM2() const;
    double getMean() const;
    double getStd() const;
    double getVariance() const;
    operator string() const;
};

ostream& operator<<(ostream& os, const AccAggregator& agg);

#endif