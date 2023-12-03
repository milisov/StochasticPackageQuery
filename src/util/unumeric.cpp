#include <fmt/core.h>
#include <algorithm>
#include <boost/math/special_functions/logsumexp.hpp>

#include "udebug.hpp"
#include "unumeric.hpp"

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

