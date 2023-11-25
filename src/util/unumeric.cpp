#include <cmath>
#include <fmt/core.h>
#include <boost/math/special_functions/logaddexp.hpp>

#include "udebug.hpp"
#include "unumeric.hpp"

string strf(const long_double& value) {
    return value.str();
}

vector<long long> divideInterval(long long start, long long end, int div){
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

