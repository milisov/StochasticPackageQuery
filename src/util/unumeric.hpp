#ifndef UNUMERIC_HPP
#define UNUMERIC_HPP

#include <vector>
#include <ostream>
#include <string>

#include "udeclare.hpp"

using std::string;
using std::vector;
using std::ostream;

/**
 * @brief Divide the interval [start, end] inclusively into div intervals
 * 
 */
vector<long long> divideInterval(long long start, long long end, int div);

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