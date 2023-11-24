#ifndef KDE_HPP
#define KDE_HPP

#include <vector>

using std::vector;

class KDE{
private:
    long double sum, m2;
    double mean, variance;
    vector<double> sup, supPdf, cumCdf, cumMean, cumCvar;
private:
    void getSupports(const vector<double>& sortedArr, const double& h);
public:
    KDE(const vector<double>& arr, bool quickMode=false);
// QuickMode's API availablity
    long double getSum() const;
    long double getM2() const;
    double getMean() const;
    double getVariance() const;
    void getQuantiles(vector<float>& sortedQuantile) const;
// NormalMode
};

#endif