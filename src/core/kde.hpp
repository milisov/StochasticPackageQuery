#ifndef KDE_HPP
#define KDE_HPP

#include <vector>
#include <utility>

using std::vector;
using std::pair;

class KDE{
private:
    long double sum, m2;
    double mean, variance;
    vector<double> sup, supPdf, cumCdf, cumMean, cumCvar;
private:
    void getSupports(const vector<double>& sortedArr, const double& h);
public:
    KDE(const vector<double>& arr, bool quickMode);
// QuickMode's API availablity
    long double getSum() const;
    long double getM2() const;
    double getMean() const;
    double getVariance() const;
    double getMin() const;
    double getMax() const;
    double getPdf(const double& x) const;
    double getQuickCdf(const double& x) const;
    double convolve(const vector<double>& quantiles, const double& v, const double& m) const;
    template <typename T> void getSortedQuantiles(vector<T>& sortedQuantiles) const;
    template <typename T> void getQuantiles(vector<T>& quantiles) const;
// NormalMode
    double getCdf(const double& x) const;
    double getCvarInv(const double& v) const;
};

/**
 * @brief Compute sorted quantiles quickly without offline cumCdf
 * 
 * @param sortedQuantile results will be in-place
 */
template<typename T>
void KDE::getSortedQuantiles(vector<T>& sortedQuantiles) const{
    double cdfV = 0, nextCdfV;
    size_t sIndex = 1, qIndex = 0, n = sortedQuantiles.size(), N = sup.size();
    while (qIndex < n){
        while (sIndex < N && (nextCdfV=cdfV+(sup[sIndex]-sup[sIndex-1])*supPdf[sIndex]) <= sortedQuantiles[qIndex]){
            cdfV = nextCdfV;
            sIndex ++;
        }
        if (sIndex == N){
            for (size_t i = qIndex; i<n; ++i) sortedQuantiles[i] = static_cast<T>(sup.back());
            break;
        }
        while (qIndex < n && sortedQuantiles[qIndex] < nextCdfV){
            sortedQuantiles[qIndex] = static_cast<T>((sortedQuantiles[qIndex]-cdfV)/supPdf[sIndex]+sup[sIndex-1]);
            qIndex ++;
        }
    }
}

template<typename T>
void KDE::getQuantiles(vector<T>& quantiles) const{
    auto n = quantiles.size();
    vector<pair<T, size_t>> sortedQuantiles (n);
    for (size_t i = 0; i < n; ++i){
        sortedQuantiles[i] = {quantiles[i], i};
    }
    sort(sortedQuantiles.begin(), sortedQuantiles.end());
    vector<T> sortedResults (n);
    for (size_t i = 0; i < n; ++i) sortedResults[i] = sortedQuantiles[i].first;
    getSortedQuantiles(sortedResults);
    for (size_t i = 0; i < n; ++i) quantiles[sortedQuantiles[i].second] = sortedResults[i];
}

#endif