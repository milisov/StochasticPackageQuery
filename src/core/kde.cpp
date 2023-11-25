#include <algorithm>
#include <cmath>
#include <numeric>

#include "kde.hpp"
#include "util/unumeric.hpp"
#include "util/udebug.hpp"

using std::min;
using std::max;
using std::partial_sum;
using std::sort;

void KDE::getSupports(const vector<double>& sortedArr, const double& h){
    auto n = sortedArr.size();
    vector<double> st (n);
    vector<double> fn (n);
    for (size_t i = 0; i < n; ++i){
        st[i] = sortedArr[i]-h;
        fn[i] = sortedArr[i]+h;
    }
    sup.reserve(2*n);
    supPdf.reserve(2*n);
    supPdf.push_back(0);
    size_t stIndex=0, fnIndex=0;
    double sumV = 0, add=1/(2*h*n);
    while (stIndex < n){
        if (st[stIndex] < fn[fnIndex]){
            sumV += add;
            stIndex ++;
            if (stIndex == n || st[stIndex] > st[stIndex-1]){
                sup.push_back(st[stIndex-1]);
                supPdf.push_back(sumV);
            }
        } else if (st[stIndex] > fn[fnIndex]){
            sumV -= add;
            fnIndex ++;
            if (fn[fnIndex] > fn[fnIndex-1]){
                sup.push_back(fn[fnIndex-1]);
                supPdf.push_back(sumV);
            }
        } else{
            double v = st[stIndex];
            while (stIndex<n && st[stIndex]==fn[fnIndex]){
                stIndex ++;
                fnIndex ++;
            }
            if (stIndex == n || st[stIndex] > v){
                sup.push_back(v);
                supPdf.push_back(sumV);
            }
        }
    }
    while (fnIndex < n){
        double prev = fn[fnIndex];
        fnIndex ++;
        if (fnIndex<n && prev==fn[fnIndex]){
            sumV -= add;
            continue;
        }
        sumV -= add;
        sup.push_back(prev);
        supPdf.push_back(sumV);
    }
    supPdf.pop_back();
}

KDE::KDE(const vector<double>& arr, bool quickMode){
    AccSet acc;
    for (auto v : arr) acc(v);
    auto n = arr.size();
    sum = ba::sum(acc);
    variance = ba::variance(acc);
    m2 = variance*n;
    mean = sum/n;
    double std = sqrt(variance);

    vector<double> sortedArr = arr;
    sort(sortedArr.begin(), sortedArr.end());
    double h = max(0.9*min(std, sortedIQR(sortedArr)/1.34)*pow(static_cast<double>(n), -0.2), MACHINE_EPS);
    getSupports(sortedArr, h);

    if (!quickMode){
        auto N = sup.size();
        vector<double> base (N, 0);
        for (size_t i = 1; i < N; ++i) base[i] = (sup[i]-sup[i-1])*supPdf[i];
        cumCdf.resize(N);
        partial_sum(base.begin(), base.end(), cumCdf.begin());
        cumCdf.back() = 1;
        for (size_t i = 1; i < N; ++i) base[i] = 0.5*base[i]*(sup[i]+sup[i-1]);
        cumMean.resize(N);
        partial_sum(base.begin(), base.end(), cumMean.begin());
        cumMean.back() = mean;
        cumCvar.resize(N); cumCvar[0] = sup[0];
        for (size_t i = 1; i < N; ++i) cumCvar[i] = cumMean[i] / cumCdf[i];
    }
}

long double KDE::getSum() const{
    return sum;
}

long double KDE::getM2() const{
    return m2;
}

double KDE::getMean() const{
    return mean;
}

double KDE::getVariance() const{
    return variance;
}

double KDE::getMin() const{
    return sup.front();
}

double KDE::getMax() const{
    return sup.back();
}
