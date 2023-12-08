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
using std::lower_bound;

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

KDE::KDE(const vector<double>& arr, const bool& quickMode){
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

double KDE::getPdf(const double& x) const{
    auto it = lower_bound(sup.begin(), sup.end(), x);
    if (it == sup.end()) return 0;
    return supPdf[it-sup.begin()];
}

double KDE::getQuickCdf(const double& x) const{
    if (x < sup[0]) return 0;
    double res = 0;
    size_t ind = 1;
    while (ind < sup.size() && x >= sup[ind]){
        res += supPdf[ind]*(sup[ind]-sup[ind-1]);
        ind ++;
    }
    if (ind == sup.size()) return 1;
    return res + supPdf[ind]*(x-sup[ind-1]);
}

double KDE::convolve(const vector<double>& quantiles, const double& v, const double& m) const{
    size_t n = sup.size();
    size_t N = quantiles.size();
    double res = 0, curV = 0, qV = 0;
    size_t vxmInd = 0, qInd = 0;
    double vxm = (v-sup.back()) / m;
    double step = 0.5/(N-1);
    while (1){
        if (vxm < quantiles[qInd]){
            double nextV = qV;
            if (qInd > 0) nextV += (vxm+quantiles[qInd-1])*step;
            if (vxmInd > 0) res += (nextV - curV)*supPdf[n-vxmInd];
            curV = nextV;
            vxmInd ++;
            if (vxmInd == n) break;
            vxm = (v-sup[n-vxmInd-1]) / m;
        } else{
            if (qInd > 0) qV += (quantiles[qInd]+quantiles[qInd-1])*step;
            if (vxmInd > 0) res += (qV - curV)*supPdf[n-vxmInd];
            curV = qV;
            qInd ++;
            if (qInd == N) break;
        }
    }
    return -res;
}

double KDE::getQuickCvarInv(const double& x) const{
    return 0;
}

// Normal mode

// double KDE::getCdf(const double& x) const{
//     auto it = lower_bound(sup.begin(), sup.end(), x);
//     auto ind = it-sup.begin();
//     if (ind == 0) return 0;
//     if (ind == sup.size()) return 1;
//     return cumCdf[ind-1] + supPdf[ind-1]*(x-sup[ind-1]);
// }