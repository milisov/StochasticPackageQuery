#ifndef MAPOP_HPP
#define MAPOP_HPP

#include <map>
#include <cmath>
#include <algorithm>

using std::map;
using std::max;

template<typename K, typename V>
map<K, V> sqrt(const map<K, V>& M1){
    map<K, V> res;
    for (const auto& p1 : M1) res[p1.first] = sqrt(p1.second);
    return res;
}

template<typename K, typename V>
V sum(const map<K, V>& M1){
    V res = 0;
    for (const auto& p1 : M1) res += p1.second;
    return res;
}

template<typename K, typename V>
map<K, V> add(const map<K, V>& M1, const V& c){
    map<K, V> res;
    for (const auto& p1 : M1) res[p1.first] = p1.second+c;
    return res;
}

template<typename K, typename V>
map<K, V> max(const map<K, V>& M1, const V& c){
    map<K, V> res;
    for (const auto& p1 : M1) res[p1.first] = max(p1.second, c);
    return res;
}

template<typename K, typename V>
map<K, V> abs(const map<K, V>& M1){
    map<K, V> res;
    for (const auto& p1 : M1) res[p1.first] = abs(p1.second);
    return res;
}

template<typename K, typename V>
map<K, V> add(const V& v1, const map<K, V>& M1, const V& v2, const map<K, V>& M2){
    map<K, V> res;
    for (const auto& p1 : M1) res[p1.first] = p1.second*v1;
    for (const auto& p2 : M2) res[p2.first] += p2.second*v2;
    return res;
}

template<typename K, typename V>
map<K, V> mul(const map<K, V>& M1, const V& c){
    map<K, V> res;
    for (const auto& p1 : M1) res[p1.first] = p1.second*c;
    return res;
}

template<typename K, typename V>
map<K, V> mul(const map<K, V>& M1, const map<K, V>& M2){
    map<K, V> res;
    for (const auto& p1 : M1){
        if (M2.count(p1.first)){
            res[p1.first] = p1.second*M2.at(p1.first);
        }
    }
    return res;
}

template<typename K, typename V>
map<K, V> div(const map<K, V>& M1, const map<K, V>& M2){
    map<K, V> res;
    for (const auto& p1 : M1){
        if (M2.count(p1.first)){
            res[p1.first] = p1.second/M2.at(p1.first);
        }
    }
    return res;
}

#endif