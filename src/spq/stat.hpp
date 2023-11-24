#ifndef STAT_HPP
#define STAT_HPP

#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

#include "util/uio.hpp"
#include "util/unumeric.hpp"

using std::unique_ptr;
using std::string;
using std::unordered_map;
using std::pair;

class Stat{
private:
    static bool rebuildStat;
    static const string statTable;
    unique_ptr<PgManager> pg;
private:
    bool isAnalyzed(const string& tableName, const string& columnName);
    void addStat(const string& tableName, const string& columnName, const long_double sum, const long_double m2, const long long count);
    void analyzeStochastic(const string& tableName, const string& columnName);
    void analyzeDeterministic(const string& tableName, const vector<string>& columnNames);
public:
    Stat();
    bool analyze(const string& tableName);
    void getStochasticMeanVar(const string& tableName, const string& columnName, unordered_map<long long, pair<double, double>>& storages);
};

#endif