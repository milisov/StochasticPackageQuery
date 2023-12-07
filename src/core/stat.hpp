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
private:
    bool isAnalyzed(const string& tableName, const string& columnName);
    void addStat(const string& tableName, const string& columnName, const long_double& sum, const long_double& m2, const long long& count);
    void analyzeStochastic(const string& tableName, const vector<string>& columnNames);
    void analyzeDeterministic(const string& tableName, const vector<string>& columnNames);
public:
    unique_ptr<PgManager> pg;
    Stat();
    bool analyze(const string& tableName);
    void getStoMeanVars(const string& tableName, const string& columnName, const string& sqlId, vector<double>& means, vector<double>& vars);
    void getDetAttrs(const string& tableName, const string& columnName,const string& sqlId, vector<double>& attrs);
    pair<double, double> getDetMeanVar(const string& tableName, const string& columnName);
    void getSamples(const string& tableName, const string& columnName, const long long& id, vector<double>& samples);
    void getSamples(const string& tableName, const string& columnName, const string& sqlId, const vector<size_t>& sampleIds, vector<vector<double>>& samples);
    void getQuantiles(const string& tableName, const string& columnName, const long long& id, vector<double>& quantiles);
};

#endif