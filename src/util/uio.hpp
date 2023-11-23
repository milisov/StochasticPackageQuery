#ifndef UIO_HPP
#define UIO_HPP

#include <libpq-fe.h>
#include <memory>
#include <string>
#include <map>
#include <vector>
#include <type_traits>
#include <stdexcept>
#include <fmt/core.h>
#include <boost/multiprecision/gmp.hpp>

#include "udeclare.hpp"
#include "udebug.hpp"

namespace mp = boost::multiprecision;

using std::shared_ptr;
using std::unique_ptr;
using std::string;
using std::map;
using std::vector;
using std::is_same;

template <typename T>
string strf(const T& value);

template<>
inline string strf(const float& value) {
    return fmt::format("{:.7f}", value);
}

template<>
inline string strf(const double& value) {
    return fmt::format("{:.16f}", value);
}

template<>
inline string strf(const long double& value) {
    return fmt::format("{:.16Lf}", value);
}

template <>
inline string strf(const long_double& value) {
    return value.str();
}

struct PGconnDeleter {
    void operator()(PGconn* conn) {
        if (conn) PQfinish(conn);
    }
};

using PGconnPtr = std::unique_ptr<PGconn, PGconnDeleter>;

void ck(PGconnPtr& conn, PGresult* res);
void ck(PGconnPtr& conn, bool failed);

class PgManager{
private:
    static string conninfo;
    static vector<vector<string>> typeGroups;
private:
    static string getConnInfo();
    static Column getColumn(string dataType);
public:
    static string id, schema;
    PGconnPtr conn;
public:
    PgManager();
    long long getTableSize(const string& tableName);
    vector<string> getTables();
    map<string, Column> getColumns(const string& tableName);
    void addColumn(const string& tableName, const string& columnName, const string& pgType);
    bool existTable(const string& tableName);
    void dropTable(const string& tableName);
};

class SingleRow{
private:
    unique_ptr<PgManager> pg;
    PGresult* res;
public:
    ~SingleRow();
    SingleRow(const string& query);
    bool fetchRow();
    long long getBigInt(int columnIndex);
    double getNumeric(int columnIndex);
    void getRealArray(int columnIndex, vector<float>& result);
};

string to_string(const vector<float>& arr);

class AsyncUpdate{
private:
    static const string stmName;
    int nParams;
    char** vals;
    unique_ptr<PgManager> pg;
public:
    ~AsyncUpdate();
    AsyncUpdate(const string& query);
    template <typename T>
    void set(int paramIndex, const T& v){
        string strV = to_string(v);
        vals[paramIndex] = new char[strV.size()+1];
        strcpy(vals[paramIndex], strV.c_str());
    }
    void send();
};

#endif