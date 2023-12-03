#ifndef UIO_HPP
#define UIO_HPP

#include <libpq-fe.h>
#include <memory>
#include <string>
#include <map>
#include <vector>
#include <type_traits>
#include <stdexcept>
#include <boost/multiprecision/gmp.hpp>
#include <boost/algorithm/string.hpp>
#include <cstdlib>

#include "udeclare.hpp"
#include "udebug.hpp"

using std::to_string;
using std::shared_ptr;
using std::unique_ptr;
using std::string;
using std::map;
using std::vector;
using std::is_same;

string strf(const float& value);
string strf(const double& value);
string strf(const long double& value);

struct PGconnDeleter {
    void operator()(PGconn* conn) {
        if (conn) PQfinish(conn);
    }
};

using PGconnPtr = std::unique_ptr<PGconn, PGconnDeleter>;

void check(PGconnPtr& conn, PGresult* res, const char* file, const int& line);
void check(PGconnPtr& conn, bool failed, const char* file, const int& line);

#define ck(conn, arg) check(conn, arg, __FILE__, __LINE__)

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
    int getColumnLength(const string& tableName, const string& columnName);
    void addColumn(const string& tableName, const string& columnName, const string& pgType);
    bool existTable(const string& tableName);
    void dropTable(const string& tableName);
};

void readArray(char* start, vector<double>& array);

class SingleRow{
private:
    unique_ptr<PgManager> pg;
    PGresult* res;
public:
    ~SingleRow();
    SingleRow(const string& query);
    bool fetchRow();
    long long getBigInt(const int& columnIndex);
    double getNumeric(const int& columnIndex);
    void getArray(const int& columnIndex, vector<double>& result);
};

template<typename T>
string to_string(const vector<T>& arr){
	vector<string> strArr (arr.size());
	for (size_t i = 0; i < arr.size(); ++i){
		strArr[i] = strf(arr[i]);
	}
	return fmt::format("{{{0}}}", boost::join(strArr, ","));
}

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
    void set(const int& paramIndex, const T& v){
        string strV = to_string(v);
        vals[paramIndex] = new char[strV.size()+1];
        strcpy(vals[paramIndex], strV.c_str());
    }
    void send();
};

class BulkCopy{
private:
    size_t maxSize;
    string tableName, data;
    char delimiter;
    unique_ptr<PgManager> pg;
public:
    ~BulkCopy();
    BulkCopy(const string& tableName, const char& delimiter='|');
    template <typename T>
    void add(const T& v){
        string strV = to_string(v);
        if constexpr (is_same<T, string>::value){
            strV = "'" + strV + "'";
        }
        data += strV + delimiter;
    }
    void send();
};

#endif