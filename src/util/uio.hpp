#ifndef UIO_HPP
#define UIO_HPP

#include <libpq-fe.h>
#include <memory>
#include <string>
#include <map>
#include <vector>

#include "udeclare.hpp"

using std::shared_ptr;
using std::unique_ptr;
using std::string;
using std::map;
using std::vector;

struct PGconnDeleter {
    void operator()(PGconn* conn) {
        if (conn) PQfinish(conn);
    }
};

using PGconnPtr = std::unique_ptr<PGconn, PGconnDeleter>;

class PgManager{
private:
    static string conninfo;
    static string id, schema;
    static vector<vector<string>> typeGroups;
private:
    static string getConnInfo();
    static Column getColumn(string dataType);
public:
    PGconnPtr conn;
public:
    PgManager();
    long long getTableSize(string tableName);
    vector<string> getTables();
    map<string, Column> getColumns(string tableName);
    bool existTable(string tableName);
    bool dropTable(string tableName);
};

#endif