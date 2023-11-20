#ifndef UIO_HPP
#define UIO_HPP

#include <libpq-fe.h>
#include <filesystem>
#include <memory>
#include <string>
#include <map>
#include <vector>

#include "udeclare.hpp"
#include "spq/spq.hpp"

namespace fs = std::filesystem;

using std::shared_ptr;
using std::unique_ptr;
using std::string;
using std::map;
using std::vector;

shared_ptr<StochasticPackageQuery> parseSpaqlFromFile(string filePath);
fs::path getProjectDir();

struct PGconnDeleter {
    void operator()(PGconn* conn) {
        if (conn) PQfinish(conn);
    }
};

using PGconnPtr = std::unique_ptr<PGconn, PGconnDeleter>;

class PgManager{
private:
    static string conninfo;
    static PGconnPtr conn;
    static string id, schema;
    static vector<vector<string>> typeGroups;
private:
    static string getConnInfo();
    static Column getColumn(string dataType);
public:
    static PGconnPtr getConn();
    static long long getTableSize(string tableName);
    static vector<string> getTables();
    static map<string, Column> getColumns(string tableName);
    static bool existTable(string tableName);
    static bool dropTable(string tableName);
};

#endif