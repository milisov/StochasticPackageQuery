#include <boost/algorithm/string.hpp>
#include <fmt/core.h>
#include <libpq-fe.h>

#include "stat.hpp"
#include "util/uconfig.hpp"
#include "util/udebug.hpp"

using boost::algorithm::to_lower_copy;
using std::make_unique;

bool Stat::rebuildStat = to_lower_copy(Config::getInstance()->pt.get<string>("build.rebuild_stat")) == "true";
string Stat::statTable = Config::getInstance()->pt.get<string>("database.stat_table");

Stat::Stat(){
    pg = make_unique<PgManager>();
    if (!pg->existTable(statTable) || rebuildStat){
        rebuildStat = false;
        if (pg->existTable(statTable)) pg->dropTable(statTable);
        string sql = fmt::format("CREATE TABLE {} (\
            table_name VARCHAR(50), \
            column_name VARCHAR(50), \
            mean REAL, \
            std REAL, \
            count BIGINT)", statTable);
        auto res = PQexec(pg->conn.get(), sql.c_str());
        PQclear(res);
    }
}

