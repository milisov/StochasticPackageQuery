#include <fmt/core.h>
#include <filesystem>
#include <boost/algorithm/string.hpp>

#include "uconfig.hpp"
#include "uio.hpp"
#include "udebug.hpp"

using boost::algorithm::to_lower;
using std::runtime_error;

string PgManager::id, PgManager::schema;
vector<vector<string>> PgManager::typeGroups (toColumn.size()-1);

string PgManager::getConnInfo(){
	auto pt = Config::getInstance()->pt;
	schema = pt.get<string>("postgres.schema");
	id = pt.get<string>("database.index_column");
	for (size_t i = 0; i < typeGroups.size(); ++i){
		string columnType = to_string(static_cast<Column>(i));
		boost::split(typeGroups[i], pt.get<string>(fmt::format("database.{}", columnType)), boost::is_any_of(","));
	}

	return fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", 
		pt.get<string>("postgres.username"), 
		pt.get<string>("postgres.hostname"),
		pt.get<string>("postgres.port"), 
		pt.get<string>("postgres.database"), 
		pt.get<string>("postgres.password"));
}

string PgManager::conninfo = PgManager::getConnInfo();

Column PgManager::getColumn(string dataType){
	to_lower(dataType);
	for (size_t i = typeGroups.size()-1; i >= 0; --i){
		for (auto token : typeGroups[i]){
			Column column = to_Column(i);
			if (dataType.find(token) != string::npos){
				if (column == Column::array_type){
					bool isNumericColumn = false;
					for (auto numericToken : typeGroups[to_index(column)]){
						if (dataType.find(numericToken) != string::npos){
							isNumericColumn = true;
							break;
						}
					}
					if (!isNumericColumn) return Column::unsupported;
				}
				return column;
			}
		}
	}
	return Column::unsupported;
}

void ck(PGconnPtr& conn, PGresult* res){
	if (PQresultStatus(res) != PGRES_TUPLES_OK && PQresultStatus(res) != PGRES_COMMAND_OK) {
		throw runtime_error(PQerrorMessage(conn.get()));
	}
}

PgManager::PgManager(){
	conn = PGconnPtr(PQconnectdb(conninfo.c_str()), PGconnDeleter());
	if (PQstatus(conn.get()) != CONNECTION_OK) throw runtime_error(PQerrorMessage(conn.get()));
}

/**
 * @brief Get the size of the table
 * 
 * @param tableName table's name
 * @return long long the size of the table, 0 if the table neither not exists nor has no rows
 */
long long PgManager::getTableSize(const string& tableName){
	string sql = fmt::format("SELECT {} FROM \"{}\" ORDER BY {} DESC LIMIT 1", id, tableName, id);
	auto res = PQexec(conn.get(), sql.c_str());
	ck(conn, res);
	long long size = 0;
	if (PQntuples(res) > 0) size = atoll(PQgetvalue(res, 0, 0));
	PQclear(res);
	return size;
}

vector<string> PgManager::getTables(){
	string sql = "SELECT tablename FROM pg_catalog.pg_tables WHERE schemaname != 'pg_catalog' AND schemaname != 'information_schema'";
	auto res = PQexec(conn.get(), sql.c_str());
	ck(conn, res);
	vector<string> tables;
	for (int i = 0; i < PQntuples(res); ++i) tables.push_back(string(PQgetvalue(res, i, 0)));
	PQclear(res);
	return tables;
}

map<string, Column> PgManager::getColumns(const string& tableName){
	string sql = fmt::format(" \
		SELECT column_name, udt_name::regtype FROM information_schema.columns \
		WHERE table_schema = '{}' AND table_name = '{}'", schema, tableName);
	auto res = PQexec(conn.get(), sql.c_str());
	ck(conn, res);
	map<string, Column> columns;
	for (int i = 0; i < PQntuples(res); i ++) columns[string(PQgetvalue(res, i, 0))] = PgManager::getColumn(string(PQgetvalue(res, i, 1)));
	PQclear(res);
	return columns;
}

void PgManager::addColumn(const string& tableName, const string& columnName, const string& pgType){
	string sql = fmt::format("ALTER TABLE \"{}\" ADD COLUMN \"{}\" {}", tableName, columnName, pgType);
	auto res = PQexec(conn.get(), sql.c_str());
	ck(conn, res);
	PQclear(res);
}

bool PgManager::existTable(const string& tableName){
	string sql = fmt::format("SELECT * FROM pg_tables WHERE tablename='{}' AND schemaname='{}'", tableName, schema);
	auto res = PQexec(conn.get(), sql.c_str());
	ck(conn, res);
	bool exists = PQntuples(res) > 0;
	PQclear(res);
	return exists;
}

void PgManager::dropTable(const string& tableName){
	string sql = fmt::format("DROP TABLE \"{}\"", tableName);
	auto res = PQexec(conn.get(), sql.c_str());
	ck(conn, res);
	PQclear(res);
}




