#include <fstream>
#include <sstream>
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/algorithm/string.hpp>

#include "spq/parser.hpp"
#include "uio.hpp"
#include "udebug.hpp"

using std::ifstream;
using std::endl;
using std::cerr;
using std::stringstream;
using boost::algorithm::to_lower;

shared_ptr<StochasticPackageQuery> parseSpaqlFromFile(string filePath){
	fs::path spaqlPath (filePath);
	if (spaqlPath.is_relative()) spaqlPath = getProjectDir() / spaqlPath;
	ifstream in(spaqlPath);
	if (!in.is_open()) {
		cerr << "Could not open file " << spaqlPath << endl;
		return nullptr;
	}
	stringstream buffer;
	buffer << in.rdbuf();
	string spaqlQuery = buffer.str();
	in.close();
	return parseSpaql(spaqlQuery);
}

fs::path getProjectDir(){
	auto currentPath = fs::current_path();
	return currentPath.parent_path();
}

string PgManager::id, PgManager::schema;
vector<vector<string>> PgManager::typeGroups (toColumn.size()-1);

string PgManager::getConnInfo(){
	fs::path configPath = getProjectDir() / "config.cfg";
	boost::property_tree::ptree pt;
	boost::property_tree::ini_parser::read_ini(configPath.string(), pt);

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
	for (size_t i = 0; i < typeGroups.size(); ++i){
		for (auto token : typeGroups[i]){
			if (dataType.find(token) != string::npos)
				return static_cast<Column>(i);
		}
	}
	return Column::unsupported;
}

PGconnPtr PgManager::getConn(){
	return PGconnPtr(PQconnectdb(conninfo.c_str()), PGconnDeleter());
}

PGconnPtr PgManager::conn = PgManager::getConn();

/**
 * @brief Get the size of the table without throwing error
 * 
 * @param tableName table's name
 * @return long long the size of the table, 0 if the table neither not exists nor has no rows
 */
long long PgManager::getTableSize(string tableName){
	string sql = fmt::format("SELECT {} FROM \"{}\" ORDER BY {} DESC LIMIT 1", id, tableName, id);
	auto res = PQexec(conn.get(), sql.c_str());
	long long size = 0;
	if (PQntuples(res) > 0) size = atoll(PQgetvalue(res, 0, 0));
	PQclear(res);
	return size;
}

vector<string> PgManager::getTables(){
	string sql = "SELECT tablename FROM pg_catalog.pg_tables WHERE schemaname != 'pg_catalog' AND schemaname != 'information_schema'";
	auto res = PQexec(conn.get(), sql.c_str());
	vector<string> tables;
	for (int i = 0; i < PQntuples(res); ++i) tables.push_back(string(PQgetvalue(res, i, 0)));
	PQclear(res);
	return tables;
}

map<string, Column> PgManager::getColumns(string tableName){
	string sql = fmt::format(" \
		SELECT column_name, data_type FROM information_schema.columns \
		WHERE table_schema = '{}' AND table_name = '{}'", schema, tableName);
	auto res = PQexec(conn.get(), sql.c_str());
	map<string, Column> columns;
	for (int i = 0; i < PQntuples(res); i ++) columns[string(PQgetvalue(res, i, 0))] = getColumn(string(PQgetvalue(res, i, 1)));
	PQclear(res);
	return columns;
}

bool PgManager::existTable(string tableName){
	string sql = fmt::format("SELECT * FROM pg_tables WHERE tablename='{}' AND schemaname='{}'", tableName, schema);
	auto res = PQexec(conn.get(), sql.c_str());
	bool exists = PQntuples(res) > 0;
	PQclear(res);
	return exists;
}

/**
 * @brief Drop the table without throwing error
 * 
 * @param tableName table's name
 * @return true if the drop was successful
 * @return false if the table did not exist
 */
bool PgManager::dropTable(string tableName){
	string sql = fmt::format("DROP TABLE \"{}\"", tableName);
	auto res = PQexec(conn.get(), sql.c_str());
	bool success = PQresultStatus(res) == PGRES_TUPLES_OK;
	PQclear(res);
	return success;
}




