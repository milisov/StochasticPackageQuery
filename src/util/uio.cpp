#include <fmt/core.h>
#include <iostream>
#include <filesystem>
#include <algorithm>
#include <boost/algorithm/string.hpp>

#include "uconfig.hpp"
#include "uio.hpp"

using boost::algorithm::to_lower;
using std::cerr;
using std::make_unique;

string PgManager::id, PgManager::schema;
vector<vector<string>> PgManager::typeGroups (toColumn.size()-1);

string PgManager::getConnInfo(){
	auto pt = Config::getInstance()->pt;
	schema = pt.get<string>("postgres.schema");
	id = pt.get<string>("pgmanager.index_column");
	for (size_t i = 0; i < typeGroups.size(); ++i){
		string columnType = to_string(static_cast<Column>(i));
		boost::split(typeGroups[i], pt.get<string>(fmt::format("pgmanager.{}", columnType)), boost::is_any_of(","));
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
	auto status = PQresultStatus(res);
	if (status != PGRES_SINGLE_TUPLE && status != PGRES_TUPLES_OK && status != PGRES_COMMAND_OK && status != PGRES_PIPELINE_SYNC) {
		cerr << PQerrorMessage(conn.get());
		PQclear(res);
		conn.reset();
		exit(1);
	}
}

void ck(PGconnPtr& conn, bool failed){
	if (failed){
		cerr << PQerrorMessage(conn.get());
		conn.reset();
		exit(1);
	}
}

PgManager::PgManager(){
	conn = PGconnPtr(PQconnectdb(conninfo.c_str()), PGconnDeleter());
	ck(conn, PQstatus(conn.get()) != CONNECTION_OK);
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

SingleRow::~SingleRow(){
	while (res){
		PQclear(res);
		res = PQgetResult(pg->conn.get());
	}
}

SingleRow::SingleRow(const string& query){
	pg = make_unique<PgManager>();
	ck(pg->conn, !PQsendQuery(pg->conn.get(), query.c_str()));
	ck(pg->conn, !PQsetSingleRowMode(pg->conn.get()));
	res = nullptr;
}

bool SingleRow::fetchRow(){
	ck(pg->conn, !PQconsumeInput(pg->conn.get()));
	if (res) PQclear(res);
	res = PQgetResult(pg->conn.get());
	if (res){
		ck(pg->conn, res);
		if (PQntuples(res)) return true;
	}
	return false;
}

long long SingleRow::getBigInt(int columnIndex){
	return atoll(PQgetvalue(res, 0, columnIndex));
}

double SingleRow::getNumeric(int columnIndex){
	return atof(PQgetvalue(res, 0, columnIndex));
}


void SingleRow::getRealArray(int columnIndex, vector<float>& result){
	char* pEnd = PQgetvalue(res, 0, columnIndex);
	char* curPtr;
	while (*pEnd != '}'){
		curPtr = pEnd+1;
		result.push_back(strtof(curPtr, &pEnd));
	}
}

string to_string(const vector<float>& arr){
	vector<string> strArr (arr.size());
	for (size_t i = 0; i < arr.size(); ++i){
		strArr[i] = strf(arr[i]);
	}
	return fmt::format("{{{0}}}", boost::join(strArr, ","));
}

const string AsyncUpdate::stmName = "async";

AsyncUpdate::~AsyncUpdate(){
	delete[] vals;
	ck(pg->conn, !PQpipelineSync(pg->conn.get()));
	while (1){
		auto res = PQgetResult(pg->conn.get());
		ck(pg->conn, res);
		if (res) PQclear(res);
		if (PQresultStatus(res) == PGRES_PIPELINE_SYNC) break;
		res = PQgetResult(pg->conn.get());
	}
	ck(pg->conn, !PQexitPipelineMode(pg->conn.get()));
}

/**
 * @brief Construct a new Async Update:: Async Update object
 * 
 * @param stmName 
 * @param query should not contain $ for other purpose except for specifying parameters
 */
AsyncUpdate::AsyncUpdate(const string& query){
	pg = make_unique<PgManager>();
	nParams = std::count(query.begin(), query.end(), '$');
	auto res = PQprepare(pg->conn.get(), stmName.c_str(), query.c_str(), nParams, NULL);
	ck(pg->conn, res);
	ck(pg->conn, PQsetnonblocking(pg->conn.get(), 1) || !PQisnonblocking(pg->conn.get()));
	ck(pg->conn, !PQenterPipelineMode(pg->conn.get()) || PQpipelineStatus(pg->conn.get()) != PQ_PIPELINE_ON);
	vals = new char*[nParams];
}

void AsyncUpdate::send(){
	ck(pg->conn, !PQsendQueryPrepared(pg->conn.get(), stmName.c_str(), nParams, vals, NULL, NULL, 0));
	for (int i = 0; i < nParams; ++i) delete[] vals[i];
}







