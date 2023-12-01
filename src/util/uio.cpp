#include <fmt/core.h>
#include <iostream>
#include <filesystem>
#include <algorithm>

#include "uconfig.hpp"
#include "uio.hpp"

using boost::algorithm::to_lower;
using std::cerr;
using std::make_unique;

string strf(const float& value) {
    return fmt::format("{:.7f}", value);
}

string strf(const double& value) {
    return fmt::format("{:.16f}", value);
}

string strf(const long double& value) {
    return fmt::format("{:.16Lf}", value);
}

string PgManager::id, PgManager::schema;
vector<vector<string>> PgManager::typeGroups (toColumn.size()-1);

string PgManager::getConnInfo(){
	auto pt = Config::getInstance()->pt;
	schema = pt.get<string>("postgres.schema");
	id = pt.get<string>("pgmanager.index_column");
	for (size_t i = 0; i < typeGroups.size(); ++i){
		string columnType = str(static_cast<Column>(i));
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

const vector<ExecStatusType> okTypes = {
	PGRES_SINGLE_TUPLE,
	PGRES_TUPLES_OK,
	PGRES_COMMAND_OK,
	PGRES_PIPELINE_SYNC,
	PGRES_COPY_IN
};

void check(PGconnPtr& conn, PGresult* res, const char* file, int line){
	auto status = PQresultStatus(res);
	bool isBad = true;
	for (auto okType : okTypes){
		if (status == okType){
			isBad = false;
			break;
		}
	}
	if (isBad) {
		cerr << "\033[0;31m" << "File " << file  \
			<< ", Line " << line << "\033[0m" << "\n";
		cerr << PQerrorMessage(conn.get());
		PQclear(res);
		conn.reset();
		exit(1);
	}
}

void check(PGconnPtr& conn, bool failed, const char* file, int line){
	if (failed){
		cerr << "\033[0;31m" << "File " << file  \
			<< ", Line " << line << "\033[0m" << "\n";
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

int PgManager::getColumnLength(const string& tableName, const string& columnName){
	int len = 0;
	string sql = fmt::format(" \
		SELECT udt_name::regtype FROM information_schema.columns \
		WHERE table_schema = '{}' AND table_name = '{}' AND column_name='{}'", schema, tableName, columnName);
	auto res = PQexec(conn.get(), sql.c_str());
	ck(conn, res);
	auto columnType = PgManager::getColumn(string(PQgetvalue(res, 0, 0)));
	PQclear(res);
	if (columnType == Column::numeric_type || columnType == Column::string_type) len=1;
	else if (columnType == Column::array_type){
		sql = fmt::format("SELECT array_length({}, 1) FROM \"{}\" LIMIT 1", columnName, tableName);
		res = PQexec(conn.get(), sql.c_str());
		ck(conn, res);
		len = atoi(PQgetvalue(res, 0, 0));
		PQclear(res);
	}
	return len;
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

void readArray(char* start, vector<double>& array){
	char* pEnd = start;
	char* curPtr;
	while (*pEnd != '}'){
		curPtr = pEnd+1;
		array.push_back(strtod(curPtr, &pEnd));
	}
}

void SingleRow::getArray(int columnIndex, vector<double>& result){
	readArray(PQgetvalue(res, 0, columnIndex), result);
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

BulkCopy::~BulkCopy(){
	ck(pg->conn, PQputCopyEnd(pg->conn.get(), NULL) != 1);
	PGresult* res;
    while ((res = PQgetResult(pg->conn.get())) != NULL) {
		ck(pg->conn, res);
        PQclear(res);
    }
}

BulkCopy::BulkCopy(const string& tableName, const char& delimiter): tableName(tableName), delimiter(delimiter){
	pg = make_unique<PgManager>();
	string sql = fmt::format("COPY \"{}\" FROM STDIN with(delimiter '{}');", tableName, delimiter);
	auto res = PQexec(pg->conn.get(), sql.c_str());
	ck(pg->conn, res);
	PQclear(res);
	data = "";
	maxSize = 0;
}

void BulkCopy::send(){
	data.back() = '\n';
	maxSize = max(maxSize, data.size());
	ck(pg->conn, PQputCopyData(pg->conn.get(), data.c_str(), data.size()) != 1);
	data.clear();
	data.reserve(maxSize);
}




