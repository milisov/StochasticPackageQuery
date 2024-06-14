#include <boost/algorithm/string/join.hpp>
#include <boost/variant.hpp>
#include <fmt/core.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>
#include <utility>

#include "util/udebug.hpp"
#include "util/uio.hpp"
#include "util/uconfig.hpp"
#include "spq.hpp"
#include "parser.hpp"

namespace fs = std::filesystem;

using boost::algorithm::join;
using std::cerr;
using std::pair;
using std::dynamic_pointer_cast;

bool StochasticPackageQuery::addVariable(const Bound& bound){
	if (bound.which() == 0){
		string variable = boost::get<string>(bound);
		if (varTable.count(variable)){
			cerr << fmt::format("Repeated bound variable '{}'\n", variable);
			return false;
		}
		varTable[variable] = nullptr;
	}
	return true;
}

string StochasticPackageQuery::substitute(const shared_ptr<Constraint>& con, const vector<double>& info){
	string res = "";
	shared_ptr<BoundConstraint> boundCon;
	if (isDeterministic(con, boundCon)){
		pair<Bound, Bound> bounds = {boundCon->lb, boundCon->ub};
		boundCon->lb = getBound(boundCon->lb);
		boundCon->ub = getBound(boundCon->ub);
		res = con->toStr(info);
		boundCon->lb = bounds.first;
		boundCon->ub = bounds.second;
		return res;
	}
	shared_ptr<ProbConstraint> probCon;
	if (isStochastic(con, probCon)){
		pair<Bound, Bound> bounds = {probCon->v, probCon->p};
		probCon->v = getBound(probCon->v);
		probCon->p = getBound(probCon->p);
		res = con->toStr(info);
		probCon->v = bounds.first;
		probCon->p = bounds.second;
		return res;
	}
	cerr << fmt::format("Constraint '{}' not supported for substitution\n", con->toStr());
	exit(1);
}

StochasticPackageQuery::StochasticPackageQuery(){
	repeat = StochasticPackageQuery::NO_REPEAT;
	isValid = false;
}

void StochasticPackageQuery::setTableName(const string& tableName){
	this->tableName = tableName;
}

void StochasticPackageQuery::setAttrList(const vector<string>& attrList){
	this->attrList = attrList;
}

void StochasticPackageQuery::setRepeat(const int& repeat){
	this->repeat = repeat;
}

void StochasticPackageQuery::addConstraint(const shared_ptr<Constraint>& con){
	cons.push_back(con);
}

void StochasticPackageQuery::setObjective(const shared_ptr<Objective>& obj){
	this->obj = obj;
}

void StochasticPackageQuery::setVariable(const string& var, const double& value){
	varTable[var] = std::move(std::make_unique<double>(value));
}

bool StochasticPackageQuery::validate(){
	if (isValid) return isValid;
	varTable.clear();
	auto pg = PgManager();
	if (!pg.existTable(tableName)){
		cerr << fmt::format("Table '{}' does not exist\n", tableName);
		return false;
	}
	auto columns = pg.getColumns(tableName);
	if (attrList.size()){
		for (const auto& attr : attrList){
			if (!columns.count(attr)){
				cerr << fmt::format("Column '{}' in table '{}' is not found\n", attr, tableName);
				return false;
			}
			if (columns[attr] == Column::unsupported){
				cerr << fmt::format("Column '{}' in table '{}' is not supported\n", attr, tableName);
				return false;
			}
		}
	} else {
		for (const auto& column : columns){
			if (column.second != Column::unsupported) attrList.push_back(column.first);
		}
	}
	if (repeat < NO_REPEAT){
		cerr << fmt::format("Invalid repeat value '{}'\n", repeat);
		return false;
	}
	for (size_t i = 0; i < cons.size(); ++i){
		const auto& con = cons[i];
		if (!con){
			cerr << "Encountered null constraint\n";
			return false;
		}
		shared_ptr<BoundConstraint> boundCon = dynamic_pointer_cast<BoundConstraint>(con);
		if (boundCon){
			if (!addVariable(boundCon->lb) || !addVariable(boundCon->ub)) return false;
			if (boundCon->lb.which() != boundCon->ub.which()){
				if ((boundCon->lb.which() == 1 && boost::get<double>(boundCon->lb) != NEG_INF)
				| (boundCon->ub.which() == 1 && boost::get<double>(boundCon->ub) != POS_INF)){
					cerr << fmt::format("Constraint '{}' has both numeric and variable bounds\n", con->toStr());
					return false;	
				}
			}
		}
		shared_ptr<ProbConstraint> probCon = dynamic_pointer_cast<ProbConstraint>(con);
		if (probCon){
			if (!addVariable(probCon->v)) return false;
			if (!hasValue(probCon->p)){
				cerr << fmt::format("Variable '{}' is not supported in constraint '{}'\n", boost::get<string>(probCon->p), con->toStr());
				return false;
			}
			double p = getValue(probCon->p);
			if (p <= 0 || p >= 1){
				cerr << fmt::format("Invalid probability '{}' in constraint '{}'\n", p, con->toStr());
				return false;
			}
			if (probCon->p.which() == 1){
				auto pVar = fmt::format("_p{}", i);
				setVariable(pVar, p);
				probCon->p = Bound(pVar);
			}
		}
		shared_ptr<AttrConstraint> attrCon = dynamic_pointer_cast<AttrConstraint>(con);
		if (attrCon){
			if (!columns.count(attrCon->attr)){
				cerr << fmt::format("No column '{}' exists in the contraint '{}'\n", attrCon->attr, con->toStr());
				return false;
			}
			if (columns[attrCon->attr] != attrCon->attrType){
				cerr << fmt::format("Column '{}' is not supported in the contraint '{}'\n", attrCon->attr, con->toStr());
				return false;
			}
		}
	}
	if (obj){
		shared_ptr<AttrObjective> attrObj = dynamic_pointer_cast<AttrObjective>(obj);
		if (attrObj){
			if (!columns.count(attrObj->obj)){
				cerr << fmt::format("No column '{}' exists in the objective '{}'\n", attrObj->obj, obj->toStr());
				return false;
			}
			if (columns[attrObj->obj] != attrObj->objType){
				cerr << fmt::format("Column '{}' is not supported in the objective '{}'\n", attrObj->obj, obj->toStr());
				return false;
			}
		}
	} else{
		cerr << fmt::format("No objective\n");
		return false;
	}
	isValid = true;
	return true;
}

int StochasticPackageQuery::getCvarCount() const{
	int res = 0;
	for (const auto& con : cons){
		if (getCvar(con)) res ++;
	}
	return res;
}

bool StochasticPackageQuery::executable() const{
	for (const auto& p : varTable){
		if (!(p.second)) return false;
	}
	return true;
}

StochasticPackageQuery::operator string(){
	string res = fmt::format("SELECT PACKAGE({}) FROM {}", strAttrList(), tableName);
	if (cons.size()) res += " SUCH THAT\n";
	vector<string> strCons;
	for (const auto& con : cons){
		if (con) strCons.push_back("\t" + substitute(con));
	}
	res += join(strCons, " AND\n") + '\n';
	if (obj) res += obj->toStr() + '\n';
	return res;
}

bool StochasticPackageQuery::hasValue(const Bound& bound) const{
	if (bound.which() == 1) return true;
	string var = boost::get<string>(bound);
	if (varTable.count(var) && varTable.at(var)) return true;
	return false;
}

bool StochasticPackageQuery::isStochasticObjective() const{
	return !isDeterministic(obj);
}

double StochasticPackageQuery::getValue(const Bound& bound) const{
	if (bound.which() == 1) return boost::get<double>(bound);
	string var = boost::get<string>(bound);
	if (varTable.count(var) && varTable.at(var)) return *(varTable.at(var));
	cerr << fmt::format("Variable {} is not set\n", var);
	exit(1);
}

Bound StochasticPackageQuery::getBound(const Bound& bound) const{
	if (bound.which() == 1) return bound;
	string var = boost::get<string>(bound);
	if (varTable.count(var) && varTable.at(var)) return Bound(*(varTable.at(var)));
	return bound;
}

string StochasticPackageQuery::strAttrList() const{
	string res = join(attrList, ", ");
	if (!res.size()) res += "*";
	return res;
}

ostream& operator<<(ostream& os, const shared_ptr<StochasticPackageQuery>& spq){
	if (spq) os << static_cast<string>(*spq);
    return os;
}

shared_ptr<StochasticPackageQuery> parseSpaqlFromFile(const string& filePath){
	fs::path spaqlPath (filePath);
	if (spaqlPath.is_relative()) spaqlPath = getProjectDir() / spaqlPath;
	std::ifstream in(spaqlPath);
	if (!in.is_open()) {
		cerr << fmt::format("Could not open file '{}'\n", spaqlPath.string());
		exit(1);
	}
	std::stringstream buffer;
	buffer << in.rdbuf();
	string spaqlQuery = buffer.str();
	in.close();
	return parseSpaql(spaqlQuery);
}
