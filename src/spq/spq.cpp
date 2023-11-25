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

bool StochasticPackageQuery::addVariable(Bound bound){
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

string StochasticPackageQuery::substitute(shared_ptr<Constraint> con){
	string res = "";
	shared_ptr<BoundConstraint> boundCon = dynamic_pointer_cast<BoundConstraint>(con);
	if (boundCon){
		pair<Bound, Bound> bounds = {boundCon->lb, boundCon->ub};
		boundCon->lb = getBound(boundCon->lb);
		boundCon->ub = getBound(boundCon->ub);
		res = static_cast<string>(*con);
		boundCon->lb = bounds.first;
		boundCon->ub = bounds.second;
		return res;
	}
	shared_ptr<ProbConstraint> probCon = dynamic_pointer_cast<ProbConstraint>(con);
	if (probCon){
		Bound bound = probCon->v;
		probCon->v = getBound(probCon->v);
		res = static_cast<string>(*con);
		probCon->v = bound;
		return res;
	}
	cerr << fmt::format("Constraint '{}' not supported for substitution\n", static_cast<string>(*con));
	exit(1);
}

StochasticPackageQuery::StochasticPackageQuery(){
	repeat = StochasticPackageQuery::NO_REPEAT;
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

void StochasticPackageQuery::addConstraint(shared_ptr<Constraint> con){
	cons.push_back(con);
}

void StochasticPackageQuery::setObjective(shared_ptr<Objective> obj){
	this->obj = obj;
}

void StochasticPackageQuery::setVariable(string var, double value){
	if (!varTable.count(var)){
		cerr << fmt::format("Variable '{}' is not in the SPQ '{}' (If not, it means you have not validated the SPQ)\n", var, operator string());
		exit(1);
	}
	varTable[var] = std::move(std::make_unique<double>(value));
}

bool StochasticPackageQuery::validate(){
	varTable.clear();
	auto pg = PgManager();
	if (!pg.existTable(tableName)){
		cerr << fmt::format("Table '{}' does not exist\n", tableName);
		return false;
	}
	auto columns = pg.getColumns(tableName);
	if (attrList.size()){
		for (auto attr : attrList){
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
		for (auto column : columns){
			if (column.second != Column::unsupported) attrList.push_back(column.first);
		}
	}
	if (repeat < NO_REPEAT){
		cerr << fmt::format("Invalid repeat value '{}'\n", repeat);
		return false;
	}
	for (auto con : cons){
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
					cerr << fmt::format("Constraint '{}' has both numeric and variable bounds\n", static_cast<string>(*con));
					return false;	
				}
			}
		}
		shared_ptr<ProbConstraint> probCon = dynamic_pointer_cast<ProbConstraint>(con);
		if (probCon){
			if (!addVariable(probCon->v)) return false;
			if (probCon->p.which() == 0){
				cerr << fmt::format("Variable '{}' is not supported in constraint '{}'\n", boost::get<string>(probCon->p), static_cast<string>(*con));
				return false;
			}
			double p = boost::get<double>(probCon->p);
			if (p <= 0 || p >= 1){
				cerr << fmt::format("Invalid probability '{}' in constraint '{}'\n", p, static_cast<string>(*con));
				return false;
			}
		}
		shared_ptr<AttrConstraint> attrCon = dynamic_pointer_cast<AttrConstraint>(con);
		if (attrCon){
			if (!columns.count(attrCon->attr)){
				cerr << fmt::format("No column '{}' exists in the contraint '{}'\n", attrCon->attr, static_cast<string>(*con));
				return false;
			}
			if (columns[attrCon->attr] != attrCon->attrType){
				cerr << fmt::format("Column '{}' is not supported in the contraint '{}'\n", attrCon->attr, static_cast<string>(*con));
				return false;
			}
		}
	}
	if (obj){
		shared_ptr<AttrObjective> attrObj = dynamic_pointer_cast<AttrObjective>(obj);
		if (attrObj){
			if (!columns.count(attrObj->obj)){
				cerr << fmt::format("No column '{}' exists in the objective '{}'\n", attrObj->obj, static_cast<string>(*obj));
				return false;
			}
			if (columns[attrObj->obj] != attrObj->objType){
				cerr << fmt::format("Column '{}' is not supported in the objective '{}'\n", attrObj->obj, static_cast<string>(*obj));
				return false;
			}
		}
	}
	return true;
}

bool StochasticPackageQuery::executable(){
	for (const auto& p : varTable){
		if (!(p.second)) return false;
	}
	return true;
}

int StochasticPackageQuery::countStochastic(){
	int res = 0;
	for (auto con : cons){
		shared_ptr<ProbConstraint> probCon = dynamic_pointer_cast<ProbConstraint>(con);
		if (probCon) res ++;
	}
	return res;
}

StochasticPackageQuery::operator string(){
	string res = fmt::format("SELECT PACKAGE({}) FROM {}", strAttrList(), tableName);
	if (cons.size()) res += " SUCH THAT\n";
	vector<string> strCons;
	for (auto con : cons){
		if (con) strCons.push_back("\t" + substitute(con));
	}
	res += join(strCons, " AND\n") + '\n';
	if (obj) res += static_cast<string>(*obj) + '\n';
	return res;
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

ostream& operator<<(ostream& os, const shared_ptr<StochasticPackageQuery> spq){
	if (spq) os << static_cast<string>(*spq);
    return os;
}

shared_ptr<StochasticPackageQuery> parseSpaqlFromFile(string filePath){
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
