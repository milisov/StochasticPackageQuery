#include <boost/algorithm/string/join.hpp>
#include <fmt/core.h>
#include <algorithm>
#include <iostream>

#include "util/udebug.hpp"
#include "util/uio.hpp"
#include "spq.hpp"

using boost::algorithm::join;
using std::transform;
using std::back_inserter;
using std::cerr;
using std::endl;

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

bool StochasticPackageQuery::validate(){
	
	return true;
}

StochasticPackageQuery::operator string() const{
	string res = fmt::format("SELECT PACKAGE({}) FROM {}", strAttrList(), tableName);
	if (cons.size()) res += " SUCH THAT\n";
	vector<string> strCons;
	for (auto con : cons) if (con) strCons.push_back("\t" + static_cast<string>(*con));
	res += join(strCons, " AND\n") + '\n';
	if (obj) res += static_cast<string>(*obj) + '\n';
	return res;
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
