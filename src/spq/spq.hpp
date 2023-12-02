#ifndef SPQ_HPP
#define SPQ_HPP

#include <vector>
#include <string>
#include <memory>
#include <ostream>
#include <map>

#include "spq/cons.hpp"
#include "spq/obj.hpp"

using std::vector;
using std::string;
using std::shared_ptr;
using std::unique_ptr;
using std::ostream;
using std::map;

/**
 * @brief Represents SPQ's information connecting database and SPQ solver.
 * @details This is a result after parsing of an SPQ in .spaql file.
 * It can contains abstract information about the bounds of the constraints.
 * These bounds can be realized after calling StochasticQueryGenerator on this class
 * with appropriate hardness and expected size.
 * @param repeat if -1 then No upper bound
 * @param obj if nullptr then no objective, i.e., feasibility test
 */
class StochasticPackageQuery{
private:
	bool isValid;
	map<string, unique_ptr<double>> varTable;
private:
	bool addVariable(Bound bound);
public:
	static const int NO_REPEAT = -1;
	string tableName;
	vector<string> attrList;
	int repeat;
	vector<shared_ptr<Constraint>> cons;
	shared_ptr<Objective> obj;
public:
	StochasticPackageQuery();
	void setTableName(const string& tableName);
	void setAttrList(const vector<string>& attrList);
	void setRepeat(const int& repeat);
	void addConstraint(shared_ptr<Constraint> con);
	void setObjective(shared_ptr<Objective> obj);
	void setVariable(string var, double value);
	bool validate();
	bool executable();
	int countStochastic();
	double getValue(const Bound& bound) const;
	Bound getBound(const Bound& bound) const;
	operator string();
	string strAttrList() const;
	string substitute(shared_ptr<Constraint> con, const vector<double>& info={});
};

ostream& operator<<(ostream& os, const shared_ptr<StochasticPackageQuery> spq);
shared_ptr<StochasticPackageQuery> parseSpaqlFromFile(string filePath);

#endif