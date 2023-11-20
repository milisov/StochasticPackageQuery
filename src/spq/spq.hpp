#ifndef SPQ_HPP
#define SPQ_HPP

#include <vector>
#include <string>
#include <memory>
#include <ostream>

#include "cons.hpp"
#include "obj.hpp"

using std::vector;
using std::string;
using std::shared_ptr;
using std::ostream;

/**
 * @brief Represents SPQ's information connecting database and SPQ solver.
 * @details This is a result after parsing of an SPQ in .spaql file.
 * It can contains abstract information about the bounds of the constraints.
 * These bounds can be realized after calling StochasticQueryGenerator on this class
 * with appropriate hardness and expected size.
 */
class StochasticPackageQuery{
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
	bool validate();
	operator string() const;
	string strAttrList() const;
};

ostream& operator<<(ostream& os, const shared_ptr<StochasticPackageQuery> spq);

#endif