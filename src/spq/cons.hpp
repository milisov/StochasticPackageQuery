#ifndef CONS_HPP
#define CONS_HPP

#include <vector>
#include <string>

#include "util/udeclare.hpp"

using std::string;
using std::vector;

class Constraint{
public:
	virtual string toStr(const vector<double>& info={}) const = 0;
    virtual bool isViolate(const vector<double>& info) const = 0;
};

class BoundConstraint : public Constraint {
public:
    Bound lb, ub;
    BoundConstraint(const Bound& lb, const Bound& ub);
    bool isViolate(const vector<double>& info) const override;
};

class ProbConstraint : public Constraint {
public:
    Bound v, p;
    Inequality vsign, psign;
    ProbConstraint(const Bound& v, const Bound& p, const string& vsign, const string& psign);
};

class AttrConstraint{
public:
    string attr;
    Column attrType;
public:
    AttrConstraint(const string& attr);
};

class CountConstraint: public BoundConstraint{
public:
    CountConstraint(const Bound& lb, const Bound& ub);
    string toStr(const vector<double>& info={}) const override;
};

class SumConstraint: public AttrConstraint, public BoundConstraint{
public:
    SumConstraint(const string& attr, const Bound& lb, const Bound& ub);
    string toStr(const vector<double>& info={}) const override;
};

class ExpectedSumConstraint: public AttrConstraint, public BoundConstraint{
public:
    ExpectedSumConstraint(const string& attr, const Bound& lb, const Bound& ub);
    string toStr(const vector<double>& info={}) const override;
};

class VarConstraint: public AttrConstraint, public ProbConstraint{
public:
    VarConstraint(const string& attr, const Bound& v, const Bound& p, const string& vsign, const string& psign);
    string toStr(const vector<double>& info={}) const override;
    bool isViolate(const vector<double>& info) const override;
};

#endif