#ifndef CONS_HPP
#define CONS_HPP

#include <string>

#include "util/udeclare.hpp"

using std::string;

class Constraint{
public:
	virtual operator string() const = 0;
};

class CountConstraint: public Constraint{
public:
    Bound lb, ub;
public:
    CountConstraint(const Bound& lb, const Bound& ub);
    operator string() const;
};

class SumConstraint: public Constraint{
public:
    string attr;
    Bound lb, ub;  
public:
    SumConstraint(const string& attr, const Bound& lb, const Bound& ub);
    operator string() const;    
};

class ExpectedSumConstraint: public Constraint{
public:
    string attr;
    Bound lb, ub;  
public:
    ExpectedSumConstraint(const string& attr, const Bound& lb, const Bound& ub);
    operator string() const;    
};

class VarConstraint: public Constraint{
public:
    string attr;
    Bound v, p;
    Ineq vsign, psign;
public:
    VarConstraint(const string& attr, const Bound& v, const Bound& p, const string& vsign, const string& psign);
    operator string() const; 
};

#endif