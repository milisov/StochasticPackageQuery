#ifndef OBJ_HPP
#define OBJ_HPP

#include "util/udeclare.hpp"

#include <string>

using std::string;

class Objective{
public:
    ObjectiveSense objSense;
public:
    Objective(const string& objSense);
	virtual operator string() const = 0;
};

class AttrObjective{
public:
    string obj;
    Column objType;
    AttrObjective(const string& obj);
};

class CountObjective: public Objective{
public:
    CountObjective(const string& objSense);
    operator string() const override;
};

class SumObjective: public Objective, public AttrObjective{
public:
    SumObjective(const string& objSense, const string& obj);
    operator string() const override;
};

class ExpectedSumObjective: public Objective, public AttrObjective{
public:
    ExpectedSumObjective(const string& objSense, const string& obj);
    operator string() const override;
};

#endif