#ifndef OBJ_HPP
#define OBJ_HPP

#include "util/udeclare.hpp"

#include <vector>
#include <string>
#include <memory>

using std::string;
using std::vector;
using std::shared_ptr;

class Objective{
public:
    ObjectiveSense objSense;
public:
    Objective(const string& objSense);
	virtual string toStr(const vector<double>& info={}) const = 0;
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
    string toStr(const vector<double>& info={}) const override;
};

class SumObjective: public Objective, public AttrObjective{
public:
    SumObjective(const string& objSense, const string& obj);
    string toStr(const vector<double>& info={}) const override;

};

class ExpectedSumObjective: public Objective, public AttrObjective{
public:
    ExpectedSumObjective(const string& objSense, const string& obj);
    string toStr(const vector<double>& info={}) const override;
};

class ProbObjective: public Objective, public AttrObjective{
public:
    Bound t;
    Inequality tsign;
    ProbObjective(const string& objSense, const string& obj, const Bound& t, const string& tsign);  
};

bool isDeterministic(const shared_ptr<Objective>& obj, shared_ptr<AttrObjective>& attrObj);
bool isDeterministic(const shared_ptr<Objective>& obj);

shared_ptr<CountObjective> getCount(const shared_ptr<Objective>& obj);

#endif