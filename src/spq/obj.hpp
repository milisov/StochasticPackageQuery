#ifndef OBJ_HPP
#define OBJ_HPP

#include "util/udeclare.hpp"

#include <string>

using std::string;

class Objective{
public:
	virtual operator string() const = 0;
};

class CountObjective: public Objective{
public:
    ObjSense objSense;
public:
    CountObjective(const string& objSense);
    operator string() const;
};

class SumObjective: public Objective{
public:
    ObjSense objSense;
    string obj;
public:
    SumObjective(const string& objSense, const string& obj);
    operator string() const;
};

class ExpectedSumObjective: public Objective{
public:
    ObjSense objSense;
    string obj;
public:
    ExpectedSumObjective(const string& objSense, const string& obj);
    operator string() const;
};

#endif