#include <fmt/core.h>
#include <boost/algorithm/string.hpp>

#include "obj.hpp"

using boost::algorithm::to_upper_copy;

CountObjective::CountObjective(const string& objSense){
    this->objSense = toObjSense[to_upper_copy(objSense)];
}

CountObjective::operator string() const{
    return fmt::format("{} COUNT(*)", to_string(objSense));
}

SumObjective::SumObjective(const string& objSense, const string& obj): obj(obj){
    this->objSense = toObjSense[to_upper_copy(objSense)];
}

SumObjective::operator string() const{
    return fmt::format("{} SUM({})", to_string(objSense), obj);
}

ExpectedSumObjective::ExpectedSumObjective(const string& objSense, const string& obj): obj(obj){
    this->objSense = toObjSense[to_upper_copy(objSense)];
}

ExpectedSumObjective::operator string() const{
    return fmt::format("{} EXPECTED SUM({})", to_string(objSense), obj);
}