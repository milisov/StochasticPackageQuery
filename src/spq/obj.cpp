#include <fmt/core.h>
#include <boost/algorithm/string.hpp>

#include "obj.hpp"

using boost::algorithm::to_upper_copy;

Objective::Objective(const string& objSense){
    this->objSense = toObjSense[to_upper_copy(objSense)];
}

AttrObjective::AttrObjective(const string& obj): obj(obj){
}

CountObjective::CountObjective(const string& objSense): Objective(objSense){
}

CountObjective::operator string() const{
    return fmt::format("{} COUNT(*)", to_string(objSense));
}

SumObjective::SumObjective(const string& objSense, const string& obj): Objective(objSense), AttrObjective(obj){
    objType = Column::numeric_type;
}

SumObjective::operator string() const{
    return fmt::format("{} SUM({})", to_string(objSense), obj);
}

ExpectedSumObjective::ExpectedSumObjective(const string& objSense, const string& obj): Objective(objSense), AttrObjective(obj){
    objType = Column::array_type;
}

ExpectedSumObjective::operator string() const{
    return fmt::format("{} EXPECTED SUM({})", to_string(objSense), obj);
}