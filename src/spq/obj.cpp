#include <fmt/core.h>
#include <boost/algorithm/string.hpp>

#include "obj.hpp"

using boost::algorithm::to_upper_copy;

Objective::Objective(const string& objSense){
    this->objSense = toObjectiveSense[to_upper_copy(objSense)];
}

AttrObjective::AttrObjective(const string& obj): obj(obj){
}

CountObjective::CountObjective(const string& objSense): Objective(objSense){
}

string CountObjective::toStr(const vector<double>& info) const{
    if (!info.size()) return fmt::format("{} COUNT(*)", str(objSense));
    return fmt::format("{} COUNT(*)={}", str(objSense), info.front());
}

SumObjective::SumObjective(const string& objSense, const string& obj): Objective(objSense), AttrObjective(obj){
    objType = Column::numeric_type;
}

string SumObjective::toStr(const vector<double>& info) const{
    if (!info.size()) return fmt::format("{} SUM({})", str(objSense), obj);
    return fmt::format("{} SUM({})={}", str(objSense), obj, info.front());
}

ExpectedSumObjective::ExpectedSumObjective(const string& objSense, const string& obj): Objective(objSense), AttrObjective(obj){
    objType = Column::array_type;
}

string ExpectedSumObjective::toStr(const vector<double>& info) const{
    if (!info.size()) return fmt::format("{} EXPECTED SUM({})", str(objSense), obj);
    return fmt::format("{} EXPECTED SUM({})={}", str(objSense), obj, info.front());
}