#include <fmt/core.h>
#include <boost/variant.hpp>

#include "util/unumeric.hpp"
#include "util/udebug.hpp"
#include "cons.hpp"

using std::dynamic_pointer_cast;
using std::to_string;

BoundConstraint::BoundConstraint(const Bound& lb, const Bound& ub) : lb(lb), ub(ub) {
}

bool BoundConstraint::isViolate(const vector<double>& info) const{
    deb(info);
    if (!info.size()) return false;
    if (isLess(info.front(), info[1])) return true;
    if (isGreater(info.front(), info[2])) return true;
    return false;
}

ProbConstraint::ProbConstraint(const Bound& v, const Bound& p, const string& vsign, const string& psign) : v(v), p(p) {
    this->vsign = toInequality[vsign];
    this->psign = toInequality[psign];
}

AttrConstraint::AttrConstraint(const string& attr): attr(attr){
}

CountConstraint::CountConstraint(const Bound& lb, const Bound& ub): BoundConstraint(lb, ub){
}

string CountConstraint::toStr(const vector<double>& info) const{
    if (!info.size()) return fmt::format("COUNT(*) {}", str({lb, ub}));
    auto color = isViolate(info) ? RED : GREEN;
    return fmt::format("COUNT(*)={}{}{} {}", color, info.front(), RESET, str({lb, ub}));
}

SumConstraint::SumConstraint(const string& attr, const Bound& lb, const Bound& ub): AttrConstraint(attr), BoundConstraint(lb, ub){
    attrType = Column::numeric_type;
}

string SumConstraint::toStr(const vector<double>& info) const{
    if (!info.size()) return fmt::format("SUM({}) {}", attr, str({lb, ub}));
    auto color = isViolate(info) ? RED : GREEN;
    return fmt::format("SUM({})={}{}{} {}", attr, color, info.front(), RESET, str({lb, ub}));
}

ExpectedSumConstraint::ExpectedSumConstraint(const string& attr, const Bound& lb, const Bound& ub): AttrConstraint(attr), BoundConstraint(lb, ub){
    attrType = Column::array_type;
}

string ExpectedSumConstraint::toStr(const vector<double>& info) const{
    if (!info.size()) return fmt::format("EXPECTED SUM({}) {}", attr, str({lb, ub}));
    auto color = isViolate(info) ? RED : GREEN;
    return fmt::format("EXPECTED SUM({})={}{}{} {}", attr, color, info.front(), RESET, str({lb, ub}));
}

VarConstraint::VarConstraint(const string& attr, const Bound& v, const Bound& p, const string& vsign, const string& psign): AttrConstraint(attr), ProbConstraint(v, p, vsign, psign){
    attrType = Column::array_type;
}

bool VarConstraint::isViolate(const vector<double>& info) const{
    deb(info);
    if (!info.size()) return false;
    if (psign == Inequality::gteq && isLess(info.front(), info[1])) return true;
    if (psign == Inequality::lteq && isGreater(info.front(), info[1])) return true;
    return false;
}

string VarConstraint::toStr(const vector<double>& info) const{
    if (!info.size()) return fmt::format("SUM({}) {} {} WITH PROBABILITY {} {}", attr, str(vsign), str(v), str(psign), str(p));
    auto color = isViolate(info) ? RED : GREEN;
    return fmt::format("SUM({}) {} {} WITH PROBABILITY={}{}{} {} {}", attr, str(vsign), str(v), color, info.front(), RESET, str(psign), str(p));
}

CvarConstraint::CvarConstraint(const string& attr, const Bound& v, const Bound& p, const string& vsign, const string& psign): AttrConstraint(attr), ProbConstraint(v, p, vsign, psign){
    attrType = Column::array_type;
}

bool CvarConstraint::isViolate(const vector<double>& info) const{
    if (!info.size()) return false;
    if (vsign == Inequality::gteq && isLess(info.front(),info[2])) return true;
    if (vsign == Inequality::lteq && isGreater(info.front(), info[2])) return true;
    return false;
}

string CvarConstraint::toStr(const vector<double>& info) const{
    string strPsign = psign == Inequality::lteq ? "LOWEST" : "HIGHEST";
    string strP;
    if (p.which() == 0) strP = "100" + boost::get<string>(p);
    else if (p.which() == 1) strP = to_string(boost::get<double>(p)*100);
    if (!info.size()) return fmt::format("EXPECTED SUM({}) {} {} IN {} {}% OF CASES", attr, str(vsign), str(v), strPsign, strP);
    auto color = isViolate(info) ? RED : GREEN;
    return fmt::format("EXPECTED SUM({})={}{}{} {} {} IN {} {}% OF CASES", attr, color, info.front(), RESET, str(vsign), str(v), strPsign, strP);
}

bool isStochastic(const shared_ptr<Constraint>& con, shared_ptr<ProbConstraint>& probCon, shared_ptr<AttrConstraint>& attrCon){
    probCon = dynamic_pointer_cast<ProbConstraint>(con);
    attrCon = dynamic_pointer_cast<AttrConstraint>(con);
    return probCon != nullptr;
}

bool isDeterministic(const shared_ptr<Constraint>& con, shared_ptr<BoundConstraint>& boundCon, shared_ptr<AttrConstraint>& attrCon){
    boundCon = dynamic_pointer_cast<BoundConstraint>(con);
    attrCon = dynamic_pointer_cast<AttrConstraint>(con);
    return boundCon != nullptr;
}

bool isStochastic(const shared_ptr<Constraint>& con, shared_ptr<ProbConstraint>& probCon){
    probCon = dynamic_pointer_cast<ProbConstraint>(con);
    return probCon != nullptr;
}

bool isDeterministic(const shared_ptr<Constraint>& con, shared_ptr<BoundConstraint>& boundCon){
    boundCon = dynamic_pointer_cast<BoundConstraint>(con);
    return boundCon != nullptr;
}

bool isStochastic(const shared_ptr<Constraint>& con){
    auto probCon = dynamic_pointer_cast<ProbConstraint>(con);
    return probCon != nullptr;
}

bool isDeterministic(const shared_ptr<Constraint>& con){
    auto boundCon = dynamic_pointer_cast<BoundConstraint>(con);
    return boundCon != nullptr;
}

shared_ptr<VarConstraint> getVar(const shared_ptr<Constraint>& con){
    return dynamic_pointer_cast<VarConstraint>(con);
}

shared_ptr<CvarConstraint> getCvar(const shared_ptr<Constraint>& con){
    return dynamic_pointer_cast<CvarConstraint>(con);
}

shared_ptr<CountConstraint> getCount(const shared_ptr<Constraint>& con){
    return dynamic_pointer_cast<CountConstraint>(con);
}