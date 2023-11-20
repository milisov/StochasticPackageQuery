#include <fmt/core.h>

#include "cons.hpp"

BoundConstraint::BoundConstraint(const Bound& lb, const Bound& ub) : lb(lb), ub(ub) {
}

ProbConstraint::ProbConstraint(const Bound& v, const Bound& p, const string& vsign, const string& psign) : v(v), p(p) {
    this->vsign = toIneq[vsign];
    this->psign = toIneq[psign];
}

AttrConstraint::AttrConstraint(const string& attr): attr(attr){
}

CountConstraint::CountConstraint(const Bound& lb, const Bound& ub): BoundConstraint(lb, ub){
}

CountConstraint::operator string() const {
    return fmt::format("COUNT(*) {}", to_string({lb, ub}));
}

SumConstraint::SumConstraint(const string& attr, const Bound& lb, const Bound& ub): AttrConstraint(attr), BoundConstraint(lb, ub){
    attrType = Column::numeric_type;
}

SumConstraint::operator string() const {
    return fmt::format("SUM({}) {}", attr, to_string({lb, ub}));
}

ExpectedSumConstraint::ExpectedSumConstraint(const string& attr, const Bound& lb, const Bound& ub): AttrConstraint(attr), BoundConstraint(lb, ub){
    attrType = Column::array_type;
}

ExpectedSumConstraint::operator string() const {
    return fmt::format("EXPECTED SUM({}) {}", attr, to_string({lb, ub}));
}

VarConstraint::VarConstraint(const string& attr, const Bound& v, const Bound& p, const string& vsign, const string& psign): AttrConstraint(attr), ProbConstraint(v, p, vsign, psign){
    attrType = Column::array_type;
}

VarConstraint::operator string() const{
    return fmt::format("SUM({}) {} {} WITH PROBABILITY {} {}", attr, to_string(vsign), to_string(v), to_string(psign), to_string(p));
}