#include <fmt/core.h>

#include "cons.hpp"

CountConstraint::CountConstraint(const Bound& lb, const Bound& ub): lb(lb), ub(ub){
}

CountConstraint::operator string() const {
    return fmt::format("COUNT(*) {}", to_string({lb, ub}));
}

SumConstraint::SumConstraint(const string& attr, const Bound& lb, const Bound& ub): attr(attr), lb(lb), ub(ub){
}

SumConstraint::operator string() const {
    return fmt::format("SUM({}) {}", attr, to_string({lb, ub}));
}

ExpectedSumConstraint::ExpectedSumConstraint(const string& attr, const Bound& lb, const Bound& ub): attr(attr), lb(lb), ub(ub){
}

ExpectedSumConstraint::operator string() const {
    return fmt::format("EXPECTED SUM({}) {}", attr, to_string({lb, ub}));
}

VarConstraint::VarConstraint(const string& attr, const Bound& v, const Bound& p, const string& vsign, const string& psign): attr(attr), v(v), p(p){
    this->vsign = toIneq[vsign];
    this->psign = toIneq[psign];
}

VarConstraint::operator string() const{
    return fmt::format("SUM({}) {} {} WITH PROBABILITY {} {}", attr, to_string(vsign), to_string(v), to_string(psign), to_string(p));
}