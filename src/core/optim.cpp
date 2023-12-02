#include "util/unumeric.hpp"
#include "util/udebug.hpp"

#include "optim.hpp"

Optim::Optim(const SolIndType& sol): sol(sol){
}

RMSprop::RMSprop(const SolIndType& sol, const double& alpha): Optim(sol), alpha(alpha){
}

SolIndType RMSprop::towards(const SolIndType& nextSol){
    auto g = add(1.0, nextSol, -1.0, sol);
    v = add(alpha, v, 1-alpha, mul(g, g));
    auto step = div(g, add(sqrt(v), MACHINE_EPS));
    sol = add(1.0, sol, 1.0, step);
    return step;
}