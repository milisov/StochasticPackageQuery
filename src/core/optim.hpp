#ifndef OPTIM_HPP
#define OPTIM_HPP

#include "core/mapop.hpp"
#include "util/udeclare.hpp"

class Optim{
protected:
    SolIndType sol;
public:
    Optim(const SolIndType& sol);
    virtual SolIndType towards(const SolIndType& nextSol) = 0;
};

class RMSprop: public Optim{
private:
    SolIndType v;
    double alpha;
public:
    RMSprop(const SolIndType& sol={}, const double& alpha=0.99);
    SolIndType towards(const SolIndType& nextSol) override;
};

#endif