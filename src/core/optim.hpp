#ifndef OPTIM_HPP
#define OPTIM_HPP

#include "core/mapop.hpp"
#include "util/udeclare.hpp"

class Optim{
protected:
    SolType sol;
public:
    Optim(const SolType& sol);
    virtual SolType towards(const SolType& nextSol) = 0;
};

class RMSprop: public Optim{
private:
    SolType v;
    double alpha;
public:
    RMSprop(const SolType& sol={}, const double& alpha=0.99);
    SolType towards(const SolType& nextSol) override;
};

#endif