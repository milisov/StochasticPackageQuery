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
public:
    double lr, alpha;
    RMSprop(const SolIndType& sol={}, const double& lr=1.0, const double& alpha=0.95);
    SolIndType towards(const SolIndType& nextSol) override;
};

class Direct: public Optim{
public:
    Direct(const SolIndType& sol={});
    SolIndType towards(const SolIndType& nextSol) override;
};

#endif