#ifndef CHECKER_HPP
#define CHECKER_HPP

#include <memory>

#include "core/stat.hpp"
#include "util/udeclare.hpp"
#include "spq/spq.hpp"

using std::shared_ptr;
using std::unique_ptr;

class Checker{
protected:
    unique_ptr<Stat> stat;
public:
    Checker();
    virtual double getObjective(const SolType& sol) const = 0;
    virtual bool feasible(const SolType& sol, double &distance) const = 0;
    virtual void display(const SolType& sol) const = 0;
};

class SPQChecker: public Checker{
protected:
    shared_ptr<StochasticPackageQuery> spq;
    string validateTableName;
private:
    double getConIndicator(const SolType& sol, shared_ptr<Constraint> con) const;
public:
    SPQChecker(shared_ptr<StochasticPackageQuery> spq);
    double getObjective(const SolType& sol) const override;
    bool feasible(const SolType& sol, double &distance) const override;
    void display(const SolType& sol) const override;
};

#endif