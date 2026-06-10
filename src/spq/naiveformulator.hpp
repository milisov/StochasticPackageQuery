#include "formulator.hpp"
#include <set>
#pragma once

class NaiveFormulator : public Formulator {
public:
    NaiveFormulator(std::shared_ptr<StochasticPackageQuery> spqPtr): Formulator(spqPtr)   // <-- this ensures Formulator's constructor is executed
    {
        //initStockExpectedProfit(DB_optim, "profit");
        // Additional derived-class initialization
    };

    // Build the model
    GRBModel formulate(shared_ptr<StochasticPackageQuery> spq, 
        FormulateOptions& formOptions) override;

    // Implementation of probabilistic constraints for Naive
    void formProbCons(GRBModel &model,
                            std::shared_ptr<Constraint> cons,
                            GRBVar *xx,
                            FormulateOptions& options, int conOrder);

    void formProbConsActiveness(GRBModel &model,
                            std::shared_ptr<Constraint> cons,
                            GRBVar *xx,
                            FormulateOptions& options, int conOrder);

    // MAD-based variance control with reduced scenarios
    void formMADVarianceControlConstr(GRBModel &model,
                                      std::shared_ptr<Constraint> cons,
                                      GRBVar *xx,
                                      FormulateOptions &options,
                                      int &conOrder);

    void randomScenarioSelection(int Z, int conOrder, std::set<int> &selectedScenarios, unsigned int seed = 42);
};