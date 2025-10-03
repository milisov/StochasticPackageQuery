#include "formulator.hpp"
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
                            FormulateOptions& options);

    void formProbConsActiveness(GRBModel &model,
                            std::shared_ptr<Constraint> cons,
                            GRBVar *xx,
                            FormulateOptions& options);

};