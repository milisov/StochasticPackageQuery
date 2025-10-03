#include "formulator.hpp"
#pragma once

class RSFormulator : public Formulator
{
public:
    GRBModel modelRS;
    GRBModel modelBestObj;
    
    RSFormulator(std::shared_ptr<StochasticPackageQuery> spqPtr) : modelRS(env), modelBestObj(env), Formulator(spqPtr) // <-- this ensures Formulator's constructor is executed
    {
        modelRS.set(GRB_IntParam_OutputFlag, 1); 
        modelBestObj.set(GRB_IntParam_OutputFlag, 1);
        populateShuffler(shuffler);
        // Additional derived-class initialization
    };

    // Build the model
    GRBModel formulate(shared_ptr<StochasticPackageQuery> spq,
                       FormulateOptions &formOptions) override;

    void formProbConsSAA(GRBModel &model,std::shared_ptr<Constraint> cons,
                      GRBVar *xx,int &probConOrder, FormulateOptions &formOptions);     

    void formulateBestObjProblem(GRBModel &model, std::shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions);
    void formulateDeterministic(GRBModel &model, std::shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions);
    void formulateRS(GRBModel &model, std::shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions);

    void removeLastProbConstraints(GRBModel &model);
    void removeMinMaxObjective(GRBModel &model);

    // Implementation of probabilistic constraints for Robust Satisficing
    void formProbCons(GRBModel &model,
                      std::shared_ptr<Constraint> cons,
                      GRBVar *xx,
                      std::vector<std::vector<std::vector<double>>> &summaries,
                      int &probConOrder, FormulateOptions &formOptions);

    void formObjCons(GRBModel &model,
                     std::shared_ptr<Objective> obj,
                     GRBVar *xx,
                     FormulateOptions &formOptions);

    void formMinMaxObjective(GRBModel &model, GRBVar *xx, FormulateOptions &formOptions);

};