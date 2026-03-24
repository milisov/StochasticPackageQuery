#include "formulator.hpp"
#pragma once

class SSFormulator : public Formulator
{
public:
    GRBModel model;
    // needed for deleting probConstraints
    std::vector<GRBVar> yy;
    std::vector<GRBGenConstr> genCon;
    std::vector<GRBConstr> probCons;

    SSFormulator(std::shared_ptr<StochasticPackageQuery> spqPtr) : Formulator(spqPtr), model(env)// <-- this ensures Formulator's constructor is executed
    {
        // model.set(GRB_IntParam_OutputFlag, 0);
        populateShuffler(shuffler);
        model.set(GRB_DoubleParam_TimeLimit, 600);
    };

    GRBModel& getModel() override {
        cout<<"I am returning the model"<<endl;
        return model;
    }

    void updateBound(GRBModel &model, double ub);

    GRBModel formulate(shared_ptr<StochasticPackageQuery> spq,
                       FormulateOptions &formOptions);

    void formProbCons(GRBModel &model, shared_ptr<Constraint> cons, GRBVar *xx, std::vector<std::vector<std::vector<double>>> &summaries,
                      int &probConOrder, FormulateOptions &formOptions);

    void formulateSAA(GRBModel &model, std::vector<std::vector<std::vector<double>>> &summaries, GRBVar *xx, FormulateOptions &formOptions);
    GRBModel formulateDeterministicLP(std::shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions);
    GRBModel formulateCVaRLP(std::shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions);
    void formProbConsStage3(GRBModel &model,
                            std::shared_ptr<Constraint> cons,
                            GRBVar *xx,
                            std::vector<std::vector<std::vector<double>>> &summaries,
                            int &probConOrder, FormulateOptions &formOptions);
};