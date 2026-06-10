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

    // needed for deleting variance control constraints
    std::vector<GRBVar> varianceY;          // auxiliary variables y_i = sigma_i * x_i
    GRBVar varianceT;                        // norm variable t
    std::vector<GRBConstr> varianceYCons;   // constraints y_i = sigma_i * x_i
    GRBGenConstr varianceNormCon;           // norm constraint t = ||y||_2
    GRBConstr varianceControlCon;           // main constraint: x^T*mu - k*t >= v
    bool hasVarianceConstr = false;

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

    void formVarianceControlConstr(GRBModel &model,
                                   std::shared_ptr<Constraint> cons,
                                   GRBVar *xx,
                                   FormulateOptions &formOptions);

    void removeVarianceControlConstr(GRBModel &model);

    // MAD-based variance control
    std::vector<GRBVar> madD;              // auxiliary variables d_k for |deviation| in each scenario
    std::vector<GRBConstr> madDConsPos;    // d_k >= deviation_k
    std::vector<GRBConstr> madDConsNeg;    // d_k >= -deviation_k
    GRBConstr madControlCon;               // main constraint: x^T*mu - k*(1/n)*sum(d_k) >= v
    bool hasMADConstr = false;

    void formMADVarianceControlConstr(GRBModel &model,
                                      std::shared_ptr<Constraint> cons,
                                      GRBVar *xx,
                                      FormulateOptions &formOptions);

    void removeMADVarianceControlConstr(GRBModel &model);
};