#include "sdrformulator.hpp"


GRBModel SDRFormulator::formulate(std::shared_ptr<StochasticPackageQuery> spq,
                                  FormulateOptions& formOptions)
{

    GRBModel model(env);
    std::unique_ptr<GRBVar[]> xx;
    int n = formOptions.reduced ? formOptions.reducedIds.size() : NTuples;
    deb(n);
    xx = std::make_unique<GRBVar[]>(n);
    DecisionVarOptions decVarOptions = formOptions.decisionVarOptions;
    for (int i = 0; i < n; i++)
    {
        decVarOptions.name = "xx[" + to_string(i) + "]";
        xx[i] = addDecisionVar(model, decVarOptions);
    }

    int numCons = spq->cons.size();
    int probConOrder = 0;
    for (int i = 0; i < numCons; i++)
    {
        formLCVaR(model, spq->cons[i], xx.get(), formOptions);
        formCountCons(model, spq->cons[i], xx.get(), formOptions);
        formSumCons(model, spq->cons[i], xx.get(), formOptions);
    }
    formSumObj(model, spq->obj, xx.get(), formOptions);
    formCntObj(model, spq->obj, xx.get(), formOptions);
    formExpSumObj(model, spq->obj, xx.get(), formOptions);
    model.update();
    return model;
}