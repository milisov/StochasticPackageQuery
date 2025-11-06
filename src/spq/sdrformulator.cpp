#include "sdrformulator.hpp"


GRBModel SDRFormulator::formulateDeterministicLP(std::shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions)
{
    DecisionVarOptions decVarOptions = formOptions.decisionVarOptions;
    std::unique_ptr<GRBVar[]> xx;
    xx = std::make_unique<GRBVar[]>(NTuples);
    for (int i = 0; i < NTuples; i++)
    {
        decVarOptions.name = "xxDet[" + to_string(i) + "]";
        xx[i] = addDecisionVar(modelDet,decVarOptions);
    }
    modelDet.update();
    cout << "Formulating Deterministic LP" << endl;
    int numCons = spq->cons.size();
    int probConOrder = 0;
    for (int i = 0; i < numCons; i++)
    {
        formCountCons(modelDet, spq->cons[i], xx.get(), formOptions); //REMOVE THE COUNT CONSTRAINT
        formSumCons(modelDet, spq->cons[i], xx.get(), formOptions);
    }
    formSumObj(modelDet, spq->obj, xx.get(), formOptions);
    formCntObj(modelDet, spq->obj, xx.get(), formOptions);
    formExpSumObj(modelDet , spq->obj, xx.get(), formOptions);
    modelDet.update();
    return modelDet;
}

GRBModel SDRFormulator::formulateCVaRLP(std::shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions)
{
    cout<<"Variables Set"<<endl;
    DecisionVarOptions decVarOptions = formOptions.decisionVarOptions;
    std::unique_ptr<GRBVar[]> xx;
    xx = std::make_unique<GRBVar[]>(NTuples);
    for (int i = 0; i < NTuples; i++)
    {
        decVarOptions.name = "xxCVaR[" + to_string(i) + "]";
        xx[i] = addDecisionVar(modelCVaR,decVarOptions);
    }
    modelCVaR.update();
    cout << "Formulating CVaR LP" << endl;
    vector<int> reducedIds;
    int numCons = spq->cons.size();
    int probConOrder = 0;
    for (int i = 0; i < numCons; i++)
    {
        formLCVaR(modelCVaR, spq->cons[i], xx.get(), formOptions);
    }
    formSumObj(modelCVaR, spq->obj, xx.get(), formOptions);
    formCntObj(modelCVaR, spq->obj, xx.get(), formOptions);
    formExpSumObj(modelCVaR, spq->obj, xx.get(), formOptions);
    modelCVaR.update();
    return modelCVaR;
}


void SDRFormulator::updateBound(GRBModel& model, double ub) {
    try {
        // Get all variables in the model
        int numVars = model.get(GRB_IntAttr_NumVars);
        GRBVar* vars = model.getVars();

        // Update upper bound for each variable
        for (int i = 0; i < numVars; ++i) {
            vars[i].set(GRB_DoubleAttr_UB, ub);
        }

        model.update();

    } catch (GRBException& e) {
        std::cerr << "Gurobi error code = " << e.getErrorCode()
                  << "\n" << e.getMessage() << std::endl;
    } catch (std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
}

GRBModel SDRFormulator::formulate(std::shared_ptr<StochasticPackageQuery> spq,
                                  FormulateOptions& formOptions)
{
    DecisionVarOptions decVarOptions = formOptions.decisionVarOptions;
    std::unique_ptr<GRBVar[]> xx;
    xx = std::make_unique<GRBVar[]>(NTuples);
    for (int i = 0; i < NTuples; i++)
    {
        decVarOptions.name = "xx[" + to_string(i) + "]";
        xx[i] = addDecisionVar(model,decVarOptions);
    }
    model.update();
    int numCons = spq->cons.size();
    int probConOrder = 0;
    for (int i = 0; i < numCons; i++)
    {
        formLCVaR(model, spq->cons[i], xx.get(), formOptions);
        //formCountCons(model, spq->cons[i], xx.get(), formOptions); //REMOVE THE COUNT CONSTRAINT
        formSumCons(model, spq->cons[i], xx.get(), formOptions);
    }
    formSumObj(model, spq->obj, xx.get(), formOptions);
    formCntObj(model, spq->obj, xx.get(), formOptions);
    formExpSumObj(model, spq->obj, xx.get(), formOptions);
    model.update();
    return model;
    
}