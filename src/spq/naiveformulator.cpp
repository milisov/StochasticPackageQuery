#include "naiveformulator.hpp"
#include <fmt/ranges.h>
#include <boost/algorithm/string/join.hpp>

using namespace std;


void NaiveFormulator::formProbCons(GRBModel &model, std::shared_ptr<Constraint> cons, GRBVar *xx, FormulateOptions& options) 
{
    deb("formulating Prob Cons here");
    shared_ptr<ProbConstraint> probCon;
    shared_ptr<AttrConstraint> attrCon;

    bool isstoch = isStochastic(cons, probCon, attrCon);
    if (!isstoch)
    {
        return;
    }

    double v = spq->getValue(probCon->v);
    //double p = spq->getValue(probCon->p);
    double p = options.p;
    double pM = ceil(p * cntScenarios);
    GRBVar y[cntScenarios];
    for (int j = 0; j < cntScenarios; j++)
    {
        y[j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y[" + to_string(j) + "]");
    }
    model.update();

    GRBLinExpr sumYz;
    vector<double>vect;
    initializeVectorForm(vect, NTuples, (double) 0.0);


    std::vector<GRBLinExpr> innerCons(cntScenarios);
    auto& scenarios = data.stochAttrs[attrCon->attr];
    if(options.reduced)
    {
        for(int i = 0; i < options.reducedIds.size(); i++)
        {
            int id = options.reducedIds[i] - 1;
            for(int j = 0; j < cntScenarios; j++)
            {
                innerCons[j] += xx[i] * scenarios[id][j];
            }
        }
    }else
    {
        for(int i = 0; i < NTuples; i++)
        {
            for(int j = 0; j < cntScenarios; j++)
            {
                innerCons[j] += xx[i] * scenarios[i][j];
            }
        }  
    }
    for(int i = 0; i < innerCons.size(); i++)
    {
        try
        {
            if (probCon->vsign == Inequality::gteq)
            {
                GRBGenConstr indicator = model.addGenConstrIndicator(y[i], 1, innerCons[i], GRB_GREATER_EQUAL, v);
            }
            else
            {
                GRBGenConstr indicator = model.addGenConstrIndicator(y[i], 1, innerCons[i], GRB_LESS_EQUAL, v);

            }
        }
        catch (GRBException &e)
        {
            deb("Error code 8 = ",e.getErrorCode());
            cout << e.getMessage() << endl;
        }
    }

    for (int j = 0; j < cntScenarios; j++)
    {
        sumYz += y[j];
    }
    try
    {
        if (probCon->psign == Inequality::gteq)
        {
            GRBConstr constr = model.addConstr(sumYz, GRB_GREATER_EQUAL, pM);
        }
        else
        {
            GRBConstr constr = model.addConstr(sumYz, GRB_LESS_EQUAL, pM);
        }
    }
    catch (GRBException &e)
    {
        cout << "Error code 9 = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
}



void NaiveFormulator::formProbConsActiveness(GRBModel &model,
                            std::shared_ptr<Constraint> cons,
                            GRBVar *xx,
                            FormulateOptions& options)
{
    shared_ptr<ProbConstraint> probCon;
    shared_ptr<AttrConstraint> attrCon;

    bool isstoch = isStochastic(cons, probCon, attrCon);
    if (!isstoch)
    {
        return;
    }

    sort(options.posActiveness.begin(), options.posActiveness.end(),
         [](const pair<int, double> &a, const pair<int, double> &b)
         {
             return a.second < b.second;
         });

    sort(options.negActiveness.begin(), options.negActiveness.end(),
         [](const pair<int, double> &a, const pair<int, double> &b)
         {
             return a.second > b.second;
         });

    double v = spq->getValue(probCon->v);
    //double p = spq->getValue(probCon->p);
    double p = options.p;

    int posTarget = 0;
    int negTarget = 0; 
    int scDimension = min(options.qSz,cntScenarios);
    if(options.posActiveness.size() >= ((int)p*cntScenarios))
    {
        posTarget = scDimension;
    }else
    {
        int negativeNeeded = ((int)p*cntScenarios) - options.posActiveness.size();
        negTarget = min(scDimension, negativeNeeded);
        posTarget = scDimension - negTarget;
    }
    std::vector<GRBLinExpr> innerCons(scDimension);
    auto &scenarios = data.stochAttrs[attrCon->attr];
    int q;
    for(int i = 0; i < options.reducedIds.size(); i++)
    {
        q = 0;
        for (int j = 0; j < options.negActiveness.size(); j++)
        {
            if (q == negTarget)
            {
                break;
            }
            int id = options.reducedIds[i] - 1;
            int scenarioId = options.negActiveness[j].first;
            innerCons[q] += xx[i] * scenarios[id][scenarioId];
            q++;
        }
    }

    for (int i = 0; i < options.reducedIds.size(); i++)
    {
        q = 0;
        for (int j = 0; j < options.posActiveness.size(); j++)
        {
            if (q == posTarget)
            {
                break;
            }
            int id = options.reducedIds[i] - 1;
            int scenarioId = options.posActiveness[j].first;
            innerCons[q] += xx[i] * scenarios[id][scenarioId];
            q++;
        }
    }
    try
    {
        if (probCon->psign == Inequality::gteq)
        {
            for(int i = 0; i < options.qSz; i++)
            {
                GRBConstr constr = model.addConstr(innerCons[i], GRB_GREATER_EQUAL, v);
            }
        }
        else
        {
            for(int i = 0; i < options.qSz; i++)
            {
                GRBConstr constr = model.addConstr(innerCons[i], GRB_LESS_EQUAL, v);
            }
        }
    }
    catch (GRBException &e)
    {
        cout << "Error code 9 = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
}




GRBModel NaiveFormulator::formulate(shared_ptr<StochasticPackageQuery> spq, FormulateOptions& formOptions) 
{
    GRBModel model(env);
    std::unique_ptr<GRBVar[]> xx;
    if(formOptions.reduced)
    {
        xx = std::make_unique<GRBVar[]>(formOptions.reducedIds.size());
        DecisionVarOptions decVarOptions = formOptions.decisionVarOptions;
        for (int i = 0; i < formOptions.reducedIds.size(); i++)
        {
            decVarOptions.name = "xx[" + to_string(i) + "]";
            xx[i] = addDecisionVar(model, decVarOptions);
        }
    }else
    {
        xx = std::make_unique<GRBVar[]>(NTuples);
        DecisionVarOptions decVarOptions = formOptions.decisionVarOptions;
        for (int i = 0; i < NTuples; i++)
        {
            decVarOptions.name = "xx[" + to_string(i) + "]";
            xx[i] = addDecisionVar(model,decVarOptions);
        }
    }
    int numCons = spq->cons.size();

    for (int i = 0; i < numCons; i++)
    {
        formCountCons(model, spq->cons[i], xx.get(), formOptions);
        formSumCons(model, spq->cons[i], xx.get(), formOptions);
        if(formOptions.reducedScenarios)
        {
            formProbConsActiveness(model, spq->cons[i], xx.get(), formOptions);
        }else
        {
            formProbCons(model, spq->cons[i], xx.get(), formOptions);
        }
        formExpCons(model, spq->cons[i], xx.get(), formOptions);
    }
    formSumObj(model, spq->obj, xx.get(), formOptions);
    formExpSumObj(model, spq->obj, xx.get(), formOptions);
    formCntObj(model, spq->obj, xx.get(), formOptions);
    deb("formulated");
    model.update();
    return model;
}