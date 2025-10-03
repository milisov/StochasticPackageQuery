#include "ssformulator.hpp"
#include <fmt/ranges.h>
#include <boost/algorithm/string/join.hpp>


void printConstraints(GRBModel &model)
{
    GRBConstr *constrs = model.getConstrs();
    int numConstrs = model.get(GRB_IntAttr_NumConstrs);

    cout << "NUMBER OF CONSTRAINTS ADDED: " << numConstrs << endl;

    // for (int i = 0; i < numConstrs; i++)
    // {
    //     int maxCoef = 0;
    //     int minCoef = 9999999;
    //     cout << endl;
    //     GRBConstr constr = constrs[i];
    //     string constrName = constr.get(GRB_StringAttr_ConstrName);
    //     GRBLinExpr expr = model.getRow(constr);
    //     double rhs = constr.get(GRB_DoubleAttr_RHS);
    //     char sense = constr.get(GRB_CharAttr_Sense);

    //     cout << "Constraint " << constrName << ": ";
    //     for (int j = 0; j < expr.size(); j++)
    //     {
    //         GRBVar var = expr.getVar(j);
    //         double coeff = expr.getCoeff(j);
    //         if (coeff > maxCoef)
    //         {
    //             maxCoef = coeff;
    //         }
    //         if (coeff < minCoef)
    //         {
    //             minCoef = coeff;
    //         }
    //         cout << coeff << "*" << var.get(GRB_StringAttr_VarName);
    //         if (j < expr.size() - 1)
    //         {
    //             cout << " + ";
    //         }
    //     }
    //     cout << " ";

    //     switch (sense)
    //     {
    //     case GRB_LESS_EQUAL:
    //         cout << "<= ";
    //         break;
    //     case GRB_EQUAL:
    //         cout << "= ";
    //         break;
    //     case GRB_GREATER_EQUAL:
    //         cout << ">= ";
    //         break;
    //     }

    //     cout << rhs << endl;
    //     cout << "MAX =" << maxCoef << " MIN = " << minCoef << endl;
    // }
    int numGenConstrs = model.get(GRB_IntAttr_NumGenConstrs);
    cout << "NUMBER OF GENERAL CONSTRAINTS ADDED: " << numGenConstrs << endl;

    GRBGenConstr *genConstrs = model.getGenConstrs();
    for (int i = 0; i < numGenConstrs; i++)
    {
        int minGen = 999999;
        int maxGen = 0;
        GRBGenConstr genc = genConstrs[i];
        int type = genc.get(GRB_IntAttr_GenConstrType);
        cout << endl;
        if (type == GRB_GENCONSTR_INDICATOR)
        {
            try
            {
                GRBVar binvar;
                int binval;
                GRBLinExpr expr;
                char sense;
                double rhs;

                model.getGenConstrIndicator(genc, &binvar, &binval, &expr, &sense, &rhs);

                cout << "Indicator Constraint " << genc.get(GRB_StringAttr_GenConstrName) << ": ";
                cout << binvar.get(GRB_StringAttr_VarName) << " == " << binval << " -> ";

                for (int j = 0; j < expr.size(); j++)
                {
                    GRBVar var = expr.getVar(j);
                    double coeff = expr.getCoeff(j);
                    if (coeff > maxGen)
                    {
                        maxGen = coeff;
                    }
                    if (coeff < minGen)
                    {
                        minGen = coeff;
                    }
                    cout << coeff << "*" << var.get(GRB_StringAttr_VarName);
                    if (j < expr.size() - 1)
                    {
                        cout << " + ";
                    }
                }
                cout << " ";

                switch (sense)
                {
                case GRB_LESS_EQUAL:
                    cout << "<= ";
                    break;
                case GRB_EQUAL:
                    cout << "= ";
                    break;
                case GRB_GREATER_EQUAL:
                    cout << ">= ";
                    break;
                }

                cout << rhs << endl;
            }
            catch (GRBException &e)
            {
                cout << "Error code 1 = " << e.getErrorCode() << endl;
                cout << e.getMessage() << endl;
            }
            cout << "GENERATOR MAX =" << maxGen << "MIN = " << minGen << endl;
        }
    }
}


// go through the constraints and delete them
void removeProbConstr(GRBModel &model, std::vector<GRBVar> &yy, std::vector<GRBGenConstr> &genCon, std::vector<GRBConstr> &sumyCon)
{
    for (int i = 0; i < genCon.size(); i++)
    {
        // cout << "deleting Indicator" << endl;
        model.remove(genCon[i]);
    }
    model.update();
    for (int i = 0; i < sumyCon.size(); i++)
    {
        // cout << "deleting SUM yy[i]" << endl;
        model.remove(sumyCon[i]);
    }
    model.update();
    for (int i = 0; i < yy.size(); i++)
    {
        // cout << "deleting a yy[i]" << endl;
        model.remove(yy[i]);
    }
    model.update();
    yy.clear();
    genCon.clear();
    sumyCon.clear();
}


void SSFormulator::formProbCons(GRBModel &model,
                                std::shared_ptr<Constraint> cons,
                                GRBVar *xx,
                                std::vector<std::vector<std::vector<double>>> &summaries,
                                int &probConOrder, FormulateOptions &formOptions)
{
    shared_ptr<ProbConstraint> probCon;
    shared_ptr<AttrConstraint> attrCon;

    bool isstoch = isStochastic(cons, probCon, attrCon);
    if (!isstoch)
    {
        return;
    }
    removeProbConstr(model, yy, genCon, sumyCon);
    model.update();
    // calculate the number of yk indicator variables
    int Z = summaries[probConOrder].size();
    double v = spq->getValue(probCon->v);
    double p = spq->getValue(probCon->p);
    double pZ = ceil(p * Z);
    GRBVar y[Z];
    double coeff_Y[Z];
    for (int z = 0; z < Z; z++)
    {
        y[z] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y[" + to_string(z) + "]");
        yy.push_back(y[z]);
        coeff_Y[z] = 1;
    }
    model.update();
    GRBLinExpr sumYz;
    std::vector<std::vector<double>> S = summaries[probConOrder];
    cout << probConOrder << " " << summaries.size() << endl;
    for (int z = 0; z < Z; z++)
    {
        GRBLinExpr innerCons;
        int sz = S[z].size();
        if (formOptions.reduced)
        {
            for (int i = 0; i < formOptions.reducedIds.size(); i++)
            {
                    int id = formOptions.reducedIds[i] - 1;
                    double coeffVal = S[z][id];
                    //innerCons += coeffVal * xx[id];
                    innerCons += coeffVal * xx[i];
            }
        }
        else
        {
            for (int i = 0; i < sz; i++)
            {
                int id = i;
                double coeffVal = S[z][id];
                innerCons += coeffVal * xx[id];
            }
        }

        try
        {
            if (probCon->vsign == Inequality::gteq)
            {
                GRBGenConstr indicator = model.addGenConstrIndicator(y[z], 1, innerCons, GRB_GREATER_EQUAL, v);
                genCon.push_back(indicator);
            }
            else
            {
                GRBGenConstr indicator = model.addGenConstrIndicator(y[z], 1, innerCons, GRB_LESS_EQUAL, v);
                genCon.push_back(indicator);
            }
        }
        catch (GRBException &e)
        {
            cout << "Error code 8 = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        }
    }

    sumYz.addTerms(coeff_Y, y, Z);
    try
    {
        if (probCon->psign == Inequality::gteq)
        {
            GRBConstr constr = model.addConstr(sumYz, GRB_GREATER_EQUAL, pZ);
            sumyCon.push_back(constr);
        }
        else
        {
            GRBConstr constr = model.addConstr(sumYz, GRB_LESS_EQUAL, pZ);
            sumyCon.push_back(constr);
        }
    }
    catch (GRBException &e)
    {
        cout << "Error code 9 = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    probConOrder += 1;
    model.update();
}

//formulate is used for the most basic forms of the algorithm, so this function shall not be used for SDR
GRBModel SSFormulator::formulate(std::shared_ptr<StochasticPackageQuery> spq,FormulateOptions &formOptions)
{
    if(formOptions.iteration == 0)
    {
        DecisionVarOptions decVarOptions = formOptions.decisionVarOptions;
        std::unique_ptr<GRBVar[]> xx;
        if(formOptions.reduced)
        {
            xx = std::make_unique<GRBVar[]>(formOptions.reducedIds.size());
            for (int i = 0; i < formOptions.reducedIds.size(); i++)
            {
                decVarOptions.name = "xx[" + to_string(i) + "]";
                xx[i] = addDecisionVar(model,decVarOptions);
            }
        }else
        {
            xx = std::make_unique<GRBVar[]>(NTuples);
            for (int i = 0; i < NTuples; i++)
            {
                decVarOptions.name = "xx[" + to_string(i) + "]";
                xx[i] = addDecisionVar(model,decVarOptions);
            }
        }
        model.update();
    }

    std::vector<std::vector<std::vector<double>>> summaries;
    if(formOptions.iteration != 0)
    {
        int conOrder = 0;
        for (int i = 0; i < spq->cons.size(); i++)
        {
            shared_ptr<ProbConstraint> probCon;
            shared_ptr<AttrConstraint> attrCon;
            bool isstoch = isStochastic(spq->cons[i], probCon, attrCon);
            if (isstoch)
            {
                std::vector<std::vector<double>> summariesCons = summarize(formOptions, probCon, attrCon, conOrder);
                summaries.push_back(summariesCons);
                conOrder++;
            }
        }   
    }
    
    GRBVar* xx = model.getVars();;
    formulateSAA(model, summaries, xx, formOptions);
    //printConstraints(model);
    return model;
}

void SSFormulator::formulateSAA(GRBModel &model,std::vector<std::vector<std::vector<double>>> &summaries,
                                GRBVar *xx,FormulateOptions &formOptions) //q here is iteration number from the paper
{
    int numCons = spq->cons.size();
    int probConOrder = 0;
    for (int i = 0; i < numCons; i++)
    {
        if (formOptions.iteration == 0)
        {
            formCountCons(model, spq->cons[i], xx, formOptions);
            formSumCons(model, spq->cons[i], xx, formOptions);
            formExpCons(model, spq->cons[i], xx, formOptions);
        }
        else
        {
            formProbCons(model, spq->cons[i], xx, summaries, probConOrder, formOptions);
            GRBConstr *constrs = model.getConstrs();
        }
    }
    if (formOptions.iteration == 0)
    {
        formSumObj(model, spq->obj, xx, formOptions);
        formExpSumObj(model, spq->obj, xx, formOptions);
        formCntObj(model, spq->obj, xx, formOptions);
    }
    model.update();
}

GRBModel SSFormulator::formulateDeterministicLP(std::shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions)
{
    cout<<"Variables Set"<<endl;
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
        //formCountCons(model, spq->cons[i], xx, formOptions); REMOVE THE COUNT CONSTRAINT
        formSumCons(modelDet, spq->cons[i], xx.get(), formOptions);
    }
    formSumObj(modelDet, spq->obj, xx.get(), formOptions);
    formCntObj(modelDet, spq->obj, xx.get(), formOptions);
    formExpSumObj(modelDet , spq->obj, xx.get(), formOptions);
    modelDet.update();
    return modelDet;
}

GRBModel SSFormulator::formulateCVaRLP(std::shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions)
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


void SSFormulator::updateBound(GRBModel& model, double ub) {
    try {
        // Get all variables in the model
        int numVars = model.get(GRB_IntAttr_NumVars);
        GRBVar* vars = model.getVars();

        // Update upper bound for each variable
        for (int i = 0; i < numVars; ++i) {
            vars[i].set(GRB_DoubleAttr_UB, ub);
        }

        // Apply changes
        model.update();

    } catch (GRBException& e) {
        std::cerr << "Gurobi error code = " << e.getErrorCode()
                  << "\n" << e.getMessage() << std::endl;
    } catch (std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
}