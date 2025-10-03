#include "rsformulator.hpp"
#include <fmt/ranges.h>
#include <fmt/ranges.h>
#include <boost/algorithm/string/join.hpp>
#include <iostream>

void printObjCons(GRBModel &model)
{
    try
    {
        GRBConstr objConstr = model.getConstrByName("objcons");

        // Get attributes
        std::string cname = objConstr.get(GRB_StringAttr_ConstrName);
        char sense = objConstr.get(GRB_CharAttr_Sense);
        double rhs = objConstr.get(GRB_DoubleAttr_RHS);

        std::cout << "Constraint name: " << cname << std::endl;
        std::cout << "Sense: " << sense << " ("
                  << (sense == GRB_LESS_EQUAL ? "<=" : sense == GRB_GREATER_EQUAL ? ">="
                                                                                  : "=")
                  << ")" << std::endl;
        std::cout << "RHS: " << rhs << std::endl;

        // Left-hand side
        GRBLinExpr lhs = model.getRow(objConstr);
        std::cout << "LHS:";
        if (lhs.getConstant() != 0.0)
            std::cout << " " << lhs.getConstant();

        for (int i = 0; i < 5; i++)
        {
            GRBVar v = lhs.getVar(i);
            double coeff = lhs.getCoeff(i);
            std::cout << " + (" << coeff << " * "
                      << v.get(GRB_StringAttr_VarName) << ")";
        }
        std::cout << std::endl;
    }
    catch (GRBException &e)
    {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
    catch (...)
    {
        std::cout << "Constraint 'objcons' not found." << std::endl;
    }
}

void printMinMaxModel(GRBModel &model)
{
    try
    {
        // Print all constraints with name containing "minmax_"
        int numConstrs = model.get(GRB_IntAttr_NumConstrs);
        std::cout << "MinMax Constraints:" << std::endl;
        int cnt = 0;
        for (int i = 0; i < numConstrs; i++)
        {
            GRBConstr constr = model.getConstr(i);
            std::string name = constr.get(GRB_StringAttr_ConstrName);
            if (name.find("minmax_") != std::string::npos)
            {
                cnt++;
                // GRBLinExpr lhs = model.getRow(constr);
                // double rhs = constr.get(GRB_DoubleAttr_RHS);
                // char sense = constr.get(GRB_CharAttr_Sense); // '<', '>', '='

                // std::cout << "  " << name << ": ";

                // // Print LHS
                // for (int j = 0; j < lhs.size(); j++)
                // {
                //     GRBVar var = lhs.getVar(j);
                //     double coeff = lhs.getCoeff(j);
                //     std::cout << coeff << "*" << var.get(GRB_StringAttr_VarName);
                //     if (j < lhs.size() - 1)
                //         std::cout << " + ";
                // }

                // Print sense and RHS
                //  std::cout << " " << sense << " " << rhs << std::endl;
            }
        }
    }
    catch (GRBException &e)
    {
        std::cerr << "Gurobi exception: " << e.getMessage() << std::endl;
    }
    catch (...)
    {
        std::cerr << "Unknown exception while printing model" << std::endl;
    }
}

void printCVaRConstraintsFromModel(GRBModel &model)
{
    int numConstrs = model.get(GRB_IntAttr_NumConstrs);
    GRBConstr *constrs = model.getConstrs();

    for (int k = 0; k < numConstrs; k++)
    {
        GRBLinExpr expr = model.getRow(constrs[k]);
        std::string name = constrs[k].get(GRB_StringAttr_ConstrName);
        if (name[0] == 'm' || name[0] == 'R')
        {
            continue;
        }
        std::cout << name << ": ";
        int numTerms = expr.size();
        for (int i = 0; i < numTerms; i++)
        {
            double coeff = expr.getCoeff(i);
            std::string varName = expr.getVar(i).get(GRB_StringAttr_VarName);

            if (i > 0 && coeff >= 0)
                std::cout << " + ";
            else if (coeff < 0)
                std::cout << " - ";

            std::cout << std::abs(coeff) << "*" << varName;
        }

        char sense = constrs[k].get(GRB_CharAttr_Sense);
        double rhs = constrs[k].get(GRB_DoubleAttr_RHS);

        switch (sense)
        {
        case GRB_LESS_EQUAL:
            std::cout << " <= ";
            break;
        case GRB_EQUAL:
            std::cout << " = ";
            break;
        case GRB_GREATER_EQUAL:
            std::cout << " >= ";
            break;
        }
        std::cout << rhs << std::endl;
    }

    delete[] constrs;
}

void RSFormulator::removeLastProbConstraints(GRBModel &model)
{
    // Remove all y[...] variables
    int numVars = model.get(GRB_IntAttr_NumVars);
    GRBVar *vars = model.getVars();
    for (int i = 0; i < numVars; i++)
    {
        std::string name = vars[i].get(GRB_StringAttr_VarName);
        if (name.rfind("y[", 0) == 0) // starts with "y["
        {
            model.remove(vars[i]);
        }
        else if (name == "beta")
        {
            model.remove(vars[i]);
        }
    }
    delete[] vars;

    // Remove all y[...] constraints and the cvar constraint
    int numConstrs = model.get(GRB_IntAttr_NumConstrs);
    GRBConstr *constrs = model.getConstrs();
    for (int i = 0; i < numConstrs; i++)
    {
        std::string name = constrs[i].get(GRB_StringAttr_ConstrName);

        // Check if constraint name starts with "y[" and ends with "Cons"
        if (name.rfind("y[", 0) == 0 &&
            name.size() >= 4 &&
            name.substr(name.size() - 4) == "Cons")
        {
            model.remove(constrs[i]);
        }
        // Remove the cvar constraint explicitly
        else if (name == "cvar")
        {
            model.remove(constrs[i]);
        }
    }
    delete[] constrs;

    model.update();
}

void RSFormulator::formulateRS(GRBModel &model, std::shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions)
{
    std::vector<std::vector<std::vector<double>>> summaries;
    int numCons = spq->cons.size();
    int probConOrder = 0;
    int conOrder = 0;

    if (formOptions.SAA == false)
    {
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
        for (int i = 0; i < spq->cons.size(); i++)
        {
            // checks inside if stochastic or not, probConOrder is increased within the function too
            formProbCons(model, spq->cons[i], model.getVars(), summaries, probConOrder, formOptions);
        }
    }
    else
    {
        for (int i = 0; i < spq->cons.size(); i++)
        {
            // checks inside if stochastic or not, probConOrder is increased within the function too
            formProbConsSAA(model, spq->cons[i], model.getVars(), probConOrder, formOptions);
        }
    }
    model.update();
}

void RSFormulator::formulateBestObjProblem(GRBModel &model, std::shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions)
{
    std::unique_ptr<GRBVar[]> xx;
    xx = std::make_unique<GRBVar[]>(NTuples);
    DecisionVarOptions decVarOptions = formOptions.decisionVarOptions;
    for (int i = 0; i < NTuples; i++)
    {
        decVarOptions.name = "xx[" + to_string(i) + "]";
        xx[i] = addDecisionVar(model, decVarOptions);
    }

    int numCons = spq->cons.size();

    for (int i = 0; i < numCons; i++)
    {
        formCountCons(model, spq->cons[i], xx.get(), formOptions);
        formSumCons(model, spq->cons[i], xx.get(), formOptions);
        formExpCons(model, spq->cons[i], xx.get(), formOptions);
    }

    formCntObj(model, spq->obj, xx.get(), formOptions);
    formSumObj(model, spq->obj, xx.get(), formOptions);
    formExpSumObj(model, spq->obj, xx.get(), formOptions);
    model.update();
}

void RSFormulator::formulateDeterministic(GRBModel &model, std::shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions)
{
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
            xx[i] = addDecisionVar(model, decVarOptions);
        }
    }

    int numCons = spq->cons.size();

    for (int i = 0; i < numCons; i++)
    {
        formCountCons(model, spq->cons[i], xx.get(), formOptions);
        formSumCons(model, spq->cons[i], xx.get(), formOptions);
        formExpCons(model, spq->cons[i], xx.get(), formOptions);
    }
    formObjCons(model, spq->obj, xx.get(), formOptions);
    formMinMaxObjective(model, xx.get(), formOptions);
    model.update();
}

void RSFormulator::formProbConsSAA(GRBModel &model, std::shared_ptr<Constraint> cons,
                                   GRBVar *xx, int &probConOrder, FormulateOptions &formOptions)
{
    shared_ptr<ProbConstraint> probCon;
    shared_ptr<AttrConstraint> attrCon;

    bool isstoch = isStochastic(cons, probCon, attrCon);
    if (!isstoch)
    {
        return;
    }
    int numCons = model.get(GRB_IntAttr_NumConstrs);
    int numVars = model.get(GRB_IntAttr_NumVars);
    removeLastProbConstraints(model);
    model.update();
    // calculate the number of yk indicator variables
    int Z = cntScenarios;
    double v = spq->getValue(probCon->v);
    double p = spq->getValue(probCon->p);
    double pZ = ceil(p * Z);
    GRBVar y[Z];

    GRBVar beta = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "beta");

    for (int i = 0; i < Z; i++)
    {
        y[i] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "y[" + to_string(i) + "]");
    }

    auto& scenarios = data.stochAttrs[attrCon->attr];
    std::vector<GRBLinExpr> innerCons(Z);
    if(formOptions.reduced)
    {
        for(int i = 0; i < formOptions.reducedIds.size(); i++)
        {
            int id = formOptions.reducedIds[i] - 1; 
            for(int j = 0; j < cntScenarios; j++)
            {
                double coefVal = -scenarios[id][j];
                //innerCons[j] += coefVal * xx[id];
                innerCons[j] += coefVal * xx[i];
            }
        }
    }
    else
    {
        for(int i = 0; i < NTuples; i++)
        {
            for(int j = 0; j < cntScenarios; j++)
            {
                double coefVal = -scenarios[i][j];
                innerCons[j] += coefVal * xx[i];
            }
        }
        
    }
    for (int z = 0; z < Z; z++)
    {
        innerCons[z] += beta;
        model.addConstr(y[z] >= innerCons[z], "y[" + std::to_string(z) + "]Cons");
    }

    GRBLinExpr cvar = beta;
    GRBLinExpr sumYs = 0.0;
    for (int i = 0; i < Z; i++)
    {
        sumYs += y[i];
    }
    cvar += (-1.0 / (Z * p)) * sumYs;
    model.addConstr(cvar >= v, "cvar");
    model.update();
    probConOrder += 1;
}

void RSFormulator::formProbCons(GRBModel &model,
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
    int numCons = model.get(GRB_IntAttr_NumConstrs);
    int numVars = model.get(GRB_IntAttr_NumVars);
    removeLastProbConstraints(model);
    model.update();
    // calculate the number of yk indicator variables
    int Z = summaries[probConOrder].size();
    double v = spq->getValue(probCon->v);
    double p = spq->getValue(probCon->p);
    double pZ = ceil(p * Z);
    GRBVar y[Z];

    GRBVar beta = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "beta");

    for (int i = 0; i < Z; i++)
    {
        y[i] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "y[" + to_string(i) + "]");
    }

    std::vector<std::vector<double>> S = summaries[probConOrder];
    for (int z = 0; z < Z; z++)
    {
        GRBLinExpr innerCons;
        int sz = S[z].size();
        innerCons = beta;
        if (formOptions.reduced)
        {
            for (int i = 0; i < formOptions.reducedIds.size(); i++)
            {
                int id = formOptions.reducedIds[i] - 1;
                double coeffVal = -S[z][id];
                //innerCons += coeffVal * xx[id];
                innerCons += coeffVal * xx[i];
            }
        }
        else
        {
            for (int i = 0; i < sz; i++)
            {
                double coeffVal = -S[z][i];
                int id = i;
                innerCons += coeffVal * xx[id];
            }
        }
        model.addConstr(y[z] >= innerCons, "y[" + std::to_string(z) + "]Cons");
    }

    GRBLinExpr cvar = beta;
    GRBLinExpr sumYs = 0.0;
    for (int i = 0; i < Z; i++)
    {
        sumYs += y[i];
    }
    // cvar += (-1.0 / (Z * (1 - p))) * sumYs;
    cvar += (-1.0 / (Z * p)) * sumYs;
    model.addConstr(cvar >= v, "cvar");
    model.update();
    probConOrder += 1;
}

void RSFormulator::formObjCons(GRBModel &model,
                               std::shared_ptr<Objective> obj,
                               GRBVar *xx,
                               FormulateOptions &formOptions)
{
    shared_ptr<AttrObjective> attrObj;
    bool isDet = isDeterministic(obj, attrObj);
    if (!isDet || attrObj->objType != array_type || formOptions.objCons == false)
    {
        return;
    }

    GRBConstr objConstr;
    bool exists = true;
    try
    {
        objConstr = model.getConstrByName("objcons"); // check if already exists
    }
    catch (...)
    {
        exists = false;
    }

    if (exists)
    {
        objConstr.set(GRB_DoubleAttr_RHS, formOptions.objValue); // update RHS
        model.update();
    }
    else
    {
        GRBLinExpr expSumObjExpr;
        if (formOptions.reduced)
        {
            for (int i = 0; i < formOptions.reducedIds.size(); i++)
            {
                int id = formOptions.reducedIds[i] - 1;
                //expSumObjExpr += xx[id] * data.stockExpectedProfit[id];
                expSumObjExpr += xx[i] * data.stockExpectedProfit[id];
            }
        }
        else
        {
            for (int i = 0; i < NTuples; i++)
            {
                int id = i;
                expSumObjExpr += xx[id] * data.stockExpectedProfit[id];
            }
        }
    
        try
        {
            if (obj->objSense == maximize)
            {
                model.addConstr(expSumObjExpr >= formOptions.objValue, "objcons");
            }
            else
            {
                model.addConstr(expSumObjExpr <= formOptions.objValue, "objcons");
            }
        }
        catch (GRBException e)
        {
            cout << "Error code 4 = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        }
    }
}

void RSFormulator::formMinMaxObjective(GRBModel &model, GRBVar *xx, FormulateOptions &formOptions)
{
    // removeMinMaxObjective(model);
    GRBLinExpr minMaxExpr;
    GRBVar k = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "k");
    if (formOptions.reduced)
    {
        // if reduced, we need to use the reduced ids
        for (int i = 0; i < formOptions.reducedIds.size(); i++)
        {
            //int id = formOptions.reducedIds[i] - 1;
            model.addConstr(xx[i] <= k, "minmax_" + std::to_string(i));
        }
    }
    else
    {
        for (int i = 0; i < NTuples; i++)
        {
            // if not reduced, we can use the original index
            model.addConstr(xx[i] <= k, "minmax_" + std::to_string(i));
        }
    }
    minMaxExpr = k;
    model.setObjective(minMaxExpr, GRB_MINIMIZE);
    model.update();
}

GRBModel RSFormulator::formulate(shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions)
{
    if (formOptions.SAA)
    {
        formulateDeterministic(modelRS, spq, formOptions);
        formulateRS(modelRS, spq, formOptions);
    }
    else if (formOptions.iteration == 0 && formOptions.Z == formOptions.Zinit)
    {
        formulateDeterministic(modelRS, spq, formOptions);
        formulateRS(modelRS, spq, formOptions);
    }
    else
    {
        formulateRS(modelRS, spq, formOptions);
    }
    return modelRS;
}