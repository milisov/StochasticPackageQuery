#include "Naive.hpp"
#include <boost/algorithm/string/join.hpp>
#include <fmt/ranges.h>

void printNaiveConstraints(GRBModel &model)
{
    GRBConstr *constrs = model.getConstrs();
    int numConstrs = model.get(GRB_IntAttr_NumConstrs);

    cout << "NUMBER OF CONSTRAINTS ADDED: " << numConstrs << endl;
    for (int i = 0; i < numConstrs; i++)
    {
        int maxCoef = 0;
        int minCoef = 9999999;
        cout << endl;
        GRBConstr constr = constrs[i];
        string constrName = constr.get(GRB_StringAttr_ConstrName);
        GRBLinExpr expr = model.getRow(constr);
        double rhs = constr.get(GRB_DoubleAttr_RHS);
        char sense = constr.get(GRB_CharAttr_Sense);

        cout << "Constraint " << constrName << ": ";
        for (int j = 0; j < expr.size(); j++)
        {
            GRBVar var = expr.getVar(j);
            double coeff = expr.getCoeff(j);
            if (coeff > maxCoef)
            {
                maxCoef = coeff;
            }
            if (coeff < minCoef)
            {
                minCoef = coeff;
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
        cout << "MAX =" << maxCoef << " MIN = " << minCoef << endl;
    }
    int numGenConstrs = model.get(GRB_IntAttr_NumGenConstrs);
    cout << "NUMBER OF GENERAL CONSTRAINTS ADDED: " << numGenConstrs << endl;

    GRBGenConstr *genConstrs = model.getGenConstrs();
    for (int i = 0; i < 3; i++)
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

void printNaiveObjective(GRBModel &model)
{
    try
    {
        int sense = model.get(GRB_IntAttr_ModelSense);
        std::string senseStr = (sense == GRB_MINIMIZE) ? "Minimize" : "Maximize";

        // Print the sense of the objective
        std::cout << senseStr << " ";

        // Get the linear part of the objective
        GRBQuadExpr quadExpr = model.getObjective();
        GRBLinExpr linearObj = quadExpr.getLinExpr();

        int numVars = linearObj.size();
        bool isFirstTerm = true;

        // Print linear terms
        for (int i = 0; i < numVars; ++i)
        {
            GRBVar var = linearObj.getVar(i);
            double coeff = linearObj.getCoeff(i);

            if (coeff != 0.0)
            {
                if (!isFirstTerm && coeff > 0)
                {
                    std::cout << " + ";
                }
                else if (coeff < 0)
                {
                    std::cout << " - ";
                    coeff = -coeff;
                }

                if (coeff != 1.0)
                {
                    std::cout << coeff << "*";
                }

                std::cout << var.get(GRB_StringAttr_VarName);
                isFirstTerm = false;
            }
        }

        std::cout << std::endl;
    }
    catch (GRBException &e)
    {
        std::cout << "Error code 2 = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
}