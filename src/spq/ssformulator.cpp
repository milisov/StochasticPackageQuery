#include "ssformulator.hpp"
#include <algorithm>
#include <cmath>
#include <fmt/ranges.h>
#include <boost/algorithm/string/join.hpp>

// go through the constraints and delete them
void removeProbConstr(GRBModel &model, std::vector<GRBVar> &yy, std::vector<GRBGenConstr> &genCon, std::vector<GRBConstr> &probCons)
{
    for (int i = 0; i < genCon.size(); i++)
    {
        // cout << "deleting Indicator" << endl;
        model.remove(genCon[i]);
    }
    model.update();
    for (int i = 0; i < probCons.size(); i++)
    {
        // cout << "deleting SUM yy[i]" << endl;
        model.remove(probCons[i]);
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
    probCons.clear();
}

void SSFormulator::formProbConsStage3(GRBModel &model,
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
    deb("Formulating-ProbCons-Stage-3");
    removeProbConstr(model, yy, genCon, probCons);
    model.update();
    // calculate the number of yk indicator variables
    int Z = summaries[probConOrder].size();
    double v = spq->getValue(probCon->v);
    double p = spq->getValue(probCon->p);
    if(formOptions.changeVaRBound)
    {
        double epsV = 1e-1;
        // double mu = formOptions.meanPerVaR[probConOrder];
        // double var = sqrt(formOptions.variancePerVaR[probConOrder]);
        // double k = sqrt(p/(1-p));
        // deb(mu, k, var);
        // double new_v = mu - k * var;
        // double potential_v = mu + k*var;
        // deb(v, new_v,potential_v);
        deb(v, formOptions.minInnerConstPerVaR[probConOrder]);
        v = formOptions.minInnerConstPerVaR[probConOrder] + epsV;
        deb(v);
    }
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
                // innerCons += coeffVal * xx[id];
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
            if (probCon->psign == Inequality::gteq)
            {
                GRBConstr constr = model.addConstr(innerCons, GRB_GREATER_EQUAL, v);
                probCons.push_back(constr);
            }
            else
            {
                GRBConstr constr = model.addConstr(innerCons, GRB_LESS_EQUAL, v);
                probCons.push_back(constr);
            }
        }
        catch (GRBException &e)
        {
            cout << "Error code 9 = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        }
    }
    probConOrder += 1;
    model.update();
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
    removeProbConstr(model, yy, genCon, probCons);
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
                // innerCons += coeffVal * xx[id];
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
            probCons.push_back(constr);
        }
        else
        {
            GRBConstr constr = model.addConstr(sumYz, GRB_LESS_EQUAL, pZ);
            probCons.push_back(constr);
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

// formulate is used for the most basic forms of the algorithm, so this function shall not be used for SDR
GRBModel SSFormulator::formulate(std::shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions)
{
    if (formOptions.iteration == 0)
    {
        DecisionVarOptions decVarOptions = formOptions.decisionVarOptions;
        std::unique_ptr<GRBVar[]> xx;
        if (formOptions.reduced)
        {
            xx = std::make_unique<GRBVar[]>(formOptions.reducedIds.size());
            for (int i = 0; i < formOptions.reducedIds.size(); i++)
            {
                decVarOptions.name = "xx[" + to_string(i) + "]";
                xx[i] = addDecisionVar(model, decVarOptions);
            }
        }
        else
        {
            xx = std::make_unique<GRBVar[]>(NTuples);
            for (int i = 0; i < NTuples; i++)
            {
                decVarOptions.name = "xx[" + to_string(i) + "]";
                xx[i] = addDecisionVar(model, decVarOptions);
            }
        }
        model.update();
    }

    std::vector<std::vector<std::vector<double>>> summaries;
    if (formOptions.iteration != 0)
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

    GRBVar *xx = model.getVars();
    formulateSAA(model, summaries, xx, formOptions);
    if(formOptions.iteration == 1)
    {
        model.write("/home/fm2288/StochasticPackageQuery/model"+std::to_string(model.get(GRB_IntAttr_NumVars))+"_"+std::to_string(formOptions.iteration)+".lp");
    }
    return model;
}

void SSFormulator::formulateSAA(GRBModel &model, std::vector<std::vector<std::vector<double>>> &summaries,
                                GRBVar *xx, FormulateOptions &formOptions) // q here is iteration number from the paper
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
            if(formOptions.finalStage && !formOptions.indicators)
            {
                formProbConsStage3(model, spq->cons[i], xx, summaries, probConOrder, formOptions);
                if(formOptions.varianceControl)
                {
                    // formVarianceControlConstr(model, spq->cons[i], xx, formOptions);
                    formMADVarianceControlConstr(model, spq->cons[i], xx, formOptions);
                }
            }else
            {
                formProbCons(model, spq->cons[i], xx, summaries, probConOrder, formOptions);
            }
        }
    }
    if (formOptions.iteration == 0)
    {
        formSumObj(model, spq->obj, xx, formOptions);
        formExpSumObj(model, spq->obj, xx, formOptions);
        formCntObj(model, spq->obj, xx, formOptions);
    }
    model.update();
    int numVars = model.get(GRB_IntAttr_NumVars);
    //model.write("/home/fm2288/StochasticPackageQuery/model"+std::to_string(numVars)+"_"+std::to_string(formOptions.iteration)+".lp");
}

// Formulates variance control constraint:
// v <= x^T * mu - sqrt(p/(1-p)) * ||x .* sigma||_2
// Using SOC reformulation with auxiliary variables:
//   1. y_i = sigma_i * x_i  (linear, since sigma_i is constant)
//   2. t = ||y||_2          (norm constraint)
//   3. sum_i(mu_i * x_i) - k*t >= v  (linear)
void SSFormulator::formVarianceControlConstr(GRBModel &model,
                                             std::shared_ptr<Constraint> cons,
                                             GRBVar *xx,
                                             FormulateOptions &formOptions)
{
    shared_ptr<ProbConstraint> probCon;
    shared_ptr<AttrConstraint> attrCon;

    bool isstoch = isStochastic(cons, probCon, attrCon);
    if (!isstoch)
    {
        return;
    }

    // Remove existing variance control constraint if present
    if (hasVarianceConstr)
    {
        removeVarianceControlConstr(model);
    }

    double v = spq->getValue(probCon->v);
    double p = spq->getValue(probCon->p);
    double k = std::sqrt(p / (1.0 - p));

    // Get expected profit (mu) and standard deviation (sigma) from data
    const std::vector<double> &mu = data.stockExpectedProfit;
    const std::vector<double> &sigma = data.stockProfitStdDev;

    // Build linear expression: sum_i(mu_i * x_i)
    GRBLinExpr expectedProfit;

    int n = formOptions.reduced ? formOptions.reducedIds.size() : NTuples;

    // Clear storage vectors
    varianceY.clear();
    varianceYCons.clear();

    // Create auxiliary variables y_i = sigma_i * x_i
    for (int i = 0; i < n; i++)
    {
        GRBVar yi = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "y_" + std::to_string(i));
        varianceY.push_back(yi);
    }
    model.update();

    // Add constraints: y_i = sigma_i * x_i
    if (formOptions.reduced)
    {
        for (int i = 0; i < n; i++)
        {
            int id = formOptions.reducedIds[i] - 1;
            expectedProfit += mu[id] * xx[i];
            GRBConstr con = model.addConstr(varianceY[i] == sigma[id] * xx[i], "y_def_" + std::to_string(i));
            varianceYCons.push_back(con);
        }
    }
    else
    {
        for (int i = 0; i < n; i++)
        {
            expectedProfit += mu[i] * xx[i];
            GRBConstr con = model.addConstr(varianceY[i] == sigma[i] * xx[i], "y_def_" + std::to_string(i));
            varianceYCons.push_back(con);
        }
    }

    // Add auxiliary variable t >= 0 for the norm ||y||_2
    varianceT = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "t_norm");
    model.update();

    // Add norm constraint: t = ||y||_2 using addGenConstrNorm
    try
    {
        varianceNormCon = model.addGenConstrNorm(varianceT, varianceY.data(), n, 2.0, "norm_constraint");
    }
    catch (GRBException &e)
    {
        cout << "Error in norm constraint: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }

    // Linear constraint: sum_i(mu_i * x_i) - k*t >= v
    try
    {
        if (probCon->vsign == Inequality::gteq)
        {
            // v <= x^T * mu - k * ||y||_2  =>  x^T * mu - k*t >= v
            varianceControlCon = model.addConstr(expectedProfit - k * varianceT, GRB_GREATER_EQUAL, v, "variance_control");
        }
        else
        {
            throw("please check this again -- we only implemented the gteq version for now");
        }
    }
    catch (GRBException &e)
    {
        cout << "Error in variance control constraint: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }

    hasVarianceConstr = true;
    model.update();
}

// Remove all variance control constraint components from the model
void SSFormulator::removeVarianceControlConstr(GRBModel &model)
{
    if (!hasVarianceConstr)
    {
        return;
    }

    // Remove the main variance control constraint
    try
    {
        model.remove(varianceControlCon);
    }
    catch (GRBException &e)
    {
        cout << "Error removing variance control constraint: " << e.getMessage() << endl;
    }

    // Remove the norm constraint
    try
    {
        model.remove(varianceNormCon);
    }
    catch (GRBException &e)
    {
        cout << "Error removing norm constraint: " << e.getMessage() << endl;
    }

    // Remove y_i definition constraints
    for (auto &con : varianceYCons)
    {
        try
        {
            model.remove(con);
        }
        catch (GRBException &e)
        {
            cout << "Error removing y constraint: " << e.getMessage() << endl;
        }
    }
    model.update();

    // Remove the t variable
    try
    {
        model.remove(varianceT);
    }
    catch (GRBException &e)
    {
        cout << "Error removing t variable: " << e.getMessage() << endl;
    }

    // Remove y_i variables
    for (auto &var : varianceY)
    {
        try
        {
            model.remove(var);
        }
        catch (GRBException &e)
        {
            cout << "Error removing y variable: " << e.getMessage() << endl;
        }
    }

    model.update();

    // Clear storage
    varianceY.clear();
    varianceYCons.clear();
    hasVarianceConstr = false;
}

// Formulates MAD-based variance control constraint:
// v <= x^T * mu - k * MAD_Z
// Z in this case is the sum(x*profit_i)
// where MAD_Z = (1/n) * sum_k |sum_i x_i * (p_i^(k) - mu_i)|
// Linearized using auxiliary variables d_k >= 0:
//   d_k >= sum_i x_i * (p_i^(k) - mu_i)   (positive deviation)
//   d_k >= -sum_i x_i * (p_i^(k) - mu_i)  (negative deviation)
//   sum_i(mu_i * x_i) - (k/n) * sum_k(d_k) >= v
void SSFormulator::formMADVarianceControlConstr(GRBModel &model,
                                                 std::shared_ptr<Constraint> cons,
                                                 GRBVar *xx,
                                                 FormulateOptions &formOptions)
{
    shared_ptr<ProbConstraint> probCon;
    shared_ptr<AttrConstraint> attrCon;

    bool isstoch = isStochastic(cons, probCon, attrCon);
    if (!isstoch)
    {
        return;
    }

    // Remove existing MAD constraint if present
    if (hasMADConstr)
    {
        removeMADVarianceControlConstr(model);
    }

    double v = spq->getValue(probCon->v);
    double p = spq->getValue(probCon->p);
    double k = std::sqrt(p / (1.0 - p)) * std::sqrt(M_PI / 2.0); // Scaling factor for MAD to approximate std dev

    // Get expected profit (mu) from data
    const std::vector<double> &mu = data.stockExpectedProfit;

    // Get scenario data for the stochastic attribute
    const std::string &attrName = attrCon->attr;
    const std::vector<std::vector<double>> &scenarios = data.stochAttrs[attrName];

    int n = formOptions.reduced ? formOptions.reducedIds.size() : NTuples;

    // Clear storage vectors
    madD.clear();
    madDConsPos.clear();
    madDConsNeg.clear();

    // Create auxiliary variables d_k for each scenario
    for (int scenarioIdx = 0; scenarioIdx < cntScenarios; scenarioIdx++)
    {
        GRBVar dk = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "d_" + std::to_string(scenarioIdx));
        madD.push_back(dk);
    }
    model.update();

    // Build linear expression: sum_i(mu_i * x_i)
    GRBLinExpr expectedProfit;

    if (formOptions.reduced)
    {
        for (int i = 0; i < n; i++)
        {
            int id = formOptions.reducedIds[i] - 1;
            expectedProfit += mu[id] * xx[i];
        }
    }
    else
    {
        for (int i = 0; i < n; i++)
        {
            expectedProfit += mu[i] * xx[i];
        }
    }

    // For each scenario k, add constraints:
    //   d_k >= sum_i x_i * (p_i^(k) - mu_i)
    //   d_k >= -sum_i x_i * (p_i^(k) - mu_i)
    for (int scenarioIdx = 0; scenarioIdx < cntScenarios; scenarioIdx++)
    {
        GRBLinExpr deviation;

        if (formOptions.reduced)
        {
            for (int i = 0; i < n; i++)
            {
                int id = formOptions.reducedIds[i] - 1;
                double diff = scenarios[id][scenarioIdx] - mu[id];
                deviation += diff * xx[i];
            }
        }
        else
        {
            for (int i = 0; i < n; i++)
            {
                double diff = scenarios[i][scenarioIdx] - mu[i];
                deviation += diff * xx[i];
            }
        }

        try
        {
            // d_k >= deviation
            GRBConstr conPos = model.addConstr(madD[scenarioIdx] >= deviation, "d_pos_" + std::to_string(scenarioIdx));
            madDConsPos.push_back(conPos);

            // d_k >= -deviation
            GRBConstr conNeg = model.addConstr(madD[scenarioIdx] >= -deviation, "d_neg_" + std::to_string(scenarioIdx));
            madDConsNeg.push_back(conNeg);
        }
        catch (GRBException &e)
        {
            cout << "Error in MAD deviation constraints: " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        }
    }

    // Build sum of d_k
    GRBLinExpr sumD;
    for (int scenarioIdx = 0; scenarioIdx < cntScenarios; scenarioIdx++)
    {
        sumD += madD[scenarioIdx];
    }

    // Main constraint: sum_i(mu_i * x_i) - sqrt(pi/2) * (k/n) * sum_k(d_k) >= v
    double scaleFactor = k / static_cast<double>(cntScenarios);
    try
    {
        if (probCon->vsign == Inequality::gteq)
        {
            madControlCon = model.addConstr(expectedProfit - scaleFactor * sumD, GRB_GREATER_EQUAL, v, "mad_control");
        }
        else
        {
            throw("please check this again -- we only implemented the gteq version for now");
        }
    }
    catch (GRBException &e)
    {
        cout << "Error in MAD control constraint: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }

    hasMADConstr = true;
    model.update();
}

// Remove all MAD variance control constraint components from the model
void SSFormulator::removeMADVarianceControlConstr(GRBModel &model)
{
    if (!hasMADConstr)
    {
        return;
    }

    // Remove the main MAD control constraint
    try
    {
        model.remove(madControlCon);
    }
    catch (GRBException &e)
    {
        cout << "Error removing MAD control constraint: " << e.getMessage() << endl;
    }

    // Remove d_k positive deviation constraints
    for (auto &con : madDConsPos)
    {
        try
        {
            model.remove(con);
        }
        catch (GRBException &e)
        {
            cout << "Error removing MAD positive constraint: " << e.getMessage() << endl;
        }
    }

    // Remove d_k negative deviation constraints
    for (auto &con : madDConsNeg)
    {
        try
        {
            model.remove(con);
        }
        catch (GRBException &e)
        {
            cout << "Error removing MAD negative constraint: " << e.getMessage() << endl;
        }
    }
    model.update();

    // Remove d_k variables
    for (auto &var : madD)
    {
        try
        {
            model.remove(var);
        }
        catch (GRBException &e)
        {
            cout << "Error removing d variable: " << e.getMessage() << endl;
        }
    }

    model.update();

    // Clear storage
    madD.clear();
    madDConsPos.clear();
    madDConsNeg.clear();
    hasMADConstr = false;
}