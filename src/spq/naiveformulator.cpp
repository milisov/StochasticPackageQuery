#include "naiveformulator.hpp"
#include <fmt/ranges.h>
#include <boost/algorithm/string/join.hpp>
#include <cmath>
#include <set>

using namespace std;


void NaiveFormulator::randomScenarioSelection(int Z, int conOrder, set<int> &selectedScenarios, unsigned int seed)
{
    Z = std::min(Z, cntScenarios);
    if (Z <= 0) return;
    std::mt19937 gen(seed + static_cast<unsigned int>(conOrder));
    std::uniform_int_distribution<int> dist(0, cntScenarios - 1);

    while ((int)selectedScenarios.size() < Z)
    {
        selectedScenarios.insert(dist(gen));
    }
}

void NaiveFormulator::formProbCons(GRBModel &model, std::shared_ptr<Constraint> cons, GRBVar *xx, FormulateOptions &options, int conOrder)
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
    double p = options.p;

    set<int> selectedScenarios;
    randomScenarioSelection(options.Z, conOrder, selectedScenarios, options.randomSeed);
    int Z = selectedScenarios.size();
    double pM = ceil(p * Z);

    auto &scenarios = data.stochAttrs[attrCon->attr];
    int n = options.reduced ? options.reducedIds.size() : NTuples;

    std::vector<GRBVar> y_vars;
    for (int scenarioId : selectedScenarios)
    {
        y_vars.push_back(model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y[" + to_string(scenarioId) + "]"));
    }
    model.update();

    GRBLinExpr sumYz;
    int idx = 0;
    for (int scenarioId : selectedScenarios)
    {
        GRBLinExpr innerCon;
        if (options.reduced)
        {
            for (int i = 0; i < n; i++)
            {
                int id = options.reducedIds[i] - 1;
                innerCon += xx[i] * scenarios[id][scenarioId];
            }
        }
        else
        {
            for (int i = 0; i < NTuples; i++)
            {
                innerCon += xx[i] * scenarios[i][scenarioId];
            }
        }
        try
        {
            if (probCon->vsign == Inequality::gteq)
            {
                model.addGenConstrIndicator(y_vars[idx], 1, innerCon, GRB_GREATER_EQUAL, v, "prob_indicator[" + to_string(scenarioId) + "]");
            }
            else
            {
                model.addGenConstrIndicator(y_vars[idx], 1, innerCon, GRB_LESS_EQUAL, v, "prob_indicator[" + to_string(scenarioId) + "]");
            }
        }
        catch (GRBException &e)
        {
            deb("Error code 8 = ", e.getErrorCode());
            cout << e.getMessage() << endl;
        }
        sumYz += y_vars[idx];
        idx++;
    }
    try
    {
        if (probCon->psign == Inequality::gteq)
        {
            model.addConstr(sumYz, GRB_GREATER_EQUAL, pM, "prob_control");
        }
        else
        {
            model.addConstr(sumYz, GRB_LESS_EQUAL, pM, "prob_control");
        }
    }
    catch (GRBException &e)
    {
        cout << "Error code 9 = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
}

GRBModel NaiveFormulator::formulate(shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions)
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
    int conOrder = 0;

    for (int i = 0; i < numCons; i++)
    {
        formCountCons(model, spq->cons[i], xx.get(), formOptions);
        formSumCons(model, spq->cons[i], xx.get(), formOptions);
        if(formOptions.varianceControl)
        {
            formMADVarianceControlConstr(model, spq->cons[i], xx.get(), formOptions, conOrder);
        }else
        {
            formProbCons(model, spq->cons[i], xx.get(), formOptions, conOrder);
        }
    }
    formSumObj(model, spq->obj, xx.get(), formOptions);
    formExpSumObj(model, spq->obj, xx.get(), formOptions);
    formCntObj(model, spq->obj, xx.get(), formOptions);
    deb("formulated");
    model.update();
    return model;
}

// Formulates MAD-based variance control constraint with reduced scenarios:
// v <= x^T * mu_all - k * MAD_reduced
// where mu_all is the mean from ALL scenarios, and MAD_reduced is computed
// only over the Z most active scenarios.
// MAD_reduced = (1/Z) * sum_{k in active} |sum_i x_i * (p_i^(k) - mu_all_i)|
void NaiveFormulator::formMADVarianceControlConstr(GRBModel &model,
                                                   std::shared_ptr<Constraint> cons,
                                                   GRBVar *xx,
                                                   FormulateOptions &options,
                                                   int &conOrder)
{
    if (!options.varianceControl)
    {
        return;
    }

    shared_ptr<ProbConstraint> probCon;
    shared_ptr<AttrConstraint> attrCon;

    bool isstoch = isStochastic(cons, probCon, attrCon);
    if (!isstoch)
    {
        return;
    }

    double v = spq->getValue(probCon->v);
    double p = spq->getValue(probCon->p);
    double k = std::sqrt(p / (1.0 - p)) * options.kMADValues[conOrder];

    // Mean from ALL scenarios (not just the reduced ones)
    const std::vector<double> &mu = data.stockExpectedProfit;
    const std::string &attrName = attrCon->attr;
    const std::vector<std::vector<double>> &scenarios = data.stochAttrs[attrName];

    int n = options.reduced ? options.reducedIds.size() : NTuples;
    set<int>activeScenarios;
    //activeScenarios = getReducedScenariosFromActiveness(options.activenessNoAbsValue[conOrder],options.Z, p);
    randomScenarioSelection(options.Z, conOrder, activeScenarios, options.randomSeed);
    int Z = activeScenarios.size();

    // Create auxiliary variables d_k for each active scenario
    std::vector<GRBVar> d_vars;
    for (int scenarioId : activeScenarios)
    {
        GRBVar dk = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS,
                                 "d_" + std::to_string(scenarioId));
        d_vars.push_back(dk);
    }
    model.update();

    // Build expected profit: sum_i(mu_i * x_i)
    GRBLinExpr expectedProfit;
    if (options.reduced)
    {
        for (int i = 0; i < n; i++)
        {
            int id = options.reducedIds[i] - 1;
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

    // For each active scenario, add: d_k >= |deviation_k|
    int dIdx = 0;
    for (int scenarioId : activeScenarios)
    {
        GRBLinExpr deviation;

        if (options.reduced)
        {
            for (int i = 0; i < n; i++)
            {
                int id = options.reducedIds[i] - 1;
                double diff = scenarios[id][scenarioId] - mu[id];
                deviation += diff * xx[i];
            }
        }
        else
        {
            for (int i = 0; i < n; i++)
            {
                double diff = scenarios[i][scenarioId] - mu[i];
                deviation += diff * xx[i];
            }
        }

        // d_k >= deviation and d_k >= -deviation (linearization of |deviation|)
        model.addConstr(d_vars[dIdx] >= deviation, "d_pos_" + std::to_string(scenarioId));
        model.addConstr(d_vars[dIdx] >= -deviation, "d_neg_" + std::to_string(scenarioId));
        dIdx++;
    }

    // Sum of d_k
    GRBLinExpr sumD;
    for (int i = 0; i < Z; i++)
    {
        sumD += d_vars[i];
    }

    // Main constraint: x^T*mu - (k/Z)*sum(d_k) >= v
    double scaleFactor = k / static_cast<double>(Z);
    deb(scaleFactor);
    if (probCon->vsign == Inequality::gteq)
    {
        model.addConstr(expectedProfit - scaleFactor * sumD, GRB_GREATER_EQUAL, v, "mad_control");
    }
    else
    {
        model.addConstr(expectedProfit + scaleFactor * sumD, GRB_LESS_EQUAL, v, "mad_control");
    }
    model.update();
    conOrder++;
}