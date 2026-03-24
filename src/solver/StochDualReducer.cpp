#include "StochDualReducer.hpp"
#include <set>

void findUnionDetCVaR(vector<double> &solDet, vector<double> &solCVaR, vector<double> &solLP1)
{
    for (int i = 0; i < solDet.size(); i++)
    {
        double value = (solDet[i] + solCVaR[i]) / 2;
        solLP1[i] = value;
    }
}


vector<double> StochDualReducer::solvelcvarRS(shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions)
{
    SDRFormulator formulator(spq);
    setDecisionVarOptions(formOptions.decisionVarOptions, 0.0, 1.0, 0.0, GrbVarType::Binary);
    GRBModel model = formulator.formulate(spq,formOptions);
    SolveOptions options;
    options.reduced = formOptions.reduced;
    options.reducedIds = formOptions.reducedIds;
    vector<double> x(NTuples, 0.0);
    solve(model, x, options);
    validate(model, x, spq, options);

    double E = calculateE(x);
    deb(E);
    // double ub = E/500;
    // deb(E, ub);

    // formulator.updateBound(model, ub);
    // vector<double> lcvarSol(NTuples, 0.0);
    // solve(model, lcvarSol, options);
    // validate(model, lcvarSol, spq, options);
    // deb(lcvarSol);
    return x;
}