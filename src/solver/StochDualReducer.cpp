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


SolutionMetadata<double> StochDualReducer::solveLCVaR(shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions, SolveOptions &solveOptions)
{
    SDRFormulator formulator(spq);
    GRBModel model = formulator.formulate(spq, formOptions);
    vector<double> x;
    initializeVector(x, NTuples, 0.0);
    solve(model, x, solveOptions);
    
    SolutionMetadata<double> sol;
    if(!x.empty())
    {
        validate(model, x, spq, solveOptions);
        sol.x = x;
        sol.w = this->W_q;
        sol.bestRk = this->r;
        sol.isFeasible = isFeasible(this->r);
    }
    return sol;
}