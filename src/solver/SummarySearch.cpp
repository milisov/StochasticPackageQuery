#include "SummarySearch.hpp"
#include <boost/algorithm/string/join.hpp>
#include <fmt/ranges.h>


void SummarySearch::guessOptimalConservativeness(std::vector<std::vector<pair<double, double>>> &history, std::vector<double> &alpha)
{
    fitter.fit_and_predict(history, alpha);
}

double calculateE(vector<double> &x)
{
    double E = 0.0;
    for (int i = 0; i < x.size(); i++)
    {
        E += x[i];
    }

    return E;
}

void findUnion(vector<double> &solLP1, vector<double> &solLP2, vector<int> &reducedIds)
{
    for (int i = 0; i < solLP1.size(); i++)
    {
        if (solLP1[i] > 0 || solLP2[i] > 0)
        {
            reducedIds.push_back(i + 1);
        }
    }
}

void findUnionDetCVaR(vector<double> &solDet, vector<double> &solCVaR, set<int> &reducedIds)
{
    for (int i = 0; i < solDet.size(); i++)
    {
        if (solCVaR[i] > 0 || solDet[i] > 0)
        {
            reducedIds.insert(i + 1);
        }
    }
}

void findUnionDetCVaR2(vector<double> &solDet, vector<double> &solCVaR, vector<double> &solLP1)
{
    for (int i = 0; i < solDet.size(); i++)
    {
        double value = (solDet[i] + solCVaR[i]) / 2;
        solLP1[i] = value;
    }
}

void populateReducedIds(vector<double> &solLP, set<int> &reducedIds)
{
    for (int i = 0; i < solLP.size(); i++)
    {
        reducedIds.insert(i + 1);
    }
}

// The idea here is that if the LP1 system is infeasible
// --> formulate separate LPDet1 and LPCVaR1, so that we get tuples that will satisfy both constraints
// i.e. narrow the number of tuples down only to tuples that are meaningful for the problem itself"""
// SolutionMetadata<int> SummarySearch::stochDualReducer(shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions, map<string, Option> &curveFitOptions)
// {
//     cout << "StochDualReducer" << endl;
//     SSFormulator formulatorDet(spq);
//     SSFormulator formulatorCVaR(spq);

//     DecisionVarOptions decVarOptionsDet;
//     setDecisionVarOptions(decVarOptionsDet, 0.0, 1.0, 0.0, GrbVarType::Continuous);
//     FormulateOptions formOptionsDet = formOptions;
//     formOptionsDet.decisionVarOptions = decVarOptionsDet;

//     DecisionVarOptions decVarOptionsCVaR;
//     setDecisionVarOptions(decVarOptionsCVaR, 0.0, 1.0, 0.0, GrbVarType::Continuous);
//     FormulateOptions formOptionsCVaR = formOptions;
//     formOptionsCVaR.decisionVarOptions = decVarOptionsCVaR;

//     vector<double> solLP1;
//     vector<double> solDetLP1;
//     vector<double> solCVaRLP1;

//     for (int i = 0; i < NTuples; i++)
//     {
//         solLP1.push_back(0.0);
//         solDetLP1.push_back(0.0);
//         solCVaRLP1.push_back(0.0);
//     }
//     // formulate and solve the deterministic and CVaR LP
//     GRBModel modelDetLP1 = formulatorDet.formulateDeterministicLP(spq, formOptionsDet);
//     GRBModel modelCVaRLP1 = formulatorCVaR.formulateCVaRLP(spq, formOptionsCVaR);
//     SolveOptions options;
//     options.reduced = formOptions.reduced;
//     options.reducedIds = formOptions.reducedIds;
//     solve(modelDetLP1, solDetLP1, options);
//     solve(modelCVaRLP1, solCVaRLP1, options);
//     findUnionDetCVaR2(solDetLP1, solCVaRLP1, solLP1);
//     // if union empty no sol
//     double E = calculateE(solLP1);
//     while (true)
//     {
//         double ub = E / formOptions.qSz;
//         vector<double> solLP2;
//         vector<double> solDetLP2;
//         vector<double> solCVaRLP2;
//         if (ub > 0)
//         {
//             for (int i = 0; i < NTuples; i++)
//             {
//                 solLP2.push_back(0.0);
//                 solDetLP2.push_back(0.0);
//                 solCVaRLP2.push_back(0.0);
//             }
//             // find the union of both solutions and store it in solLP1 in solve2 we ubdate the upper bound
//             formulatorDet.updateBound(modelDetLP1, ub);
//             formulatorCVaR.updateBound(modelCVaRLP1, ub);
//             SolveOptions options;
//             options.reduced = formOptions.reduced;
//             options.reducedIds = formOptions.reducedIds;
//             solve(modelDetLP1, solDetLP2, options);
//             solve(modelCVaRLP1, solCVaRLP2, options);
//             findUnionDetCVaR2(solDetLP2, solCVaRLP2, solLP2);
//         }
//         else
//         {
//             // if the solution is empty even on deterministic LP then it's infeasible
//             vector<int> v;
//             SolutionMetadata<int> sol(v, -1.0, -1.0, 0);
//             return sol;
//         }

//         vector<int> reducedIds;
//         findUnion(solLP1, solLP2, reducedIds);

//         int steps = 0;
//         if (reducedIds.size() < formOptions.qSz)
//         {
//             double low = ub;
//             double high = 1.0;
//             double eps = 1e-6;
//             while (high - low > eps)
//             {
//                 steps += 1;
//                 reducedIds.clear();
//                 double mid = low + (high - low) / 2;
//                 cout << "Low = " << low << endl;
//                 cout << "High = " << high << endl;
//                 cout << "Mid = " << mid << endl;
//                 vector<double> solDetLP2;
//                 vector<double> solCVaRLP2;
//                 for (int i = 0; i < NTuples; i++)
//                 {
//                     solDetLP2.push_back(0.0);
//                     solCVaRLP2.push_back(0.0);
//                 }
//                 // find the union of both solutions and store it in solLP1 in solve2 we ubdate the upper bound
//                 formulatorDet.updateBound(modelDetLP1, mid);
//                 formulatorCVaR.updateBound(modelCVaRLP1, mid);
//                 SolveOptions options;
//                 options.reduced = formOptions.reduced;
//                 options.reducedIds = formOptions.reducedIds;
//                 solve(modelDetLP1, solDetLP2, options);
//                 solve(modelCVaRLP1, solCVaRLP2, options);
//                 findUnionDetCVaR2(solDetLP2, solCVaRLP2, solLP2);

//                 findUnion(solLP1, solLP2, reducedIds);

//                 if (reducedIds.size() >= formOptions.qSz)
//                 {
//                     break;
//                 }

//                 if (reducedIds.size() < formOptions.qSz && solLP2.size() > 0)
//                 {
//                     high = mid;
//                 }

//                 if (solLP2.size() == 0)
//                 {
//                     low = mid;
//                 }
//             }
//         }

//         SSFormulator formulator(spq);
//         formOptions.reducedIds = reducedIds;
//         formOptions.reduced = true;
//         SolutionMetadata<int> sol = summarySearch<int>(this->spq, formulator, formOptions, curveFitOptions, 1);
//         if (sol.isFeasible || formOptions.qSz == NTuples)
//         {
//             sol.binarySearchSteps = steps;
//             sol.qSz = formOptions.qSz;
//             return sol;
//         }
//         else
//         {
//             formOptions.qSz = min(2 * formOptions.qSz, NTuples);
//         }
//     }
// }