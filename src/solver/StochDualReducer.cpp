#include "StochDualReducer.hpp"
#include <set>

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

void findUnionDetCVaR(vector<double> &solDet, vector<double> &solCVaR, vector<double> &solLP1)
{
    for (int i = 0; i < solDet.size(); i++)
    {
        double value = (solDet[i] + solCVaR[i]) / 2;
        solLP1[i] = value;
    }
}

vector<double> StochDualReducer::stage1(shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions, GRBModel &modelDetLP, GRBModel &modelCVaRLP)
{
    vector<double> solLP1;
    vector<double> solDetLP1;
    vector<double> solCVaRLP1;

    for (int i = 0; i < NTuples; i++)
    {
        solLP1.push_back(0.0);
        solDetLP1.push_back(0.0);
        solCVaRLP1.push_back(0.0);
    }
    SolveOptions options;
    options.reduced = formOptions.reduced;
    options.reducedIds = formOptions.reducedIds;
    solve(modelDetLP, solDetLP1, options);
    solve(modelCVaRLP, solCVaRLP1, options);
    findUnionDetCVaR(solDetLP1, solCVaRLP1, solLP1);
    return solLP1;
}

vector<int> StochDualReducer::stage2(std::shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions,
                                     SDRFormulator &formulatorDet,
                                     SDRFormulator &formulatorCVaR,
                                     GRBModel &modelDetLP,
                                     GRBModel &modelCVaRLP,
                                     vector<double> &solLP1)
{
    // if union empty no sol
    double E = calculateE(solLP1);
    double ub = E / formOptions.qSz;
    vector<double> solLP2;
    vector<double> solDetLP2;
    vector<double> solCVaRLP2;
    if (ub > 0)
    {
        for (int i = 0; i < NTuples; i++)
        {
            solLP2.push_back(0.0);
            solDetLP2.push_back(0.0);
            solCVaRLP2.push_back(0.0);
        }
        // find the union of both solutions and store it in solLP1 in solve2 we ubdate the upper bound
        formulatorDet.updateBound(modelDetLP, ub);
        formulatorCVaR.updateBound(modelCVaRLP, ub);
        SolveOptions options;
        options.reduced = formOptions.reduced;
        options.reducedIds = formOptions.reducedIds;
        solve(modelDetLP, solDetLP2, options);
        solve(modelCVaRLP, solCVaRLP2, options);
        findUnionDetCVaR(solDetLP2, solCVaRLP2, solLP2);
    }

    vector<int> reducedIds;
    findUnion(solLP1, solLP2, reducedIds);
    deb(reducedIds.size());
    if (reducedIds.size() < formOptions.qSz)
    {
        double low = ub;
        double high = 1.0;
        double eps = 1e-6;
        while (high - low > eps)
        {
            reducedIds.clear();
            double mid = low + (high - low) / 2;
            deb(low, high, mid);
            vector<double> solDetLP2;
            vector<double> solCVaRLP2;
            for (int i = 0; i < NTuples; i++)
            {
                solDetLP2.push_back(0.0);
                solCVaRLP2.push_back(0.0);
            }
            // find the union of both solutions and store it in solLP1 in solve2 we ubdate the upper bound
            formulatorDet.updateBound(modelDetLP, mid);
            formulatorCVaR.updateBound(modelCVaRLP, mid);
            SolveOptions options;
            options.reduced = formOptions.reduced;
            options.reducedIds = formOptions.reducedIds;
            solve(modelDetLP, solDetLP2, options);
            solve(modelCVaRLP, solCVaRLP2, options);
            findUnionDetCVaR(solDetLP2, solCVaRLP2, solLP2);
            findUnion(solLP1, solLP2, reducedIds);
            deb(reducedIds.size());

            if (reducedIds.size() >= formOptions.qSz) //if there's a solution return
            {
                break;
            }

            if (solLP2.size() > 0) //if the system was feasible decrease the upper bound
            {
                high = mid;
            }

            if (solLP2.size() == 0) //if the system was infeasible increase the upper bound
            {
                low = mid;
            }
        }
    }
    return reducedIds;
}

SolutionMetadata<int> StochDualReducer::stage3(shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions, map<string, Option> &curveFitOptions, int z, 
                                               const std::chrono::steady_clock::time_point& start_time,
                                               double timeout_seconds)
{
    // we have the reduced ids, now we can formulate the summary search
    SummarySearch SS(this->M, spq, this->epsilon);
    SSFormulator formulator(spq);
    SolutionMetadata<int> sol = SS.summarySearch<int>(this->spq, formulator, formOptions, curveFitOptions, z, start_time, timeout_seconds);
    return sol;
}

// The idea here is that if the LP1 system is infeasible
// --> formulate separate LPDet1 and LPCVaR1, so that we get tuples that will satisfy both constraints
// i.e. narrow the number of tuples down only to tuples that are meaningful for the problem itself"""
SolutionMetadata<int> StochDualReducer::stochDualReducer(shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions, map<string, Option> &curveFitOptions, 
                                               const std::chrono::steady_clock::time_point& start_time,
                                               double timeout_seconds)
{
    SDRFormulator formulatorDet(spq);
    SDRFormulator formulatorCVaR(spq);
    SDRFormulator lpformulator(spq);

    DecisionVarOptions decVarOptions;
    setDecisionVarOptions(decVarOptions, 0.0, 1.0, 0.0, GrbVarType::Continuous);
    formOptions.decisionVarOptions = decVarOptions;

    // formulate the deterministic and CVaR LP
    // GRBModel modelDetLP = formulatorDet.formulateDeterministicLP(spq, formOptions);
    // GRBModel modelCVaRLP = formulatorCVaR.formulateCVaRLP(spq, formOptions);

    // modelDetLP.write("/home/fm2288/StochasticPackageQuery/src/modelDetLP.lp");
    // modelCVaRLP.write("/home/fm2288/StochasticPackageQuery/src/modelCVaRLP.lp");
    GRBModel modelLP = lpformulator.formulateDeterministicLP(spq, formOptions);
    // -------------- STAGE 1 ----------------

    vector<double> solLP1 = stage1(spq, formOptions, modelLP, modelLP);
    int cnt = 0;
    for(int i = 0; i < solLP1.size(); i++)
    {
        if(solLP1[i] > 0)
        {
            cnt += 1;
        }
    }
    deb(cnt);

    // -------------- STAGE 2 and 3 ----------------
    //vector<int> reducedIds = stage2(spq, formOptions, formulatorDet, formulatorCVaR, modelDetLP, modelCVaRLP, solLP1);
    //formOptions.reducedIds = reducedIds;
    formOptions.reduced = true;
    SolutionMetadata<int> sol; //= stage3(spq, formOptions, curveFitOptions, formOptions.Z, start_time, timeout_seconds);
    //sol.qSz = reducedIds.size();
    return sol;
}

    // this is the original SDR approach
    // while(true)
    // {
    //     auto current_time = std::chrono::steady_clock::now();
    //     if(current_time - start_time > std::chrono::seconds(600))
    //     {
    //         cout << "Time limit reached, exiting..." << endl;
    //         break;
    //     }
    //     vector<int> reducedIds = stage2(spq, formOptions, formulatorDet, formulatorCVaR, modelDetLP, modelCVaRLP, solLP1);
    //     SolutionMetadata<int> sol = stage3(spq, formOptions, curveFitOptions, formOptions.Z, start_time, 600);
    //     if (sol.isFeasible || formOptions.qSz == NTuples)
    //     {
    //         sol.qSz = reducedIds.size();
    //         return sol;
    //     }
    //     else
    //     {
    //         formOptions.qSz = min(2 * formOptions.qSz, NTuples);
    //     }
    // }