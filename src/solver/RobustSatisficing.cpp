#include "RobustSatisficing.hpp"
#include "StochDualReducer.hpp"

std::vector<int> randomTupleSelect(
    int NTuples,
    int qTarget,
    vector<int> &reducedIds)
{
    std::unordered_set<int> used;

    for (int x : reducedIds)
        used.insert(x);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(1, NTuples);

    while (reducedIds.size() < qTarget)
    {
        int x = dist(gen);
        if (used.insert(x).second)
        {
            reducedIds.push_back(x);
        }
    }

    return reducedIds;
}

void RobustSatisficing::populateMapNonZero(map<int, double> &reducedIdsMap, const vector<double> &sol)
{
    for (int i = 0; i < sol.size(); i++)
    {
        if (sol[i] > 0.0)
        {
            int id = i + 1;
            if (reducedIdsMap.find(id) == reducedIdsMap.end())
            {
                reducedIdsMap[id] = data.stockExpectedProfit[id]; // we store the id of the tuple
            }
        }
    }
}

void RobustSatisficing::populateMapFromVector(map<int, double> &reducedIdsMap, const vector<int> &reduced)
{
    for (int i = 0; i < reduced.size(); i++)
    {
        int id = reduced[i];
        if (reducedIdsMap.find(id) == reducedIdsMap.end())
        {
            reducedIdsMap[id] = data.stockExpectedProfit[id];
        }
    }
}

void RobustSatisficing::getReduced(vector<int> &reducedIds, vector<double> &solDet, vector<double> &solStage1, int qTarget)
{
    vector<pair<int, double>> solutions;
    for (int i = 0; i < solDet.size(); i++)
    {
        double max_coef = max(solStage1[i], solDet[i]);
        if (max_coef > 0.0)
        {
            solutions.push_back(make_pair(i + 1, max_coef)); // only add solutions that contain a positive val
        }
    }

    sort(solutions.begin(), solutions.end(), [](const pair<int, double> &a, const pair<int, double> &b)
         {
             return a.second > b.second; // Sort in descending order based on the second element // the more positive the better
         });

    int included = 0;
    for (int i = 0; i < solutions.size(); i++)
    {
        reducedIds.push_back(solutions[i].first);
        included++;
        if (included == qTarget)
        {
            break;
        }
    }
}

SolutionMetadata<int> RobustSatisficing::solveDeterministic(std::shared_ptr<StochasticPackageQuery> spq, SolveOptions &solveOptions)
{
    SSFormulator formulator(spq);
    DecisionVarOptions decVarOptions;
    setDecisionVarOptions(decVarOptions, 0.0, 1.0, 0.0, GrbVarType::Binary);
    FormulateOptions formOptions;
    formOptions.decisionVarOptions = decVarOptions;
    formOptions.iteration = 0;

    GRBModel model = formulator.formulate(spq, formOptions);
    vector<int> x(NTuples, 0);
    solveOptions.reduced = false;
    solve(model, x, solveOptions);
    if(x.size() > 0)
    {
        validate(model, x, spq, solveOptions);
        vector<double> feasibilities;
        vector<double> surpluses;
        SolType res;
        for (int i = 0; i < NTuples; i++)
        {
            if (x[i] > 0)
            {
                res[i + 1] = x[i];
            }
        }
        bool feas = checker->feasible(res, feasibilities, surpluses);
        double validObj = checker->getObjective(res);
    
        gpro.clock("effectiveRuntime");
    
        bool currentIsBetter = isCurrentBetterThanBest(r, this->W_q, bestSolGlobal, spq->obj->objSense);
        if (currentIsBetter)
        {
            bestSolGlobal.x = x;
            bool isFeas = isFeasible(r);
            deb(isFeasible(r));
            bestSolGlobal.isFeasible = isFeas;
            bestSolGlobal.bestRk = r;
            bestSolGlobal.w = this->W_q;
            bestSolGlobal.solutions.push_back(make_tuple(gpro.getTime("effectiveRuntime"), satisfiedScenarios, r, this->W_q, feasibilities, surpluses, validObj, 0, 0));
        }

    }
    gpro.stop("effectiveRuntime");

    return bestSolGlobal;
}

SolutionMetadata<int> RobustSatisficing::stochasticDualReducer(std::shared_ptr<StochasticPackageQuery> spq, SolveOptions &solveOptions)
{
    SolutionMetadata<int> solDeterministic = solveDeterministic(spq, solveOptions);
    if(solDeterministic.x.size() == 0)
    {
        return solDeterministic;
    }
    
    if (!solveOptions.timeBudgeted)
    {
        if (solDeterministic.isFeasible)
        {
            cout << "Deterministic Solution is Feasible" << endl;
            return solDeterministic;
        }
    }
    else
    {
        gpro.stop("effectiveRuntime");
        double usedTime = gpro.getTime("effectiveRuntime") / 1000;
        if (solDeterministic.isFeasible && usedTime >= solveOptions.timeBudget)
        {
            cout << "Deterministic Solution is Feasible" << endl;
            return solDeterministic;
        }
        gpro.clock("effectiveRuntime");
    }

    int qTarget = 500;
    RSFormulator formulator(spq);
    DecisionVarOptions decVarOptions;
    setDecisionVarOptions(decVarOptions, 0.0, 1.0, 0.0, GrbVarType::Continuous);
    FormulateOptions formOptions;
    formOptions.decisionVarOptions = decVarOptions;
    formOptions.iteration = 0;

    formulator.formulateBestObjProblem(formulator.modelBestObj, spq, formOptions);
    std::vector<double> solDet;
    initializeVector(solDet, NTuples, 0.0);

    solveOptions.reduced = formOptions.reduced;
    solveOptions.reducedIds = formOptions.reducedIds;
    solve(formulator.modelBestObj, solDet, solveOptions);
    double Z0 = formulator.modelBestObj.get(GRB_DoubleAttr_ObjVal);
    validate(formulator.modelBestObj, solDet, spq, solveOptions);
    // //--------------------- STAGE 1 -----------------------------//
    // this map contains the reducedIds and the expected profit for each id
    // formOptions.stage1 = true;
    map<int, double> reducedIdsMap;
    formOptions.innerConstraints = this->innerConstraints;
    formOptions.Z = 1;
    formOptions.Zinit = formOptions.Z;

    //vector<int> reducedIds = reduceTuplesStageNoObjCons(spq, formOptions, solveOptions, solDet, qTarget);

    vector<int> reducedIds = reduceTuplesStageNoObjConsNoObj(spq, formOptions, solveOptions, solDet, qTarget);
    cout<<"REDUCEDIDS SIZE = "<<reducedIds.size()<<endl;

    vector<int> finalReducedIds;
    populateMapFromVector(reducedIdsMap, reducedIds);
    finalReduce(reducedIdsMap, finalReducedIds, qTarget);
    deb(finalReducedIds, finalReducedIds.size(), reducedIds.size());
    formOptions.qSz = finalReducedIds.size();
    //--------------------- STAGE 2 -----------------------------//
    formOptions.stage1 = false;
    formOptions.Z = min(cntScenarios, qTarget);
    formOptions.Zinit = formOptions.Z;
    // formOptions.computeActiveness = true;
    RSFormulator formulatorSAA(spq);
    formOptions.low = 0.0;
    formOptions.high = 1.0;
    formOptions.objCons = true;
    formOptions.reducedIds = finalReducedIds;
    formOptions.reduced = true;
    solveOptions.reduced = formOptions.reduced;
    solveOptions.reducedIds = formOptions.reducedIds;
    solveOptions.activenessAbsoluteValue = true;
    setDecisionVarOptions(decVarOptions, 0.0, 1.0, 0.0, GrbVarType::Continuous);
    formOptions.decisionVarOptions = decVarOptions;
    formOptions.iteration = 0;
    formOptions.objValue = NEG_INF;

    if (cntScenarios > qTarget)
    {
        double bestEps = findBestObjectiveStage(spq, formOptions, solveOptions, Z0);
        deb(bestEps);
    }
    else
    {
        cout << "SKIPPING STAGE 2" << " " << cntScenarios << " " << qTarget << endl;
    }
    // // // ------------------- FINAL ----------------------------//
    int qInc = 500;
    int exp = 2;
    int qTargetScenarios = qTarget;
    while (true)
    {
        SSFormulator formulatorStage3(spq);
        formOptions.iteration = 0;
        formOptions.indicators = false;
        formOptions.Z = min(cntScenarios, qTargetScenarios);
        deb(qTarget, formOptions.Z);
        formOptions.reduced = true;
        formOptions.reducedIds = finalReducedIds;
        solveOptions.reduced = formOptions.reduced;
        solveOptions.reducedIds = formOptions.reducedIds;
        formOptions.activeness = this->activeness; // if stage 2 is skipped, the values are filled from stage 1
        formOptions.finalStage = true;
        // formOptions.partitionSpreadActiveness = true;
        formOptions.activenessPartition = true;
        formOptions.partitionMostActive = true;
        setDecisionVarOptions(decVarOptions, 0.0, 1.0, 0.0, GrbVarType::Binary);
        formOptions.decisionVarOptions = decVarOptions;
        SummarySearch<int> SS(M, spq, 0.20, DB_valid);
        formOptions.partitionMaximums.clear();
        SolutionMetadata<int> sol = SS.summarySearchStage3(spq, formulatorStage3, formOptions, solveOptions);
        cout << "RETURNED TO WHILE LOOP" << endl;

        
        bool isCurrentBetter = false;
        if(sol.x.size() > 0)
        {
            isCurrentBetter = isCurrentBetterThanBest(sol.bestRk, sol.w, bestSolGlobal, spq->obj->objSense);
        }
        if (isCurrentBetter)
        {
            bestSolGlobal.bestRk = sol.bestRk;
            bestSolGlobal.w = sol.w;
            bestSolGlobal.isFeasible = sol.isFeasible;
            bestSolGlobal.x = sol.x;
        }
        bestSolGlobal.solutions.insert(bestSolGlobal.solutions.end(), sol.solutions.begin(), sol.solutions.end());

        bestSolGlobal.terminate = (qTarget == NTuples) && (qTargetScenarios == cntScenarios);
        gpro.stop("effectiveRuntime");
        bool timeOut = (gpro.getTime("effectiveRuntime") / 1000 > solveOptions.timeout_seconds);
        gpro.clock("effectiveRuntime");

        deb(solveOptions.timeBudgeted, timeOut);
        if (solveOptions.timeBudgeted)
        {
            if (timeOut || bestSolGlobal.terminate)
            {
                cout << "EXITING HERE" << endl;
                bestSolGlobal.qSz = finalReducedIds.size();
                bestSolGlobal.qScenarios = formOptions.Z;
                deb(bestSolGlobal.qSz);
                return bestSolGlobal;
            }
            else
            {
                cout << "TUKAAA" << endl;
                qTarget += qInc;
                qTargetScenarios += qInc;
                qTarget = min(NTuples, qTarget);
                qTargetScenarios = min(cntScenarios, qTargetScenarios);
                qInc *= exp;
                randomTupleSelect(NTuples, qTarget, finalReducedIds);
                deb(finalReducedIds.size());
            }
        }
        else
        {
            if (bestSolGlobal.isFeasible || bestSolGlobal.terminate || timeOut)
            {
                cout << "EXITING HERE" << endl;
                bestSolGlobal.qSz = finalReducedIds.size();
                bestSolGlobal.qScenarios = formOptions.Z;
                deb(bestSolGlobal.qSz);
                return bestSolGlobal;
            }
            else
            {
                cout << "TUKAAA" << endl;
                qTarget += qInc;
                qTargetScenarios += qInc;
                qTarget = min(NTuples, qTarget);
                qTargetScenarios = min(cntScenarios, qTargetScenarios);
                qInc *= exp;
                randomTupleSelect(NTuples, qTarget, finalReducedIds);
                deb(finalReducedIds.size());
            }
        }
    }
    // SolutionMetadata<int> sol;
    // sol.qSz = finalReducedIds.size();
    // deb(sol.qSz);
    return bestSolGlobal;
}

std::vector<int> RobustSatisficing::reduceTuplesStageNoObjCons(std::shared_ptr<StochasticPackageQuery> spq,
                                                               FormulateOptions &formOptions,
                                                               SolveOptions &solveOptions, vector<double> &solDet, int qTarget)
{
    vector<int> reducedIds;
    double epsilonStage1 = 1e7;
    formOptions.objCons = false;
    formOptions.iteration = 0;
    RSFormulator formulator(spq);
    SummarySearch<double> SS(M, spq, epsilonStage1, DB_valid);
    SolutionMetadata<double> sol = SS.summarySearchRS(SS.spq, formulator, formOptions, solveOptions);
    findNonzero(reducedIds, sol.x);
    vector<int> unioned;
    getReduced(unioned, solDet, sol.x, qTarget);
    deb(reducedIds.size(), unioned.size());
    return unioned;
}

vector<int> RobustSatisficing::reduceTuplesStageNoObjConsNoObj(std::shared_ptr<StochasticPackageQuery> spq,
                                                               FormulateOptions &formOptions,
                                                               SolveOptions &solveOptions, vector<double> &solDet, int qTarget)
{
    formOptions.includeObjectiveFunction = false;
    solveOptions.includeObjectiveFunction = false;
    double low = 0.0;
    double high = 1.0; 
    double epsBinSearch = 1e-3;

    formOptions.partitionMaximums.clear();
    SolutionMetadata<double> best;
    while(high - low > epsBinSearch)
    {
        double mid = low + (high - low) / 2;
        double epsilonStage1 = 1e7;
        formOptions.objCons = false;
        formOptions.iteration = 0;
        DecisionVarOptions decVarOptions;
        setDecisionVarOptions(decVarOptions, 0.0, mid, 0.0, GrbVarType::Continuous);
        deb(low,mid, high);
        formOptions.decisionVarOptions = decVarOptions;
        formOptions.iteration = 0;
        RSFormulator formulator(spq);
        SummarySearch<double> SS(M, spq, epsilonStage1, DB_valid);
        SolutionMetadata<double> sol = SS.summarySearchRS(SS.spq, formulator, formOptions, solveOptions);

        if(sol.x.empty())
        {
            low = mid;
        }else
        if(sol.isFeasible)
        {
            high = mid;
            int sense = GRB_MAXIMIZE;
            bool currentIsBetter = isCurrentBetterThanBest(sol.bestRk, sol.w, best, sense);
            if(currentIsBetter)
            {
                best.bestRk = sol.bestRk;
                best.w = sol.w;
                best.isFeasible = sol.isFeasible;
                best.x = sol.x;
            }
        }else
        {
            low = mid;
            int sense = GRB_MAXIMIZE;
            bool currentIsBetter = isCurrentBetterThanBest(sol.bestRk, sol.w, best, sense);
            if(currentIsBetter)
            {
                best.bestRk = sol.bestRk;
                best.w = sol.w;
                best.isFeasible = sol.isFeasible;
                best.x = sol.x;
            }
        }
    }
    cout<<"I FINISHED BINARY SEARCH"<<endl;
    // deb(best.x);
    vector<int> reducedIds;
    findNonzero(reducedIds, best.x);
    vector<int> unioned;
    getReduced(unioned, solDet, best.x, qTarget);
    deb(reducedIds.size(), unioned.size());
    formOptions.includeObjectiveFunction = true;
    solveOptions.includeObjectiveFunction = true;
    return unioned;
}

double RobustSatisficing::findBestObjectiveStage(std::shared_ptr<StochasticPackageQuery> spq,
                                                 FormulateOptions &formOptions,
                                                 SolveOptions &solveOptions, double Z0)
{
    cout << "Stage 2" << endl;
    int totalSystems = 0;
    int totalSolveTime = 0;

    double low = formOptions.low;
    double high = formOptions.high;
    double eps = 1e-5;
    vector<double> bestRk;
    double best = -1; // a dummy for epsilon,
    while (high - low > eps)
    {
        double mid = low + (high - low) / 2; // <- this is our new epsilon
        deb(low, mid, high);
        double Z = (1 - mid) * Z0; // <- this is the new Objective value
        double epsilonStage1 = 1e7;
        formOptions.objValue = Z; // <- this is the new objective value
        formOptions.iteration = 0;
        RSFormulator formulator(spq);
        SummarySearch<double> SS(M, spq, epsilonStage1, DB_valid);
        formOptions.partitionMaximums.clear();
        SolutionMetadata<double> sol = SS.summarySearchRS(SS.spq, formulator, formOptions, solveOptions);

        if (sol.isFeasible) // if feasible this is our current best
        {
            best = mid;
            high = mid;
            this->activeness = sol.bestActiveness;
            bestRk = sol.bestRk;
        }
        else
        {
            if (sol.x.size() > 0) // it's worth comparing only if there is a solution
            {
                double worstRkSol = 1e7;
                double worstRkBestSol = 1e7;
                if (bestRk.empty()) // if there's no best so far, set a best
                {
                    this->activeness = sol.bestActiveness;
                    bestRk = sol.bestRk;
                    best = mid;
                }
                else
                {
                    // compare the worst surpluses
                    for (int i = 0; i < r.size(); i++)
                    {
                        if (sol.bestRk[i] < worstRkSol)
                        {
                            worstRkSol = r[i];
                        }
                        if (bestRk[i] < worstRkBestSol)
                        {
                            worstRkBestSol = bestRk[i];
                        }
                    }
                    // get only if the current worst surplus is better than the best so far
                    if (worstRkSol > worstRkBestSol)
                    {
                        this->activeness = sol.bestActiveness;
                        bestRk = sol.bestRk;
                        best = mid;
                    }
                }
            }
            low = mid;
        }
    }
    if (best == -1)
    {
        best = low + (high - low) / 2;
    }
    double interval = high - low;
    deb(totalSystems, totalSolveTime);
    return best;
}

vector<int> RobustSatisficing::finalReduce(map<int, double> &reducedIdsMap, vector<int> &reducedIds, int q)
{
    if (reducedIdsMap.size() < q)
    {
        for (int i = 0; i < data.stockExpectedProfitSorted.size(); i++)
        {
            int id = data.stockExpectedProfitSorted[i].first + 1;
            if (reducedIdsMap.find(id) == reducedIdsMap.end())
            {
                reducedIdsMap[id] = data.stockExpectedProfitSorted[i].second;
                if (reducedIdsMap.size() >= q)
                {
                    break;
                }
            }
        }
    }
    else if (reducedIdsMap.size() > q)
    {
        vector<pair<int, double>> tmp(reducedIdsMap.begin(), reducedIdsMap.end());

        // Sort ascending by value
        sort(tmp.begin(), tmp.end(), [](const auto &a, const auto &b)
             { return a.second < b.second; });

        // Remove the smallest (current_size - qRight) elements
        int removeCount = reducedIdsMap.size() - q;
        for (int i = 0; i < removeCount; i++)
        {
            reducedIdsMap.erase(tmp[i].first);
        }
    }
    for (auto &p : reducedIdsMap)
    {
        reducedIds.push_back(p.first);
    }
    return reducedIds;
}