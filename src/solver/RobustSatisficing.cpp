#include "RobustSatisficing.hpp"
#include "StochDualReducer.hpp"

std::vector<int> randomTupleSelect(
    int NTuples,
    int qTarget,
    vector<int> &reducedIds,
    unsigned int seed = 42)
{
    std::unordered_set<int> used;

    for (int x : reducedIds)
        used.insert(x);

    std::mt19937 gen(seed + reducedIds.size()); // offset by current size so repeated calls diverge
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

void RobustSatisficing::applyWarmstart(GRBModel &model, const std::vector<int> &vbasis, const std::vector<int> &cbasis)
{
    if (vbasis.empty() || cbasis.empty())
    {
        return;
    }

    try
    {
        GRBVar *vars = model.getVars();
        int numVars = model.get(GRB_IntAttr_NumVars);
        int vbasisSize = std::min(static_cast<int>(vbasis.size()), numVars);
        for (int i = 0; i < vbasisSize; ++i)
        {
            vars[i].set(GRB_IntAttr_VBasis, vbasis[i]);
        }
        delete[] vars;

        GRBConstr *constrs = model.getConstrs();
        int numConstrs = model.get(GRB_IntAttr_NumConstrs);
        int cbasisSize = std::min(static_cast<int>(cbasis.size()), numConstrs);
        for (int i = 0; i < cbasisSize; ++i)
        {
            constrs[i].set(GRB_IntAttr_CBasis, cbasis[i]);
        }
        delete[] constrs;
    }
    catch (GRBException &e)
    {
        std::cerr << "Warmstart apply failed: " << e.getMessage() << std::endl;
    }
}

void RobustSatisficing::saveWarmstart(GRBModel &model, std::vector<int> &vbasis, std::vector<int> &cbasis)
{
    try
    {
        GRBVar *vars = model.getVars();
        int numVars = model.get(GRB_IntAttr_NumVars);
        vbasis.resize(numVars);
        for (int i = 0; i < numVars; ++i)
        {
            vbasis[i] = vars[i].get(GRB_IntAttr_VBasis);
        }
        delete[] vars;

        GRBConstr *constrs = model.getConstrs();
        int numConstrs = model.get(GRB_IntAttr_NumConstrs);
        cbasis.resize(numConstrs);
        for (int i = 0; i < numConstrs; ++i)
        {
            cbasis[i] = constrs[i].get(GRB_IntAttr_CBasis);
        }
        delete[] constrs;
    }
    catch (GRBException &e)
    {
        std::cerr << "Warmstart save failed: " << e.getMessage() << std::endl;
    }
}

void RobustSatisficing::getReduced(vector<int> &reducedIds, vector<double> &solDet, vector<double> &solStage1, int qTarget) //this function returns at most qTarget reduced ids
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

void RobustSatisficing::getReducedFromRS(vector<int> &reducedIds, vector<double> &solStage1, int qTarget)
{
    vector<pair<int, double>> solutions;
    for (int i = 0; i < solStage1.size(); i++)
    {
        if (solStage1[i] > 0.0)
        {
            solutions.push_back(make_pair(i + 1, solStage1[i])); // only add solutions that contain a positive val
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
    spq->setTableName(DB_valid+"_validate");
    SSFormulator formulator(spq);
    DecisionVarOptions decVarOptions;
    setDecisionVarOptions(decVarOptions, 0.0, 1.0, 0.0, GrbVarType::Binary);
    FormulateOptions formOptions;
    formOptions.decisionVarOptions = decVarOptions;
    formOptions.iteration = 0;

    GRBModel model = formulator.formulate(spq, formOptions);
    SolutionMetadata<int>sol;
    vector<int> x(NTuples, 0);
    solveOptions.reduced = false;
    solve(model, x, solveOptions);
    if (x.size() > 0)
    {
        validate(model, x, spq, solveOptions);
        sol.x = x;
        sol.bestRk = this->r;
        sol.w = this->W_q;

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
        bool deterFeasible, probFeasible;
        bool feas = checker->feasible(res, feasibilities, surpluses, deterFeasible, probFeasible);
        double validObj = checker->getObjective(res);
        deb(feas, deterFeasible, probFeasible, validObj, sol.bestRk, sol.w, surpluses, feasibilities);
    }
    return sol;
}

SolutionMetadata<int> RobustSatisficing::stochasticDualReducer(std::shared_ptr<StochasticPackageQuery> spq, SolveOptions &solveOptions, map<string, bool> ablate)
{
    int qTarget = solveOptions.hyperParams["reducedTuples"];
    RSFormulator formulator(spq);
    DecisionVarOptions decVarOptions;
    setDecisionVarOptions(decVarOptions, 0.0, 1.0, 0.0, GrbVarType::Continuous);
    FormulateOptions formOptions;
    formOptions.ablate = ablate;
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
    formOptions.innerConstraints = this->innerConstraints;
    formOptions.Z = 1;
    formOptions.Zinit = formOptions.Z;
    if(formOptions.ablate["ablate"] && formOptions.ablate["stage1lcvar"])
    {
        reduceTuplesLCVaR(spq, formOptions, solveOptions, solDet, qTarget);
    }else
    if(formOptions.ablate["ablate"] && formOptions.ablate["stage1random"])
    {
        vector<int>finalReducedIds;
        randomTupleSelect(NTuples, qTarget, finalReducedIds, solveOptions.randomSeed);
        formOptions.qSz = finalReducedIds.size();
        formOptions.reducedIds = finalReducedIds;
        formOptions.reduced = true;
    }else
    {
        reduceTuplesAndScenarios(spq, formOptions, solveOptions, solDet, qTarget);
    }
    
    if(formOptions.ablate["stage2"])
    {
        return runNaiveFinalStage(spq, formOptions, solveOptions, qTarget);
    }
    return runNaiveVarianceControl(spq, formOptions, solveOptions, qTarget);
}

SolutionMetadata<int> RobustSatisficing::runNaiveVarianceControl(std::shared_ptr<StochasticPackageQuery> spq,
                                                                  FormulateOptions &formOptions,
                                                                  SolveOptions &solveOptions,
                                                                  int qTarget)
{
    DecisionVarOptions decVarOptions;
    int qTargetScenarios = solveOptions.hyperParams["reducedScenarios"];
    int exp = 2;

    // Enable variance control
    setDecisionVarOptions(decVarOptions, 0.0, 1.0, 0.0, GrbVarType::Binary);
    formOptions.decisionVarOptions = decVarOptions;
    solveOptions.enableRelaxation = true;
    formOptions.varianceControl = true;
    formOptions.includeObjectiveFunction = true;
    solveOptions.includeObjectiveFunction = true;

    while (true)
    {
        Naive naiveSolver(spq, this->DB_valid);
        formOptions.Z = qTargetScenarios;
        formOptions.randomSeed = solveOptions.randomSeed;
        solveOptions.reduced = formOptions.reduced;
        solveOptions.reducedIds = formOptions.reducedIds;
        solveOptions.foundSolutionUsingGurobiRelax = false;
        // Formulate and solve
        SolutionMetadata<int> sol = naiveSolver.formulateAndSolve<int>(spq, formOptions, solveOptions);

        gpro.stop("effectiveRuntime");
        double currentTime = gpro.getTime("effectiveRuntime");
        gpro.clock("effectiveRuntime");

        cout << "Naive iteration: qTuples=" << qTarget << ", time=" << currentTime / 1000 << "s" << endl;

        if (sol.x.size() > 0)
        {
            bestSolGlobal.solutions.insert(bestSolGlobal.solutions.end(), sol.solutions.begin(), sol.solutions.end());
            bool isCurrentBetter = isCurrentBetterThanBest(sol.bestRk, sol.w, bestSolGlobal, spq->obj->objSense);
            if(isCurrentBetter)
            {
                // Last solution is always the best
                bestSolGlobal.x = sol.x;
                bestSolGlobal.w = sol.w;
                bestSolGlobal.bestRk = sol.bestRk;
                bestSolGlobal.isFeasible = sol.isFeasible;
                bestSolGlobal.qSz = formOptions.reducedIds.size();
                bestSolGlobal.qScenarios = formOptions.Z;
            }
        }
        if(solveOptions.hyperParams["benchmark"] == 2 || formOptions.ablate["stage1"])
        {
            bestSolGlobal.terminate = (qTarget >= solveOptions.hyperParams["reducedTuples"]);
        }else
        {
            bestSolGlobal.terminate = (qTarget >= NTuples);
        }
        gpro.stop("effectiveRuntime");
        bool timeOut = (gpro.getTime("effectiveRuntime") / 1000 > solveOptions.timeout_seconds);
        
        if (timeOut || bestSolGlobal.terminate)
        {
            cout << "EXITING NAIVE VARIANCE CONTROL" << endl;
            return bestSolGlobal;
        }
        gpro.clock("effectiveRuntime");
        
        qTarget *= exp;
        qTarget = min(NTuples, qTarget);
        randomTupleSelect(NTuples, qTarget, formOptions.reducedIds, solveOptions.randomSeed);
    }
    return bestSolGlobal;
}

void RobustSatisficing::reduceTuplesAndScenarios(std::shared_ptr<StochasticPackageQuery> spq,
                                                  FormulateOptions &formOptions,
                                                  SolveOptions &solveOptions, vector<double> &solDet, int qTarget)
{
    formOptions.includeObjectiveFunction = false;
    formOptions.lcvarFormulation = true;
    solveOptions.includeObjectiveFunction = false;
    solveOptions.enableRelaxation = true;
    double low = 0.0;
    double high = 1.0;
    double epsBinSearch = 1e-5;

    formOptions.partitionMaximums.clear();
    formOptions.clearLCVaRCache = true;
    SolutionMetadata<double> best;

    // Warmstart basis vectors
    std::vector<int> vbasis;
    std::vector<int> cbasis;

    while (high - low > epsBinSearch)
    {
        double mid = low + (high - low) / 2;
        deb(low, mid, high);

        // Fresh model each iteration with UB = mid.
        // LCVaR coefficients are cached statically in RSFormulator, so the O(N×M log M)
        formOptions.objCons = false;
        formOptions.iteration = 0;
        DecisionVarOptions decVarOptions;
        setDecisionVarOptions(decVarOptions, 0.0, mid, 0.0, GrbVarType::Continuous);
        formOptions.decisionVarOptions = decVarOptions;

        RSFormulator formulator(spq);
        GRBModel model = formulator.formulate(spq, formOptions);

        applyWarmstart(model, vbasis, cbasis);

        vector<double> x;
        initializeVector(x, NTuples, 0.0);
        solve(model, x, solveOptions);

        SolutionMetadata<double> sol;
        if (!x.empty())
        {
            if (solveOptions.foundSolutionUsingGurobiRelax)
            {
                vbasis.clear();
                cbasis.clear();
                solveOptions.foundSolutionUsingGurobiRelax = false;
            }
            else
            {
                saveWarmstart(model, vbasis, cbasis);
            }

            validate(model, x, spq, solveOptions);
            sol.x = x;
            sol.w = this->W_q;
            sol.bestRk = this->r;
            sol.isFeasible = isFeasible(this->r);
            sol.bestActiveness = this->activeness;
            sol.bestActivenessNoAbsValue = this->activenessNoAbsValue;

            if (sol.isFeasible)
            {
                high = mid;
                int sense = GRB_MAXIMIZE;
                bool currentIsBetter = isCurrentBetterThanBest(sol.bestRk, sol.w, best, sense);
                if (currentIsBetter)
                {
                    best.bestRk = sol.bestRk;
                    best.w = sol.w;
                    best.isFeasible = sol.isFeasible;
                    best.bestActivenessNoAbsValue = sol.bestActivenessNoAbsValue;
                    best.bestActiveness = sol.bestActiveness;
                    best.x = sol.x;
                }
            }
            else
            {
                low = mid;
                int sense = GRB_MAXIMIZE;
                bool currentIsBetter = isCurrentBetterThanBest(sol.bestRk, sol.w, best, sense);
                if (currentIsBetter)
                {
                    best.bestRk = sol.bestRk;
                    best.w = sol.w;
                    best.isFeasible = sol.isFeasible;
                    best.x = sol.x;
                    best.bestActivenessNoAbsValue = sol.bestActivenessNoAbsValue;
                    best.bestActiveness = sol.bestActiveness;
                }
            }
        }
        else
        {
            low = mid;
        }
    }

    vector<int> finalReducedIds;
    getReducedFromRS(finalReducedIds, best.x, qTarget);
    if(finalReducedIds.size() < qTarget)
    {
        randomTupleSelect(NTuples, qTarget, finalReducedIds, solveOptions.randomSeed);
    }
    formOptions.qSz = finalReducedIds.size();
    formOptions.reducedIds = finalReducedIds;
    formOptions.reduced = true;
    formOptions.activenessNoAbsValue = best.bestActivenessNoAbsValue;
    formOptions.activeness = best.bestActiveness;
    gpro.stop("fetchingToEndOfFirstStage");
}


void RobustSatisficing::reduceTuplesLCVaR(std::shared_ptr<StochasticPackageQuery> spq,
                                           FormulateOptions &formOptions,
                                           SolveOptions &solveOptions, vector<double> &solDet, int qTarget)
{
    formOptions.includeObjectiveFunction = false;
    formOptions.lcvarFormulation = true;
    solveOptions.includeObjectiveFunction = false;
    solveOptions.enableRelaxation = true;
    double low = 0.0;
    double high = 1.0;
    double epsBinSearch = 1e-5;

    formOptions.partitionMaximums.clear();
    formOptions.clearLCVaRCache = true;
    SolutionMetadata<double> best;

    // Warmstart basis vectors
    std::vector<int> vbasis;
    std::vector<int> cbasis;

    while (high - low > epsBinSearch)
    {
        double mid = low + (high - low) / 2;
        deb(low, mid, high);

        formOptions.objCons = false;
        formOptions.iteration = 0;
        DecisionVarOptions decVarOptions;
        setDecisionVarOptions(decVarOptions, 0.0, mid, 0.0, GrbVarType::Continuous);
        formOptions.decisionVarOptions = decVarOptions;

        SDRFormulator formulator(spq);
        GRBModel model = formulator.formulate(spq, formOptions);

        applyWarmstart(model, vbasis, cbasis);

        vector<double> x;
        initializeVector(x, NTuples, 0.0);
        solve(model, x, solveOptions);

        SolutionMetadata<double> sol;
        if (!x.empty())
        {
            if (solveOptions.foundSolutionUsingGurobiRelax)
            {
                vbasis.clear();
                cbasis.clear();
                solveOptions.foundSolutionUsingGurobiRelax = false;
            }
            else
            {
                saveWarmstart(model, vbasis, cbasis);
            }

            validate(model, x, spq, solveOptions);
            sol.x = x;
            sol.w = this->W_q;
            sol.bestRk = this->r;
            sol.isFeasible = isFeasible(this->r);
            sol.bestActiveness = this->activeness;
            sol.bestActivenessNoAbsValue = this->activenessNoAbsValue;

            if (sol.isFeasible)
            {
                high = mid;
                int sense = GRB_MAXIMIZE;
                bool currentIsBetter = isCurrentBetterThanBest(sol.bestRk, sol.w, best, sense);
                if (currentIsBetter)
                {
                    best.bestRk = sol.bestRk;
                    best.w = sol.w;
                    best.isFeasible = sol.isFeasible;
                    best.bestActivenessNoAbsValue = sol.bestActivenessNoAbsValue;
                    best.bestActiveness = sol.bestActiveness;
                    best.x = sol.x;
                }
            }
            else
            {
                low = mid;
                int sense = GRB_MAXIMIZE;
                bool currentIsBetter = isCurrentBetterThanBest(sol.bestRk, sol.w, best, sense);
                if (currentIsBetter)
                {
                    best.bestRk = sol.bestRk;
                    best.w = sol.w;
                    best.isFeasible = sol.isFeasible;
                    best.x = sol.x;
                    best.bestActivenessNoAbsValue = sol.bestActivenessNoAbsValue;
                    best.bestActiveness = sol.bestActiveness;
                }
            }
        }
        else
        {
            low = mid;
        }
    }

    vector<int> finalReducedIds;
    getReducedFromRS(finalReducedIds, best.x, qTarget);
    if(finalReducedIds.size() < qTarget)
    {
        randomTupleSelect(NTuples, qTarget, finalReducedIds, solveOptions.randomSeed);
    }

    formOptions.qSz = finalReducedIds.size();
    formOptions.reducedIds = finalReducedIds;
    formOptions.reduced = true;
    formOptions.activenessNoAbsValue = best.bestActivenessNoAbsValue;
    formOptions.activeness = best.bestActiveness;
    gpro.stop("fetchingToEndOfFirstStage");
}

SolutionMetadata<int> RobustSatisficing::runNaiveFinalStage(std::shared_ptr<StochasticPackageQuery> spq,
                                                       FormulateOptions &formOptions,
                                                       SolveOptions &solveOptions,
                                                       int qTarget)
{
    DecisionVarOptions decVarOptions;
    int qTargetScenarios = solveOptions.hyperParams["reducedScenarios"];
    int exp = 2;

    // Enable variance control
    setDecisionVarOptions(decVarOptions, 0.0, 1.0, 0.0, GrbVarType::Binary);
    formOptions.decisionVarOptions = decVarOptions;
    solveOptions.enableRelaxation = true;
    formOptions.varianceControl = false;
    formOptions.includeObjectiveFunction = true;
    solveOptions.includeObjectiveFunction = true;

    while (true)
    {
        Naive naiveSolver(spq, this->DB_valid);
        formOptions.Z = qTargetScenarios;
        formOptions.randomSeed = solveOptions.randomSeed;
        solveOptions.reduced = formOptions.reduced;
        solveOptions.reducedIds = formOptions.reducedIds;
        solveOptions.foundSolutionUsingGurobiRelax = false;
        // Formulate and solve
        SolutionMetadata<int> sol = naiveSolver.naiveSolve<int>(spq, formOptions, solveOptions);

        gpro.stop("effectiveRuntime");
        double currentTime = gpro.getTime("effectiveRuntime");
        gpro.clock("effectiveRuntime");

        cout << "Naive iteration: qTuples=" << qTarget << ", time=" << currentTime / 1000 << "s" << endl;

        if (sol.x.size() > 0)
        {
            bestSolGlobal.solutions.insert(bestSolGlobal.solutions.end(), sol.solutions.begin(), sol.solutions.end());
            bool isCurrentBetter = isCurrentBetterThanBest(sol.bestRk, sol.w, bestSolGlobal, spq->obj->objSense);
            if(isCurrentBetter)
            {
                // Last solution is always the best
                bestSolGlobal.x = sol.x;
                bestSolGlobal.w = sol.w;
                bestSolGlobal.bestRk = sol.bestRk;
                bestSolGlobal.isFeasible = sol.isFeasible;
                bestSolGlobal.qSz = formOptions.reducedIds.size();
                bestSolGlobal.qScenarios = formOptions.Z;
            }
        }


        bestSolGlobal.terminate = (qTarget >= NTuples);
        gpro.stop("effectiveRuntime");
        bool timeOut = (gpro.getTime("effectiveRuntime") / 1000 > solveOptions.timeout_seconds);
        
        if (timeOut || bestSolGlobal.terminate)
        {
            return bestSolGlobal;
        }
        gpro.clock("effectiveRuntime");
        
        qTarget *= exp;
        qTarget = min(NTuples, qTarget);
        randomTupleSelect(NTuples, qTarget, formOptions.reducedIds, solveOptions.randomSeed);
    }
    return bestSolGlobal;    
}




SolutionMetadata<int> RobustSatisficing::runFinalStage(std::shared_ptr<StochasticPackageQuery> spq,
                                                       FormulateOptions &formOptions,
                                                       SolveOptions &solveOptions,
                                                       vector<int> &finalReducedIds,
                                                       int qTarget)
{
    DecisionVarOptions decVarOptions;
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
        if (sol.x.size() > 0)
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
                randomTupleSelect(NTuples, qTarget, finalReducedIds, solveOptions.randomSeed);
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
                randomTupleSelect(NTuples, qTarget, finalReducedIds, solveOptions.randomSeed);
                deb(finalReducedIds.size());
            }
        }
    }
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

void RobustSatisficing::reduceTuplesStageNoObjConsNoObj(std::shared_ptr<StochasticPackageQuery> spq,
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
    while (high - low > epsBinSearch)
    {
        double mid = low + (high - low) / 2;
        double epsilonStage1 = 1e7;
        formOptions.objCons = false;
        formOptions.iteration = 0;
        DecisionVarOptions decVarOptions;
        setDecisionVarOptions(decVarOptions, 0.0, mid, 0.0, GrbVarType::Continuous);
        deb(low, mid, high);
        formOptions.decisionVarOptions = decVarOptions;
        formOptions.iteration = 0;
        RSFormulator formulator(spq);
        SummarySearch<double> SS(M, spq, epsilonStage1, DB_valid);
        SolutionMetadata<double> sol = SS.summarySearchRS(SS.spq, formulator, formOptions, solveOptions);

        if (sol.x.empty())
        {
            low = mid;
        }
        else if (sol.isFeasible)
        {
            high = mid;
            int sense = GRB_MAXIMIZE;
            bool currentIsBetter = isCurrentBetterThanBest(sol.bestRk, sol.w, best, sense);
            if (currentIsBetter)
            {
                best.bestRk = sol.bestRk;
                best.w = sol.w;
                best.isFeasible = sol.isFeasible;
                best.bestActivenessNoAbsValue = sol.bestActivenessNoAbsValue;
                best.x = sol.x;
            }
        }
        else
        {
            low = mid;
            int sense = GRB_MAXIMIZE;
            bool currentIsBetter = isCurrentBetterThanBest(sol.bestRk, sol.w, best, sense);
            if (currentIsBetter)
            {
                best.bestRk = sol.bestRk;
                best.w = sol.w;
                best.isFeasible = sol.isFeasible;
                best.x = sol.x;
                best.bestActivenessNoAbsValue = sol.bestActivenessNoAbsValue;
            }
        }
    }
    cout << "I FINISHED BINARY SEARCH" << endl;
    vector<int> reducedIds;
    findNonzero(reducedIds, best.x);
    vector<int> unioned;
    getReduced(unioned, solDet, best.x, qTarget);
    deb(reducedIds.size(), unioned.size());
    // formOptions.includeObjectiveFunction = true;
    // solveOptions.includeObjectiveFunction = true;


    cout << "REDUCEDIDS SIZE = " << reducedIds.size() << endl;
    map<int, double> reducedIdsMap;
    vector<int> finalReducedIds;
    populateMapFromVector(reducedIdsMap, reducedIds);
    finalReduce(reducedIdsMap, finalReducedIds, qTarget, solveOptions);
    deb(finalReducedIds, finalReducedIds.size(), reducedIds.size());
    formOptions.qSz = finalReducedIds.size();
    formOptions.reducedIds = finalReducedIds;
    formOptions.reduced = true;
    formOptions.activenessNoAbsValue = best.bestActivenessNoAbsValue;
    formOptions.activeness = best.bestActiveness;
}


vector<vector<pair<int, double>>> RobustSatisficing::findBestObjectiveStage(std::shared_ptr<StochasticPackageQuery> spq,
                                                 FormulateOptions &formOptions,
                                                 SolveOptions &solveOptions, double Z0)
{
    formOptions.Z = cntScenarios;
    cout << "Stage 2" << endl;
    double low = formOptions.low;
    double high = formOptions.high;
    double eps = 1e-5;
    SolutionMetadata<double> bestSolGlobal;
    bestSolGlobal.isFeasible = false;

    // Warmstart basis vectors
    std::vector<int> vbasis;
    std::vector<int> cbasis;

    while (high - low > eps)
    {
        double mid = low + (high - low) / 2; // <- this is our new epsilon
        deb(low, mid, high);
        double Z = (1 - mid) * Z0; // <- this is the new Objective value
        double epsilonStage1 = 1e7;
        formOptions.objValue = Z; // <- this is the new objective value
        formOptions.iteration = 0;

        RSFormulator formulator(spq);
        GRBModel model = formulator.formulate(spq, formOptions);

        // Apply warmstart from previous iteration
        applyWarmstart(model, vbasis, cbasis);

        vector<double> x;
        initializeVector(x, NTuples, 0.0);
        SolutionMetadata<double> sol;
        solve(model, x, solveOptions);

        if (!x.empty())
        {
            // Save warmstart for next iteration
            saveWarmstart(model, vbasis, cbasis);

            validate(model, x, spq, solveOptions);
            sol.x = x;
            sol.w = this->W_q;
            sol.bestRk = this->r;
            sol.isFeasible = isFeasible(this->r);
            sol.bestActiveness = this->activeness;
            sol.bestActivenessNoAbsValue = this->activenessNoAbsValue;

            if (sol.isFeasible) // if feasible this is our current best
            {
                high = mid;
                bestSolGlobal = sol;
            }
            else
            {
                bool isCurrentBetter = isCurrentBetterThanBest(sol.bestRk, sol.w, bestSolGlobal, spq->obj->objSense);
                if (isCurrentBetter)
                {
                    bestSolGlobal = sol;
                }
                low = mid;
            }
        }else
        {
            low = mid;
        }
    }
    return bestSolGlobal.bestActivenessNoAbsValue;
}

vector<int> RobustSatisficing::finalReduce(map<int, double> &reducedIdsMap, vector<int> &reducedIds, int q, SolveOptions &solveOptions)
{
    if (reducedIdsMap.size() < q)
    {
        for (auto &p : reducedIdsMap)
        {
            reducedIds.push_back(p.first);
        }
        randomTupleSelect(NTuples, q, reducedIds, solveOptions.randomSeed);
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
        for (auto &p : reducedIdsMap)
        {
            reducedIds.push_back(p.first);
        }
    }else
    {
        for (auto &p : reducedIdsMap)
        {
            reducedIds.push_back(p.first);
        }
    }
    return reducedIds;
}

SolutionMetadata<int> RobustSatisficing::runFinalStageVarianceControl(std::shared_ptr<StochasticPackageQuery> spq,
                                                                       FormulateOptions &formOptions,
                                                                       SolveOptions &solveOptions,
                                                                       vector<int> &finalReducedIds,
                                                                       int qTarget)
{
    DecisionVarOptions decVarOptions;
    int qInc = 500;
    int exp = 2;
    int qTargetScenarios = qTarget;

    // Enable variance control constraint
    formOptions.varianceControl = true;

    while (true)
    {
        SSFormulator formulatorStage3(spq);
        formOptions.iteration = 0;
        formOptions.Z = min(cntScenarios, qTargetScenarios);
        formOptions.reduced = true;
        formOptions.reducedIds = finalReducedIds;
        solveOptions.reduced = formOptions.reduced;
        solveOptions.reducedIds = formOptions.reducedIds;
        formOptions.activeness = this->activeness;
        formOptions.finalStage = true;
        formOptions.activenessPartition = true;
        formOptions.partitionMostActive = true;
        setDecisionVarOptions(decVarOptions, 0.0, 1.0, 0.0, GrbVarType::Binary);
        formOptions.decisionVarOptions = decVarOptions;
        SummarySearch<int> SS(M, spq, 0.20, DB_valid);
        formOptions.partitionMaximums.clear();
        SolutionMetadata<int> sol = SS.summarySearchStage3(spq, formulatorStage3, formOptions, solveOptions);
        cout << "RETURNED TO WHILE LOOP (VarianceControl)" << endl;

        if (sol.x.size() > 0)
        {
            bool isCurrentBetter = isCurrentBetterThanBest(sol.bestRk, sol.w, bestSolGlobal, spq->obj->objSense);
            if (isCurrentBetter)
            {
                bestSolGlobal.bestRk = sol.bestRk;
                bestSolGlobal.w = sol.w;
                bestSolGlobal.isFeasible = sol.isFeasible;
                bestSolGlobal.x = sol.x;
            }
            bestSolGlobal.solutions.insert(bestSolGlobal.solutions.end(), sol.solutions.begin(), sol.solutions.end());
            double rhs = computeVarianceControlBound(sol.x, 0.95);
            deb(rhs);
        }

        deb(sol.isFeasible, qTarget, qTargetScenarios);

        // Check termination condition: q == NTuples && qScenarios == cntScenarios
        bestSolGlobal.terminate = (qTarget == NTuples) && (qTargetScenarios == cntScenarios);
        gpro.stop("effectiveRuntime");
        bool timeOut = (gpro.getTime("effectiveRuntime") / 1000 > solveOptions.timeout_seconds);
        gpro.clock("effectiveRuntime");

        if (timeOut || bestSolGlobal.terminate)
        {
            cout << "EXITING VARIANCE CONTROL STAGE" << endl;
            bestSolGlobal.qSz = finalReducedIds.size();
            bestSolGlobal.qScenarios = formOptions.Z;
            deb(bestSolGlobal.qSz, bestSolGlobal.qScenarios);
            return bestSolGlobal;
        }

        // Increase q and qScenarios for next iteration
        qTarget += qInc;
        qTargetScenarios += qInc;
        qTarget = min(NTuples, qTarget);
        qTargetScenarios = min(cntScenarios, qTargetScenarios);
        qInc *= exp;
        randomTupleSelect(NTuples, qTarget, finalReducedIds, solveOptions.randomSeed);
        deb(finalReducedIds.size(), qTargetScenarios);
    }
    return bestSolGlobal;
}