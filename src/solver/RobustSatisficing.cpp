#include "RobustSatisficing.hpp"

void RobustSatisficing::populateMapNonZero(map<int, double> &reducedIdsMap, const vector<double> &sol)
{
    for (int i = 0; i < sol.size(); i++)
    {
        if (sol[i] > 0.0)
        {
            int id = i+1;
            if(reducedIdsMap.find(id) == reducedIdsMap.end())
            {
                reducedIdsMap[id] = data.stockExpectedProfit[id]; // we store the id of the tuple
            }
        }
    }
}

void RobustSatisficing::populateMapFromVector(map<int, double> &reducedIdsMap, const vector<int> &reduced)
{
    for(int i = 0; i < reduced.size(); i++)
    {
        int id = reduced[i];
        if(reducedIdsMap.find(id) == reducedIdsMap.end())
        {
            reducedIdsMap[id] = data.stockExpectedProfit[id];
        }
    }
}

void printActiveConstraints(GRBModel& model) {
    double tol = 1e-40;
    int numConstrs = model.get(GRB_IntAttr_NumConstrs);
    double cnt = 0;
    GRBConstr* constrs = model.getConstrs();

    for (int i = 0; i < numConstrs; i++) {
        double slack = constrs[i].get(GRB_DoubleAttr_Slack);
        std::string name = constrs[i].get(GRB_StringAttr_ConstrName);

        if (std::fabs(slack) < tol) {
            //std::cout << name << " is ACTIVE" << std::endl;
            cnt += 1;
        }
    }
    deb(numConstrs, cnt);
    delete[] constrs;
}

void saveVector(const vector<int> &vec, const string &filename)
{
    ofstream out(filename);
    for (double v : vec)
    {
        out << v << "\n";
    }
    out.close();
}

// Load vector<double> from file
vector<int> loadVector(const string &filename)
{
    vector<int> vec;
    ifstream in(filename);
    int v;
    while (in >> v)
    {
        vec.push_back(v);
    }
    return vec;
}

// Save objective
void saveObjective(double obj, const string &filename)
{
    ofstream out(filename);
    out << obj << "\n";
    out.close();
}

// Load objective
double loadObjective(const string &filename)
{
    ifstream in(filename);
    double obj;
    in >> obj;
    return obj;
}

// SolutionMetadata<int> RobustSatisficing::stochasticDualReducer(std::shared_ptr<StochasticPackageQuery> spq,
//                                                                std::map<std::string, Option> &curveFitOptions, double nudge)
// {
//     RSFormulator formulator(spq);
//     DecisionVarOptions decVarOptions;
//     setDecisionVarOptions(decVarOptions, 0.0, 1.0, 0.0, GrbVarType::Continuous);
//     FormulateOptions formOptions;
//     formOptions.decisionVarOptions = decVarOptions;
//     formOptions.iteration = 0;
//     formulator.formulateBestObjProblem(formulator.modelBestObj, spq, formOptions);
//     std::vector<double> solDet;
//     initializeVector(solDet, NTuples, 0.0);
//     solve(formulator.modelBestObj, solDet);
//     double Z0 = formulator.modelBestObj.get(GRB_DoubleAttr_ObjVal);
//     validate(formulator.modelBestObj, solDet, spq);

//     // //--------------------- STAGE 1 -----------------------------//
//     formOptions.innerConstraints = this->innerConstraints;
//     formOptions.Z = 1;
//     formOptions.Zinit = formOptions.Z;

//     Profiler gpro;
//     string labelStage1 = "Stage1";
//     gpro.clock(labelStage1);

//     vector<int> reducedIds = reduceTuplesStage(spq, formOptions, curveFitOptions, Z0);
//     //saveVector(reducedIds,"/home/fm2288/StochasticPackageQuery/reducedId4.txt");

//     gpro.stop(labelStage1);
//     double totalTimeStage1 = gpro.getTime(labelStage1);

//     formOptions.qSz = reducedIds.size();

//     // //--------------------- STAGE 2 -----------------------------//
//     formOptions.reducedIds = reducedIds;
//     setDecisionVarOptions(decVarOptions, 0.0, 1.0, 0.0, GrbVarType::Continuous);
//     formOptions.decisionVarOptions = decVarOptions;
//     formOptions.iteration = 0;
//     formOptions.innerConstraints = this->innerConstraints;
//     formOptions.Z = 10000;
//     formOptions.Zinit = formOptions.Z;

//     Profiler stopwatchStage2;
//     string labelStage2 = "Stage2";
//     stopwatchStage2.clock(labelStage2);

//     double bestEps = findBestObjectiveStage(spq, formOptions, curveFitOptions, Z0);
//     bestEps *= nudge;

//     stopwatchStage2.stop(labelStage2);
//     double totalTimeStage2 = stopwatchStage2.getTime(labelStage2);
//     //------------------- FINAL ----------------------------//
//     //vector<int>reducedIds = loadVector("/home/fm2288/StochasticPackageQuery/reducedIds.txt");
//     // double bestEps = loadObjective("/home/fm2288/StochasticPackageQuery/bestObj.txt");
//     formOptions.reducedIds = reducedIds;
//     setDecisionVarOptions(decVarOptions, 0.0, 1.0, 0.0, GrbVarType::Binary);
//     formOptions.decisionVarOptions = decVarOptions;
//     formOptions.iteration = 0;
//     formOptions.innerConstraints = this->innerConstraints;
//     formOptions.Z = 10000;
//     formOptions.Zinit = formOptions.Z;
//     formOptions.objValue = (1 - bestEps) * Z0;
//     formOptions.reduced = true;
//     formOptions.RS = false;
//     SummarySearch SS(M, spq, 1);
//     SolutionMetadata<int> sol = SS.summarySearchRS<int>(SS.spq, formulator, formOptions, curveFitOptions, 1);
//     shared_ptr<AttrObjective> attrObj;
//     bool isDet = isDeterministic(spq->obj, attrObj);
//     double obj;
//     if (sol.x.size() > 0)
//     {
//         obj = calculateExpSumObj(sol.x, attrObj);
//     }
//     sol.qSz = formOptions.qSz;
//     sol.objConsValue = obj;
//     sol.timeStage1 = totalTimeStage1;
//     sol.timeStage2 = totalTimeStage2;
//     sol.bestEps = bestEps;
//     return sol;
// }

SolutionMetadata<int> RobustSatisficing::stochasticDualReducer(std::shared_ptr<StochasticPackageQuery> spq,
                                                               std::map<std::string, Option> &curveFitOptions, double nudge)
{
    RSFormulator formulator(spq);
    DecisionVarOptions decVarOptions;
    setDecisionVarOptions(decVarOptions, 0.0, 1.0, 0.0, GrbVarType::Continuous);
    FormulateOptions formOptions;
    formOptions.decisionVarOptions = decVarOptions;
    formOptions.iteration = 0;

    formulator.formulateBestObjProblem(formulator.modelBestObj, spq, formOptions);
    std::vector<double> solDet;
    initializeVector(solDet, NTuples, 0.0);
    SolveOptions options;
    options.reduced = formOptions.reduced;
    options.reducedIds = formOptions.reducedIds;
    solve(formulator.modelBestObj, solDet, options);
    double Z0 = formulator.modelBestObj.get(GRB_DoubleAttr_ObjVal);
    validate(formulator.modelBestObj, solDet, spq, options);
    
    //this map contains the reducedIds and the expected profit for each id
    map<int, double> reducedIdsMap;
    populateMapNonZero(reducedIdsMap, solDet);
    deb(reducedIdsMap, reducedIdsMap.size());

    // //--------------------- STAGE 1 -----------------------------//
    formOptions.innerConstraints = this->innerConstraints;
    formOptions.Z = 1;
    formOptions.Zinit = formOptions.Z;

    string labelStage1 = "Stage1";
    gpro.clock(labelStage1);

    //vector<int> reducedIds = reduceTuplesStage(spq, formOptions, curveFitOptions, Z0);
    vector<int> reducedIds = reduceTuplesStageNoObjCons(spq, formOptions, curveFitOptions);

    populateMapFromVector(reducedIdsMap, reducedIds);
    deb(reducedIds,reducedIdsMap.size());
    vector<int> finalReducedIds;
    map<int, double> reducedIdsMap2 = reducedIdsMap; 
    finalReduce(reducedIdsMap, finalReducedIds, 500);
    gpro.stop(labelStage1);

    for(auto &it : reducedIdsMap2)
    {
        if(reducedIdsMap.find(it.first) == reducedIdsMap.end())
        {
            deb("FOUND ID NOT IN REDUCED IDS: " + to_string(it.first));
        }
    }
    formOptions.cbasis.clear();
    formOptions.vbasis.clear();
    formOptions.qSz = finalReducedIds.size();
    //--------------------- STAGE 2 -----------------------------//
    formOptions.Z = min(cntScenarios, formOptions.qSz);
    formOptions.Zinit = formOptions.Z;
    formOptions.computeActiveness = true;
    RSFormulator formulatorSAA(spq);
    formOptions.low = 0.0;
    formOptions.high = 1.0;
    formOptions.SAA = false;
    formOptions.RS = true;
    formOptions.objCons = true;
    formOptions.reducedIds = finalReducedIds;
    formOptions.reduced = true;
    setDecisionVarOptions(decVarOptions, 0.0, 1.0, 0.0, GrbVarType::Continuous);
    formOptions.decisionVarOptions = decVarOptions;
    formOptions.iteration = 0;
    formOptions.objValue = NEG_INF;
    string labelStage2 = "Stage2";
    gpro.clock(labelStage2);
    // GRBModel modelSAA = formulatorSAA.formulate(spq, formOptions);
    // double bestEps = findBestObjectiveStageSAA(modelSAA, spq, formulator, formOptions, curveFitOptions, Z0);

    double bestEps = findBestObjectiveStage(spq, formOptions, curveFitOptions, Z0);
    deb(bestEps);
    //deb(this->bestActivenessRS);
    gpro.stop(labelStage2);
    // ------------------- FINAL ----------------------------//
    formOptions.reduced = true;
    formOptions.reducedScenarios = true;
    formOptions.posActiveness = this->bestPosActivenessRS;
    formOptions.negActiveness = this->bestNegActivenessRS;
    formOptions.computeActiveness = true;
    formOptions.reducedIds = finalReducedIds;
    setDecisionVarOptions(decVarOptions, 0.0, 1.0, 0.0, GrbVarType::Binary);
    formOptions.decisionVarOptions = decVarOptions;
    Naive naiveSolver(spq);
    SolutionMetadata<int> sol = naiveSolver.solveNaive<int>(spq, formOptions);
    sol.qSz = formOptions.qSz;
    sol.timeStage1 = 0;
    sol.timeStage2 = 0;
    gpro.print();
    sol.bestEps = bestEps;
    return sol;

    // RSFormulator formulatorSAAILP(spq);
    // formOptions.SAA = true;
    // formOptions.objCons = true;
    // formOptions.reduced = true;
    // formOptions.reducedIds = finalReducedIds;
    // setDecisionVarOptions(decVarOptions, 0.0, 1.0, 0.0, GrbVarType::Binary);
    // formOptions.decisionVarOptions = decVarOptions;
    // formOptions.iteration = 0;
    // formOptions.objValue = NEG_INF;
    // formOptions.innerConstraints = this->innerConstraints;
    // GRBModel modelILP = formulatorSAAILP.formulate(spq, formOptions);
    // SolutionMetadata<int> sol = finalStageILP(modelILP, spq, formulatorSAAILP, formOptions, bestEps, Z0);
    // sol.qSz = formOptions.qSz;
    // sol.timeStage1 = totalTimeStage1;
    // sol.timeStage2 = totalTimeStage2;
    // sol.bestEps = bestEps;
    // return sol;
}


std::vector<int> RobustSatisficing::reduceTuplesStageNoObjCons(std::shared_ptr<StochasticPackageQuery> spq,
                                                         FormulateOptions &formOptions,
                                                         std::map<std::string, Option> &curveFitOptions)
{
    vector<int>reducedIds;
    double epsilonStage1 = 1e7;
    formOptions.RS = true;
    formOptions.objCons = false;
    formOptions.iteration = 0;
    RSFormulator formulator(spq);
    SummarySearch SS(M, spq, epsilonStage1);
    SolutionMetadata<double> sol = SS.summarySearchRS<double>(SS.spq, formulator, formOptions, curveFitOptions, 1);
    findNonzero(reducedIds, sol.x);
    return reducedIds;
}

std::vector<int> RobustSatisficing::reduceTuplesStage(std::shared_ptr<StochasticPackageQuery> spq,
                                                      FormulateOptions &formOptions,
                                                      std::map<std::string, Option> &curveFitOptions, double Z0)
{
    int qLeft = 400;
    int qRight = 600;
    // doing binary search on epsilon, i.e left = 0 and right = 1
    double low = 0.0;
    double high = 1.0;
    double eps = 1e-5;
    vector<int> reducedIds;
    
    int totalSystems = 0;
    int totalSolveTime = 0;

    auto summarySearchStage1 = [&](double eps)
    {
        double Z = (1 - eps) * Z0; // <- this is the new Objective value
        double epsilonStage1 = 1e7;
        formOptions.RS = true;
        formOptions.objValue = Z; // <- this is the new objective value
        formOptions.iteration = 0;

        
        RSFormulator formulator(spq);
        SummarySearch SS(M, spq, epsilonStage1);
        SolutionMetadata<double> sol = SS.summarySearchRS<double>(SS.spq, formulator, formOptions, curveFitOptions, 1);
        totalSystems += (formOptions.iteration+1);
        totalSolveTime += SS.timeSolve;
        return sol;
    };

    SolutionMetadata<double> sol = summarySearchStage1(high);
    while (true)
    {
        shared_ptr<AttrObjective> attrObj;
        bool isDet = isDeterministic(spq->obj, attrObj);
        double obj = calculateExpSumObj(sol.x, attrObj);
        if (obj > formOptions.objValue)
        {
            high = 1 - obj / Z0;
            break;
        }
        low = high;
        high *= 2.0;
        sol = summarySearchStage1(high);
    }
    formOptions.low = low;
    formOptions.high = high;
    while (high - low > eps)
    {
        reducedIds.clear();
        double mid = low + (high - low) / 2; // <- this is our new epsilon
        SolutionMetadata<double> sol = summarySearchStage1(mid);
        findNonzero(reducedIds, sol.x);
        if(sol.x.size() > 0)
        {
            high = mid;   
        }else
        {
            low = mid;
        }
        // if (reducedIds.size() < qLeft)
        // {
        //     low = mid;
        // }
        // else if (reducedIds.size() > qRight)
        // {
        //     high = mid;
        // }
        // else
        // {
        //     deb(reducedIds.size());
        //     return reducedIds; // <- this is the reduced ids
        // }
    }
    // deb(totalSystems, timeSolve);
    return reducedIds; //<--the best we could do
}

double RobustSatisficing::findBestObjectiveStage(std::shared_ptr<StochasticPackageQuery> spq,
                                                 FormulateOptions &formOptions,
                                                 std::map<std::string, Option> &curveFitOptions, double Z0)
{
    int totalSystems = 0;
    int totalSolveTime = 0;

    double low = formOptions.low;
    double high = formOptions.high;
    double eps = 1e-5;
    double best = -1; // a dummy for epsilon,
    double bestRk = 0;
    while (high - low > eps)
    {
        double mid = low + (high - low) / 2; // <- this is our new epsilon
        double Z = (1 - mid) * Z0;           // <- this is the new Objective value
        double epsilonStage1 = 1e7;
        formOptions.RS = true;
        formOptions.reduced = true;
        formOptions.objValue = Z; // <- this is the new objective value
        formOptions.iteration = 0;
        RSFormulator formulator(spq);
        SummarySearch SS(M, spq, epsilonStage1);
        SolutionMetadata<double> sol = SS.summarySearchRS<double>(SS.spq, formulator, formOptions, curveFitOptions, 1);
        deb(bestRk, sol.bestRk, posActivenessRS.size(), negActivenessRS.size());
        if(sol.isFeasible)
        {
            best = mid;
            high = mid;
            this->bestPosActivenessRS = sol.bestPosActivenessRS;
            this->bestNegActivenessRS = sol.bestNegActivenessRS;
            bestRk = sol.bestRk;
        }else
        if (sol.bestRk > bestRk)
        {
            this->bestPosActivenessRS = sol.bestPosActivenessRS;
            this->bestNegActivenessRS = sol.bestNegActivenessRS;
            bestRk = sol.bestRk;
        }
        else
        {
            low = mid;
        }
    }
    deb(bestRk, bestPosActivenessRS.size(), bestNegActivenessRS.size());
    if (best == -1)
    {
        best = low + (high - low) / 2;
    }
    return best;
}

double RobustSatisficing::findBestObjectiveStageSAA(GRBModel &model, std::shared_ptr<StochasticPackageQuery> spq,
                                                    RSFormulator &formulator,
                                                    FormulateOptions &formOptions,
                                                    std::map<std::string, Option> &curveFitOptions, double Z0)
{
    model.write("/home/fm2288/StochasticPackageQuery/try2.lp");
    deb("STAGE 2");
    double low = formOptions.low;
    double high = formOptions.high;
    double eps = 1e-5;
    double best = -1; // a dummy for epsilon,

    Profiler stopwatch;
    string label = "solver";
    double totalTimeSolver = 0.0;


    double Z = (1 - high) * Z0;
    formOptions.RS = true;
    formOptions.reduced = true;
    formOptions.objValue = Z;
    formOptions.iteration = 0;
    GRBVar *xx = model.getVars();
    formulator.formObjCons(model, spq->obj, xx, formOptions);
    gpro.clock(label);
    SolutionMetadata<double> sol = solveSAA<double>(model, formOptions);
    gpro.stop(label);
    shared_ptr<AttrObjective> attrObj;
    bool isDet = isDeterministic(spq->obj, attrObj);
    double obj = calculateExpSumObj(sol.x, attrObj);


    while (true)
    {
        shared_ptr<AttrObjective> attrObj;
        bool isDet = isDeterministic(spq->obj, attrObj);
        double obj = calculateExpSumObj(sol.x, attrObj);
        if (obj > formOptions.objValue)
        {
            high = 1 - obj / Z0;
            break;
        }
        low = high;
        high *= 2.0;
        Z = (1 - high) * Z0;      // <- this is the new Objective value
        formOptions.objValue = Z; // <- this is the new objective value
        formOptions.iteration = 0;
        formulator.formObjCons(model, spq->obj, xx, formOptions);
        sol = solveSAA<double>(model, formOptions);
    }

    while (high - low > eps)
    {
        double mid = low + (high - low) / 2; // <- this is our new epsilon
        deb(low, mid, high);
        double Z = (1 - mid) * Z0; // <- this is the new Objective value
        formOptions.RS = true;
        formOptions.reduced = true;
        formOptions.objValue = Z; // <- this is the new objective value
        formOptions.iteration = 0;
        GRBVar *xx = model.getVars();
        formulator.formObjCons(model, spq->obj, xx, formOptions);
        SolutionMetadata<double> sol = solveSAA<double>(model, formOptions);
        if (sol.x.size() > 0)
        {
            for (int i = 0; i < NTuples; i++)
            {
                xx[i].set(GRB_DoubleAttr_Start, sol.x[i]);
            }
        }
        if (sol.isFeasible)
        {
            //printActiveConstraints(model);
            best = mid;
            high = mid;
        }
        else
        {
            low = mid;
        }
    }
    if (best == -1)
    {
        best = low + (high - low) / 2;
    }
    deb("Total time solver: ", totalTimeSolver);
    //printActiveConstraints(model);
    double interval = high - low;
    return best;
}

SolutionMetadata<int> RobustSatisficing::finalStageILP(GRBModel &model, std::shared_ptr<StochasticPackageQuery> spq,
                                                       RSFormulator &formulator,
                                                       FormulateOptions &formOptions, double bestEps, double Z0)
{
    deb("Final Stage");
    deb(bestEps);
    double low = bestEps;
    double high = low;
    formOptions.SAA = true;
    formOptions.reduced = true;
    GRBVar *xx = model.getVars();
    SolutionMetadata<int> bestSol;
    while (true)
    {
        formOptions.objValue = (1 - high) * Z0;
        formulator.formObjCons(model, spq->obj, xx, formOptions);
        SolutionMetadata<int> sol = solveSAA<int>(model, formOptions);
        // if there is solution, we need to stop the exponential search
        // go for binary search on the new interval
        if (sol.x.size() > 0)
        {
            for (int i = 0; i < NTuples; i++)
            {
                xx[i].set(GRB_DoubleAttr_Start, sol.x[i]);
            }
            if (high == low)
            {
                shared_ptr<AttrObjective> attrObj;
                bool isDet = isDeterministic(spq->obj, attrObj);
                double obj = calculateExpSumObj(sol.x, attrObj);
                bestSol = sol;
                bestSol.objConsValue = obj;
            }
            break;
        }
        low = high;
        high *= 10;
    }

    double eps = 1e-5;
    while (high - low > eps)
    {
        double mid = low + (high - low) / 2; // <- this is our new epsilon
        double Z = (1 - mid) * Z0;           // <- this is the new Objective value
        formOptions.objValue = Z;
        int numCon = model.get(GRB_IntAttr_NumConstrs);
        deb(low, mid, high, numCon, Z);

        formulator.formObjCons(model, spq->obj, xx, formOptions);
        SolutionMetadata<int> sol = solveSAA<int>(model, formOptions);
        if (sol.x.size() > 0)
        {
            for (int i = 0; i < NTuples; i++)
            {
                xx[i].set(GRB_DoubleAttr_Start, sol.x[i]);
            }
        }
        if (sol.isFeasible) // if there is a solution
        {
            shared_ptr<AttrObjective> attrObj;
            bool isDet = isDeterministic(spq->obj, attrObj);
            double obj = calculateExpSumObj(sol.x, attrObj);
            bestSol = sol;
            bestSol.objConsValue = obj;
            high = mid;
        }
        else
        {
            low = mid;
        }
    }
    return bestSol;
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
    else if  (reducedIdsMap.size() > q)
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