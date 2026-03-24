#include "formulator.hpp"
#include <fmt/ranges.h>
#include <boost/algorithm/string/join.hpp>

using namespace std;

Formulator::Formulator() : data(Data::getInstance()) {env.set(GRB_IntParam_OutputFlag, 0);}

Formulator::Formulator(shared_ptr<StochasticPackageQuery> spqPtr) : env(), spq(spqPtr), data(Data::getInstance())
{
    this->DB_optim = spq->tableName;
    this->DB_valid = fmt::format("{}_{}", DB_optim, "validate");
    this->NTuples = pg.getTableSize(spq->tableName);
    this->cntScenarios = pg.getColumnLength(spq->tableName, "profit");
    env.set(GRB_IntParam_OutputFlag, 0);
    env.start();
}

void setDecisionVarOptions(DecisionVarOptions &options, double lb, double ub, double obj, GrbVarType GrbVarType)
{
    options.lb = lb;
    options.ub = ub;
    options.obj = obj;
    options.varType = GrbVarType;
}

GRBVar Formulator::addDecisionVar(GRBModel &model, DecisionVarOptions &options)
{
    char vtype;
    switch (options.varType)
    {
    case GrbVarType::Binary:
        vtype = GRB_BINARY;
        break;
    case GrbVarType::Integer:
        vtype = GRB_INTEGER;
        break;
    case GrbVarType::Continuous:
        vtype = GRB_CONTINUOUS;
        break;
    default:
        vtype = GRB_CONTINUOUS;
    }

    return model.addVar(options.lb, options.ub, options.obj, vtype, options.name);
}

/*
Count Constraint is basically Sum (xi) with bounds
*/
void Formulator::formCountCons(GRBModel &model, shared_ptr<Constraint> cons, GRBVar *xx, FormulateOptions &options)
{
    // cout<<"Formulating a Count Constraint"<<endl;
    shared_ptr<CountConstraint> cntCons = getCount(cons);
    double ub_eps = 1e30;
    double lb_eps = -1e30;

    if (!cntCons)
    {
        return;
    }

    GRBLinExpr cntConsExpr;
    if (options.reduced)
    {
        double coeffval = 1;
        for (int i = 0; i < options.reducedIds.size(); i++)
        {
            // int id = options.reducedIds[i] - 1;
            int id = i;
            cntConsExpr += coeffval * xx[id];
        }
    }
    else
    {
        double coeff[NTuples];
        for (int i = 0; i < NTuples; i++)
        {
            coeff[i] = 1;
        }
        cntConsExpr.addTerms(coeff, xx, NTuples);
    }

    double lb = spq->getValue(cntCons->lb);
    double ub = spq->getValue(cntCons->ub);

    try
    {
        if (ub < ub_eps)
        {
            model.addConstr(cntConsExpr <= ub);
        }
        if (lb > lb_eps)
        {
            model.addConstr(cntConsExpr >= lb);
        }
    }
    catch (GRBException e)
    {
        cout << "Error code 6 = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    // cout<<"Formulated a Count Constraint"<<endl;
}

// select the attr from the DB, and then add the values of each tuple to the corresponding coeff[i]
// formulate GRBLinExpr with addTerms.(coeff,xx,NTuples)
// based on the bounds formulate the constraints ignore the bounds that are not set
void Formulator::formSumCons(GRBModel &model, shared_ptr<Constraint> cons, GRBVar *xx, FormulateOptions &options)
{
    shared_ptr<BoundConstraint> boundCon = dynamic_pointer_cast<BoundConstraint>(cons);
    shared_ptr<AttrConstraint> attrCon = dynamic_pointer_cast<AttrConstraint>(cons);

    double ub_eps = 1e30;
    double lb_eps = -1e30;

    if (!boundCon || !attrCon || attrCon->attrType != numeric_type)
    {
        return;
    }
    GRBLinExpr sumConsExpr;

    auto &detAttr = data.detAttrs[attrCon->attr];
    if (options.reduced)
    {
        for (int i = 0; i < options.reducedIds.size(); i++)
        {
            int id = options.reducedIds[i] - 1;
            double coeffVal = detAttr[id];
            // sumConsExpr += coeffVal * xx[id];
            sumConsExpr += coeffVal * xx[i];
        }
    }
    else
    {
        for (int i = 0; i < NTuples; i++)
        {
            int id = i;
            double coeffVal = detAttr[id];
            sumConsExpr += coeffVal * xx[id];
        }
    }

    double lb = spq->getValue(boundCon->lb);
    double ub = spq->getValue(boundCon->ub);
    // cout << lb << " " << ub << endl;
    try
    {
        if (ub < ub_eps)
        {
            model.addConstr(sumConsExpr, GRB_LESS_EQUAL, ub);
        }
        if (lb > lb_eps)
        {
            model.addConstr(sumConsExpr, GRB_GREATER_EQUAL, lb);
        }
    }
    catch (GRBException &e)
    {
        cout << "Error code 7 = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    // cout << "formulated sum constraint" << endl;
}

/*
check if constraint has bounds, an attribute and is of array_type -> expCons
we can make use of the percomputed table for tuple-wise mean -> determine the table by concatenating "summary"
select cols = concatenate "mean"
assign coeff[NTuples] the tuple-wise means, and then formulate the linearExpr with the coeff and the xx
add constraint based on the bounds
*/
void Formulator::formExpCons(GRBModel &model, shared_ptr<Constraint> cons, GRBVar *xx, FormulateOptions &options)
{
    // cout << "Formulating Exp Constraint" << endl;
    shared_ptr<BoundConstraint> boundCon = dynamic_pointer_cast<BoundConstraint>(cons);
    shared_ptr<AttrConstraint> attrCon = dynamic_pointer_cast<AttrConstraint>(cons);
    if (!boundCon || !attrCon || attrCon->attrType != array_type)
    {
        return;
    }
    // cout << "Expectation sum constraint" << endl;
    double ub_eps = 1e30;
    double lb_eps = -1e30;

    GRBLinExpr expSumConsExpr;

    if (options.reduced)
    {
        for (int i = 0; i < options.reducedIds.size(); i++)
        {
            int id = options.reducedIds[i] - 1;
            // expSumConsExpr += xx[id] * data.stockExpectedProfit[id];
            expSumConsExpr += xx[i] * data.stockExpectedProfit[id];
        }
    }
    else
    {
        for (int i = 0; i < NTuples; i++)
        {
            int id = i;
            expSumConsExpr += xx[id] * data.stockExpectedProfit[id];
        }
    }
    double lb = spq->getValue(boundCon->lb);
    double ub = spq->getValue(boundCon->ub);
    // cout << lb << " " << ub << endl;
    try
    {
        if (ub < ub_eps)
        {
            model.addConstr(expSumConsExpr, GRB_LESS_EQUAL, ub);
        }
        if (lb > lb_eps)
        {
            model.addConstr(expSumConsExpr, GRB_GREATER_EQUAL, lb);
        }
    }
    catch (GRBException &e)
    {
        cout << "Error code 10 = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    // cout << "Formulated Expectation Constraint" << endl;
}

/*fetch the attr from the DB and add the corresponding value to the coeff[i]
then formulate the linearExpression and based on the objSense add Maximize or Minimize
*/
void Formulator::formSumObj(GRBModel &model, shared_ptr<Objective> obj, GRBVar *xx, FormulateOptions &options)
{
    shared_ptr<AttrObjective> attrObj;
    bool isDet = isDeterministic(obj, attrObj);
    if (!isDet || attrObj->objType != numeric_type)
    {
        return;
    }

    auto &detAttr = data.detAttrs[attrObj->obj];
    GRBLinExpr sumObjExpr;

    if (options.reduced)
    {
        for (int i = 0; i < options.reducedIds.size(); i++)
        {
            int id = options.reducedIds[i] - 1;
            // int id = i;
            double coeffVal = detAttr[id];
            sumObjExpr += coeffVal * xx[i];
        }
    }
    else
    {
        for (int i = 0; i < NTuples; i++)
        {
            int id = i;
            double coeffVal = detAttr[id];
            sumObjExpr += coeffVal * xx[id];
        }
    }

    try
    {
        if (obj->objSense == minimize)
        {
            model.setObjective(sumObjExpr, GRB_MINIMIZE);
        }
        else
        {
            model.setObjective(sumObjExpr, GRB_MAXIMIZE);
        }
    }
    catch (GRBException e)
    {
        cout << "Error code 3 = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
}

/*
we can make use of the percomputed table for tuple-wise mean -> determine the table by concatenating "summary"
select cols = concatenate "mean"
assign coeff[NTuples] the tuple-wise means, and then formulate the linearExpr with the coeff and the xx
add objective based on the ObjSense
*/
void Formulator::formExpSumObj(GRBModel &model, shared_ptr<Objective> obj, GRBVar *xx, FormulateOptions &options)
{
    shared_ptr<AttrObjective> attrObj;
    bool isDet = isDeterministic(obj, attrObj);
    if (!isDet || attrObj->objType != array_type)
    {
        return;
    }

    GRBLinExpr expSumObjExpr;
    if (options.reduced)
    {
        for (int i = 0; i < options.reducedIds.size(); i++)
        {
            int id = options.reducedIds[i] - 1;
            // int id = i;
            // expSumObjExpr += xx[id] * data.stockExpectedProfit[id];
            expSumObjExpr += xx[i] * data.stockExpectedProfit[id];
        }
    }
    else
    {
        for (int i = 0; i < NTuples; i++)
        {
            int id = i;
            expSumObjExpr += xx[id] * data.stockExpectedProfit[id];
        }
    }

    try
    {
        if (obj->objSense == minimize)
        {
            model.setObjective(expSumObjExpr, GRB_MINIMIZE);
        }
        else
        {
            model.setObjective(expSumObjExpr, GRB_MAXIMIZE);
        }
    }
    catch (GRBException e)
    {
        cout << "Error code 4 = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
}

/*
check if it's a count objective
assign 1's to each of the coeff and formulate the linearExpression
check what is the objSense and based on that determine whether it's Minimize or Maximize
*/
void Formulator::formCntObj(GRBModel &model, shared_ptr<Objective> obj, GRBVar *xx, FormulateOptions &options)
{
    shared_ptr<CountObjective> cntObj = getCount(obj);
    if (!cntObj)
    {
        return;
    }

    GRBLinExpr cntObjExpr;
    if (options.reduced)
    {
        double coeffval = 1;
        for (int i = 0; i < options.reducedIds.size(); i++)
        {
            // int id = options.reducedIds[i] - 1;
            int id = i;
            cntObjExpr += coeffval * xx[id];
        }
    }
    else
    {
        double coeff[NTuples];
        for (int i = 0; i < NTuples; i++)
        {
            coeff[i] = 1;
        }
        cntObjExpr.addTerms(coeff, xx, NTuples);
    }

    try
    {
        if (obj->objSense == minimize)
        {
            model.setObjective(cntObjExpr, GRB_MINIMIZE);
        }
        else
        {
            model.setObjective(cntObjExpr, GRB_MAXIMIZE);
        }
    }
    catch (GRBException e)
    {
        cout << "Error code 5 = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
}

int getQuantileIdx(double p, int n)
{
    int idx = n * p;
    return idx;
}

// void Formulator::formLCVaR(GRBModel &model, shared_ptr<Constraint> cons, GRBVar *xx, FormulateOptions &options)
// {
//     shared_ptr<ProbConstraint> probCon;
//     shared_ptr<AttrConstraint> attrCon;

//     bool isstoch = isStochastic(cons, probCon, attrCon);
//     if (!isstoch)
//     {
//         return;
//     }
//     // cout << "Formulating CVAR" << endl;
//     double v = spq->getValue(probCon->v);
//     double p = spq->getValue(probCon->p);

//     string selectcols = fmt::format("{},{}_{}", "id", attrCon->attr, "quantiles");
//     string table = fmt::format("{}_{}", DB_optim, "summary");
//     int qtileNumber = pg.getColumnLength(table, fmt::format("{}_{}", attrCon->attr, "quantiles"));
//     // cout << "qtileNumber" << qtileNumber << endl;
//     int qtileIdx = getQuantileIdx(1 - p, qtileNumber);

//     string sql = fmt::format(
//         "SELECT id, "
//         "AVG(unnested_value) AS avg_quantiles "
//         "FROM ( "
//         "  SELECT id, unnest(profit_quantiles[1:{}]) AS unnested_value "
//         "  FROM \"{}\" "
//         ") AS unnested_table "
//         "GROUP BY id",
//         qtileIdx, table);

//     SingleRow sr = SingleRow(sql);
//     GRBLinExpr CVarExpr;

//     if (probCon->psign == Inequality::gteq && probCon->vsign == Inequality::gteq)
//     {
//         Profiler timer;
//         timer.clock("fetching");
//         double coeff[NTuples];
//         while (sr.fetchRow())
//         {
//             int id = sr.getBigInt(0) - 1;
//             double coeffVal = sr.getNumeric(1);
//             coeff[id] = coeffVal;
//         }
//         timer.stop("fetching");
//         fetchRuntime += timer.getTime("fetching");

//         GRBLinExpr CVar_m;
//         CVar_m.addTerms(coeff, xx, NTuples);
//         model.addConstr(CVar_m, GRB_GREATER_EQUAL, v);
//     }
//     else
//     {
//         cout << "Currently There's no Implementation for this combination psign and vsign" << endl;
//     }
//     // cout << "formulated CVaR" << endl;
// }

// figure out which cvar corresponds to each var
void Formulator::formLCVaR(GRBModel &model, shared_ptr<Constraint> cons, GRBVar *xx, FormulateOptions &options)
{
    shared_ptr<ProbConstraint> probCon;
    shared_ptr<AttrConstraint> attrCon;

    bool isstoch = isStochastic(cons, probCon, attrCon);
    if (!isstoch)
    {
        return;
    }
    GRBLinExpr lcvarExpr;
    // cout << "Formulating CVAR" << endl;
    // double v = spq->getValue(probCon->v);
    //double p = spq->getValue(probCon->p);
    double p = 1 - 0.28;
    double v = -2875.0167797102013;


    vector<double> coeffs;
    int quantile = (1 - p) * cntScenarios;
    auto &scenarios = data.stochAttrs[attrCon->attr];
    for (int i = 0; i < NTuples; i++)
    {
        std::vector<double> tmp = scenarios[i];
        std::nth_element(tmp.begin(), tmp.begin() + quantile, tmp.end());
        int j = 0;
        double coeff = 0.0;
        while (true)
        {
            if (j < quantile)
            {
                coeff += tmp[j];
                coeffs.push_back(coeff);
                j++;
            }
            else
            {
                coeff /= j;
                lcvarExpr += coeff * xx[i];
                break;
            }
        }
    }
    model.addConstr(lcvarExpr, GRB_GREATER_EQUAL, v);
}

// sort helper functions for sorting by descending or ascending depending on sign
bool sortbysecDESC(const pair<double, double> &a, const pair<double, double> &b)
{
    return (a.second > b.second);
}

bool sortbysecASC(const pair<double, double> &a, const pair<double, double> &b)
{
    return (a.second < b.second);
}

// function that counts the number of probConstraints given an spq
int countProbConst(shared_ptr<StochasticPackageQuery> spq)
{
    int cnt = 0;
    int cons_num = spq->cons.size();

    for (int i = 0; i < cons_num; i++)
    {
        shared_ptr<ProbConstraint> probCon;
        shared_ptr<AttrConstraint> attrCon;
        bool isStoch = isStochastic(spq->cons[i], probCon, attrCon);
        if (isStoch)
        {
            cnt++;
        }
    }
    return cnt;
}

void Formulator::reshuffleShuffler(std::vector<int> &shuffler)
{
    // cout << "Shuffling" << endl;
    srand(unsigned(time(0)));
    random_shuffle(shuffler.begin(), shuffler.end());
}

void Formulator::populateShuffler(std::vector<int> &v)
{
    for (int j = 0; j < cntScenarios; j++)
    {
        v.push_back(j);
    }
}

void Formulator::createPartitions(shared_ptr<StochasticPackageQuery> spq, FormulateOptions &formOptions)
{
    int cons_num = spq->cons.size();
    int order = 0;
    this->partitions.clear();
    for (int i = 0; i < cons_num; i++)
    {
        shared_ptr<ProbConstraint> probCon;
        shared_ptr<AttrConstraint> attrCon;
        bool isStoch = isStochastic(spq->cons[i], probCon, attrCon);
        if (isStoch)
        {
            double p = spq->getValue(probCon->p);
            if (formOptions.partitionMostActive)
            {
                formOptions.posNegActiveness = false;
                stage3Partition(formOptions, p, order);
            }else if(formOptions.partitionSpreadActiveness)
            {
                formOptions.posNegActiveness = false;
                activenessPartition(formOptions, formOptions.activeness[order], formOptions.innerConstraints[order], p);
            }
            else
            {
                partition(formOptions.Z, formOptions.innerConstraints[order], shuffler);
            }
            order++;
        }
    }
    deb("created partitions for constraints");
}


void Formulator::activenessPartition(FormulateOptions& options, vector<pair<int, double>> &activeness, std::vector<pair<int, double>> &innerConstraints, double p)
{
    int Z = options.Z;
    // 1. Sort activeness in DESCENDING order by the second element (double)
    std::sort(activeness.begin(), activeness.end(), 
        [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
            return a.second > b.second;
        });

    int totalLength = activeness.size();
    int baseSize = totalLength / Z;  // Minimum size per partition
    int remainder = totalLength % Z; // The "+ 1 additional" distribution

    std::vector<std::vector<std::pair<int, double>>> partitionsForConstraint;

    int idx = 0;
    for (int z = 0; z < Z; z++)
    {
        int currentSize = baseSize + (z < remainder ? 1 : 0);
        
        std::vector<std::pair<int, double>> partition;
        partition.reserve(currentSize);
        
        for (int j = 0; j < currentSize; j++)
        {
            int id = activeness[idx].first;
            partition.push_back(innerConstraints[id]);
            idx++;
        }
        
        partitionsForConstraint.push_back(std::move(partition));
    }
    
    // Append to the class member
    partitions.push_back(std::move(partitionsForConstraint));
    deb(partitions[0][0]);
}

// get the vector pairs with realizations and scores
// shuffle them and then find what is partitionSize
// start from beginning, count partitionSize numbers and add them into a partition
// once >= partitionSize reset cnt and start new partition add the prev into partitions
void Formulator::partition(int Z, std::vector<pair<int, double>> &innerConstraints, std::vector<int> &shuffler)
{
    int totalLength = innerConstraints.size();
    int baseSize = totalLength / Z;  // minimum size per partition
    int remainder = totalLength % Z; // number of partitions that get +1 element

    std::vector<vector<pair<int, double>>> partitionsForConstraint;
    reshuffleShuffler(shuffler);    
    int idx = 0;
    for (int z = 0; z < Z; z++)
    {
        int currentSize = baseSize + (z < remainder ? 1 : 0); // distribute remainder
        std::vector<std::pair<int, double>> partition;
        partition.reserve(currentSize);
        for (int j = 0; j < currentSize; j++)
        {
            partition.push_back(innerConstraints[shuffler[idx++]]);
        }
        partitionsForConstraint.push_back(std::move(partition));
    }
    partitions.push_back(partitionsForConstraint);
}

std::set<int> Formulator::getMostActiveScenarios(vector<pair<int, double>> &activeness, int Z, double p)
{
    set<int> mostActiveScenarios;
    sort(activeness.begin(), activeness.end(),
         [](const pair<int, double> &a, const pair<int, double> &b)
         {
             return a.second < b.second;
         });


    int scDimension = min(Z, cntScenarios);
    if (scDimension == cntScenarios)
    {
        for (int i = 0; i < activeness.size(); i++)
        {
            mostActiveScenarios.insert(activeness[i].first);
        }
        return mostActiveScenarios;
    }

    for(int i = 0; i < scDimension; i++)
    {
        mostActiveScenarios.insert(activeness[i].first);
    }
    return mostActiveScenarios;
}

std::set<int> Formulator::getMostActiveScenarios(vector<pair<int, double>> &posActiveness, vector<pair<int, double>> &negActiveness, int Z, double p)
{
    set<int> mostActiveScenarios;
    sort(posActiveness.begin(), posActiveness.end(),
         [](const pair<int, double> &a, const pair<int, double> &b)
         {
             return a.second < b.second;
         });

    sort(negActiveness.begin(), negActiveness.end(),
         [](const pair<int, double> &a, const pair<int, double> &b)
         {
             return a.second > b.second;
         });

    // deb(posActiveness, negActiveness, posActiveness.size(), negActiveness.size());

    int posTarget = 0;
    int negTarget = 0;
    int scDimension = min(Z, cntScenarios);

    if (scDimension == cntScenarios)
    {
        for (int i = 0; i < posActiveness.size(); i++)
        {
            mostActiveScenarios.insert(posActiveness[i].first);
        }
        for (int i = 0; i < negActiveness.size(); i++)
        {
            mostActiveScenarios.insert(negActiveness[i].first);
        }
        return mostActiveScenarios;
    }

    if (posActiveness.size() >= (int)(p * cntScenarios))
    {
        posTarget = scDimension;
    }
    else
    {
        int negativeNeeded = (int)(p * cntScenarios) - posActiveness.size();
        deb(negativeNeeded, scDimension);
        negTarget = min(scDimension, negativeNeeded);
        posTarget = scDimension - negTarget;
    }

    deb(posTarget, negTarget);
    for (int i = 0; i < posTarget; i++)
    {
        mostActiveScenarios.insert(posActiveness[i].first);
    }
    for (int i = 0; i < negTarget; i++)
    {
        mostActiveScenarios.insert(negActiveness[i].first);
    }
    return mostActiveScenarios;
}


std::vector<int> Formulator::getMostActiveScenariosPerPartition(int Z, std::vector<std::pair<int, double>> &sortedActiveness)
{
    std::vector<int> representativeIds;
    representativeIds.reserve(Z);

    int totalLength = sortedActiveness.size();
    int baseSize = totalLength / Z;  // Minimum size per partition
    int remainder = totalLength % Z; // The "+ 1 additional" distribution

    int idx = 0;
    
    for (int z = 0; z < Z; z++)
    {
        int currentSize = baseSize + (z < remainder ? 1 : 0);
        double minAbsValue = std::numeric_limits<double>::max();
        int bestId = -1;

        for (int j = 0; j < currentSize; j++)
        {
            double currentAbs = std::abs(sortedActiveness[idx].second);
            
            if (currentAbs < minAbsValue)
            {
                minAbsValue = currentAbs;
                bestId = sortedActiveness[idx].first;
            }
            idx++;
        }
        representativeIds.push_back(bestId);
    }

    return representativeIds;
}

void Formulator::stage3Partition(FormulateOptions& options,double p, int conOrder)
{
    int Z = options.Z;
    int totalLength = options.innerConstraints[conOrder].size();
    int baseSize = totalLength / Z;  // minimum size per partition
    int remainder = totalLength % Z; // number of partitions that get +1 element

    int idx = 0;
    std::vector<vector<pair<int, double>>> partitionsForConstraint;

    std::set<int> mostActiveScenarios; 
    if(options.posNegActiveness)
    {
        mostActiveScenarios = getMostActiveScenarios(options.posActiveness[conOrder], options.negActiveness[conOrder], Z, p);
    }else
    {
        mostActiveScenarios = getMostActiveScenarios(options.activeness[conOrder], Z, p);
    }
    for (int scenario : mostActiveScenarios)
    {
        std::vector<std::pair<int, double>> partition;
        partition.push_back(options.innerConstraints[conOrder][scenario]);
        partitionsForConstraint.push_back(std::move(partition));
    }
    // deb(mostActiveScenarios, mostActiveScenarios.size());
    // deb(partitionsForConstraint);
    vector<int> shufflerStage3;
    for (int i = 0; i < cntScenarios; i++)
    {
        if (mostActiveScenarios.find(i) != mostActiveScenarios.end())
        {
            continue;
        }
        else
        {
            shufflerStage3.push_back(i);
        }
    }
    reshuffleShuffler(shufflerStage3);
    for (int z = 0; z < Z; z++)
    {
        int currentSize = baseSize + (z < remainder ? 1 : 0); // distribute remainder
        for (int j = 1; j < currentSize; j++)
        {
            partitionsForConstraint[z].push_back(options.innerConstraints[conOrder][shufflerStage3[idx++]]);
        }
    }
    partitions.push_back(partitionsForConstraint);
}

std::vector<std::vector<double>> Formulator::summarize(FormulateOptions &formOptions,
                                                       std::shared_ptr<ProbConstraint> probCon,
                                                       std::shared_ptr<AttrConstraint> attrCon,
                                                       int conOrder)
{
    cout << "summarizing" << endl;
    gpro.clock("sort");
    // we have a vector of doubles, each representing the summary for each partition
    std::vector<std::vector<double>> summaries;
    std::vector<double> summary;
    initializeVectorForm(summary, NTuples, 0.0);
    initializeVectorForm(summaries, formOptions.Z, summary);

    std::set<int> mostActiveScenarios;
    vector<int> mostActiveScenariosIdx;
    if (formOptions.finalStage)
    {
        if(formOptions.posNegActiveness)
        {
            mostActiveScenarios = getMostActiveScenarios(formOptions.posActiveness[conOrder], formOptions.negActiveness[conOrder], formOptions.Z, spq->getValue(probCon->p));
            for (int idx : mostActiveScenarios)
            {
                mostActiveScenariosIdx.push_back(idx);
            }
        }else
        if(formOptions.partitionMostActive)
        {
            mostActiveScenarios = getMostActiveScenarios(formOptions.activeness[conOrder], formOptions.Z, spq->getValue(probCon->p));
            for (int idx : mostActiveScenarios)
            {
                mostActiveScenariosIdx.push_back(idx);
            }
        }else
        if(formOptions.partitionSpreadActiveness)
        {
            mostActiveScenariosIdx = getMostActiveScenariosPerPartition(formOptions.Z, formOptions.activeness[conOrder]);
            cout<<"MOST ACTIVE SCENARIOS PER PARTITION"<<endl;
            deb(mostActiveScenariosIdx);
        }
    }

    // update each partition with the new inner constraint values
    for (int i = 0; i < partitions[conOrder].size(); i++)
    {
        for (int j = 0; j < partitions[conOrder][i].size(); j++)
        {
            int id = partitions[conOrder][i][j].first;
            // abuse that in innerConstraints id = 0 is always first
            partitions[conOrder][i][j].second = formOptions.innerConstraints[conOrder][id].second;
        }
    }
    cout << "updated innerconstraints" << endl;

    for (int i = 0; i < formOptions.Z; i++)
    {
        if (probCon->vsign == Inequality::gteq) //>= min summary
        {
            sort(partitions[conOrder][i].begin(), partitions[conOrder][i].end(), sortbysecDESC);
        }
        else
        {
            sort(partitions[conOrder][i].begin(), partitions[conOrder][i].end(), sortbysecASC);
        }
    }
    gpro.stop("sort");
    auto &scenarios = data.stochAttrs[attrCon->attr];
    gpro.clock("calculate");

    if (formOptions.reduced)
    {
        int reducedSize = formOptions.reducedIds.size();
        gpro.clock("f1");
        for (int j = 0; j < reducedSize; j++)
        {
            gpro.clock("g1");
            // int id = formOptions.reducedIds[i] - 1;
            for (int i = 0; i < partitions[conOrder].size(); i++)
            {
                gpro.clock("s1");
                int id = formOptions.reducedIds[j] - 1;
                double summary;
                // decide when ASC when DESC
                if (probCon->vsign == Inequality::gteq) //>= min summary
                {
                    summary = calculateSummarySmooth(scenarios[id], partitions[conOrder][i], false, formOptions.alpha[conOrder]);
                    if (!mostActiveScenariosIdx.empty())
                    {
                        // summary = POS_INF;
                        double mostActiveCoef = scenarios[id][mostActiveScenariosIdx[i]];
                        if (mostActiveCoef < summary)
                        {
                            summary = mostActiveCoef;
                        }
                    }
                }
                else
                {
                    summary = calculateSummarySmooth(scenarios[id], partitions[conOrder][i], true, formOptions.alpha[conOrder]);
                    if (!mostActiveScenariosIdx.empty())
                    {
                        // summary = NEG_INF;
                        double mostActiveCoef = scenarios[id][mostActiveScenariosIdx[i]];
                        if (mostActiveCoef > summary)
                        {
                            summary = mostActiveCoef;
                        }
                    }
                }
                summaries[i][id] = summary;
                gpro.stop("s1");
            }
            gpro.stop("g1");
        }
        gpro.stop("f1");
    }
    else
    {
        gpro.clock("f2");
        for (int i = 0; i < formOptions.Z; i++) // partitions
        {
            for (int j = 0; j < NTuples; j++) // tuples
            {
                gpro.clock("s2");
                double summary;
                // decide when ASC when DESC
                if (probCon->vsign == Inequality::gteq) //>= min summary
                {
                    summary = calculateSummarySmooth(scenarios[j], partitions[conOrder][i], false, formOptions.alpha[conOrder]);
                }
                else
                {
                    summary = calculateSummarySmooth(scenarios[j], partitions[conOrder][i], true, formOptions.alpha[conOrder]);
                }
                summaries[i][j] = summary;
                gpro.stop("s2");
            }
        }
        gpro.stop("f2");
    }
    if(formOptions.reduced)
    {
        int id = formOptions.reducedIds[0] - 1;
        deb(summaries[0][formOptions.reducedIds[0] - 1],summaries[0][formOptions.reducedIds[1] - 1]);
    }
    cout << "summarized" << endl;
    gpro.stop("calculate");
    return summaries;
}