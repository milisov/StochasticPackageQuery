#include "formulator.hpp"
#include <fmt/ranges.h>
#include <boost/algorithm/string/join.hpp>

using namespace std;

Formulator::Formulator() : data(Data::getInstance()) {}

Formulator::Formulator(shared_ptr<StochasticPackageQuery> spqPtr) : env(), spq(spqPtr), data(Data::getInstance())
{
    this->DB_optim = spq->tableName;
    this->DB_valid = fmt::format("{}_{}", DB_optim, "validate");
    this->NTuples = pg.getTableSize(spq->tableName);
    this->cntScenarios = pg.getColumnLength(spq->tableName, "profit");
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
            //int id = options.reducedIds[i] - 1;
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
            //sumConsExpr += coeffVal * xx[id];
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
            //expSumConsExpr += xx[id] * data.stockExpectedProfit[id];
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
            //int id = i;
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
            //int id = i;
            //expSumObjExpr += xx[id] * data.stockExpectedProfit[id];
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
            //int id = options.reducedIds[i] - 1;
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
    // cout << "Formulating CVAR" << endl;
    double v = spq->getValue(probCon->v);
    double p = spq->getValue(probCon->p);

    string selectcols = fmt::format("{},{}_{}", "id", attrCon->attr, "quantiles");
    string table = fmt::format("{}_{}", DB_optim, "summary");
    int qtileNumber = pg.getColumnLength(table, fmt::format("{}_{}", attrCon->attr, "quantiles"));
    // cout << "qtileNumber" << qtileNumber << endl;
    int qtileIdx = getQuantileIdx(1 - p, qtileNumber);

    string sql = fmt::format(
        "SELECT id, "
        "AVG(unnested_value) AS avg_quantiles "
        "FROM ( "
        "  SELECT id, unnest(profit_quantiles[1:{}]) AS unnested_value "
        "  FROM \"{}\" "
        ") AS unnested_table "
        "GROUP BY id",
        qtileIdx, table);

    SingleRow sr = SingleRow(sql);
    GRBLinExpr CVarExpr;

    if (probCon->psign == Inequality::gteq && probCon->vsign == Inequality::gteq)
    {
        Profiler timer;
        timer.clock("fetching");
        double coeff[NTuples];
        while (sr.fetchRow())
        {
            int id = sr.getBigInt(0) - 1;
            double coeffVal = sr.getNumeric(1);
            coeff[id] = coeffVal;
        }
        timer.stop("fetching");
        fetchRuntime += timer.getTime("fetching");

        GRBLinExpr CVar_m;
        CVar_m.addTerms(coeff, xx, NTuples);
        model.addConstr(CVar_m, GRB_GREATER_EQUAL, v);
    }
    else
    {
        cout << "Currently There's no Implementation for this combination psign and vsign" << endl;
    }
    // cout << "formulated CVaR" << endl;
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

// get the vector pairs with realizations and scores
// shuffle them and then find what is partitionSize
// start from beginning, count partitionSize numbers and add them into a partition
// once >= partitionSize reset cnt and start new partition add the prev into partitions
void Formulator::partition(int Z, std::vector<pair<int, double>> &innerConstraints, std::vector<int> &shuffler, std::vector<std::vector<pair<int, double>>> &partitions)
{
    int totalLength = innerConstraints.size();
    int baseSize = totalLength / Z;  // minimum size per partition
    int remainder = totalLength % Z; // number of partitions that get +1 element

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
        partitions.push_back(std::move(partition));
    }
}

std::vector<std::vector<double>> Formulator::summarize(FormulateOptions &formOptions,
                                                       std::shared_ptr<ProbConstraint> probCon,
                                                       std::shared_ptr<AttrConstraint> attrCon,
                                                       int conOrder)
{
    gpro.clock("sort");
    // we have a vector of doubles, each representing the summary for each partition
    std::vector<std::vector<double>> summaries;
    std::vector<double> summary;
    initializeVectorForm(summary, NTuples, 0.0);
    initializeVectorForm(summaries, formOptions.Z, summary);
    //std::vector<std::vector<pair<int, double>>> partitions;
    for (int i = 0; i < formOptions.Z; i++)
    {
        if (probCon->vsign == Inequality::gteq) //>= min summary
        {
            sort(partitions[i].begin(), partitions[i].end(), sortbysecDESC);
        }
        else
        {
            sort(partitions[i].begin(), partitions[i].end(), sortbysecASC);
        }
    }
    gpro.stop("sort");
    bool init;
    double timeNotFetching = 0.0;
    auto &scenarios = data.stochAttrs[attrCon->attr];
    gpro.clock("calculate");
    
    if(formOptions.reduced)
    {
        int reducedSize = formOptions.reducedIds.size();
        gpro.clock("f1");
        for (int j = 0; j < reducedSize; j++)
        {
            int id = formOptions.reducedIds[j] - 1;
            gpro.clock("g1");
            //int id = formOptions.reducedIds[i] - 1;
            for (int i = 0; i < partitions.size(); i++)
            {
                gpro.clock("s1");
                double summary;
                // decide when ASC when DESC
                if (probCon->vsign == Inequality::gteq) //>= min summary
                {
                    summary = calculateSummary(scenarios[id],partitions[i], false, formOptions.alpha[conOrder]);
                }
                else
                {
                    summary = calculateSummary(scenarios[id],partitions[i], true, formOptions.alpha[conOrder]);
                }
                summaries[i][id] = summary;
                gpro.stop("s1");
            }
            gpro.stop("g1");
        } 
        gpro.stop("f1");
    }else
    {
        gpro.clock("f2");
        for (int i = 0; i < formOptions.Z; i++) //partitions
        {
            for (int j = 0; j < NTuples; j++) //tuples
            {
                gpro.clock("s2");
                double summary;
                // decide when ASC when DESC
                if (probCon->vsign == Inequality::gteq) //>= min summary
                {
                    // deb(i, scenarios[i], partitions[j], formOptions.alpha[conOrder]);
                    summary = calculateSummary(scenarios[j],partitions[i], false, formOptions.alpha[conOrder]);
                }
                else
                {
                    summary = calculateSummary(scenarios[j],partitions[i], true, formOptions.alpha[conOrder]);
                }
                summaries[i][j] = summary;
                gpro.stop("s2");
            }
        }
        gpro.stop("f2");
    }
    gpro.stop("calculate");
    // deb(summaries);
    return summaries;
}