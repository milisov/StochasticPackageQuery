#include "SummarySearch.hpp"
#include <boost/algorithm/string/join.hpp>
#include <fmt/ranges.h>

// sort helper functions for sorting by descending or ascending depending on sign
bool sortbysecDESC(const pair<double, double> &a, const pair<double, double> &b)
{
    return (a.second > b.second);
}

bool sortbysecASC(const pair<double, double> &a, const pair<double, double> &b)
{
    return (a.second < b.second);
}

void printConstraints(GRBModel &model)
{
    GRBConstr *constrs = model.getConstrs();
    int numConstrs = model.get(GRB_IntAttr_NumConstrs);

    cout << "NUMBER OF CONSTRAINTS ADDED: " << numConstrs << endl;

    for (int i = 0; i < numConstrs; i++)
    {
        int maxCoef = 0;
        int minCoef = 9999999;
        cout << endl;
        GRBConstr constr = constrs[i];
        string constrName = constr.get(GRB_StringAttr_ConstrName);
        GRBLinExpr expr = model.getRow(constr);
        double rhs = constr.get(GRB_DoubleAttr_RHS);
        char sense = constr.get(GRB_CharAttr_Sense);

        cout << "Constraint " << constrName << ": ";
        for (int j = 0; j < expr.size(); j++)
        {
            GRBVar var = expr.getVar(j);
            double coeff = expr.getCoeff(j);
            if (coeff > maxCoef)
            {
                maxCoef = coeff;
            }
            if (coeff < minCoef)
            {
                minCoef = coeff;
            }
            cout << coeff << "*" << var.get(GRB_StringAttr_VarName);
            if (j < expr.size() - 1)
            {
                cout << " + ";
            }
        }
        cout << " ";

        switch (sense)
        {
        case GRB_LESS_EQUAL:
            cout << "<= ";
            break;
        case GRB_EQUAL:
            cout << "= ";
            break;
        case GRB_GREATER_EQUAL:
            cout << ">= ";
            break;
        }

        cout << rhs << endl;
        cout << "MAX =" << maxCoef << " MIN = " << minCoef << endl;
    }
    int numGenConstrs = model.get(GRB_IntAttr_NumGenConstrs);
    cout << "NUMBER OF GENERAL CONSTRAINTS ADDED: " << numGenConstrs << endl;

    GRBGenConstr *genConstrs = model.getGenConstrs();
    for (int i = 0; i < numGenConstrs; i++)
    {
        int minGen = 999999;
        int maxGen = 0;
        GRBGenConstr genc = genConstrs[i];
        int type = genc.get(GRB_IntAttr_GenConstrType);
        cout << endl;
        if (type == GRB_GENCONSTR_INDICATOR)
        {
            try
            {
                GRBVar binvar;
                int binval;
                GRBLinExpr expr;
                char sense;
                double rhs;

                model.getGenConstrIndicator(genc, &binvar, &binval, &expr, &sense, &rhs);

                cout << "Indicator Constraint " << genc.get(GRB_StringAttr_GenConstrName) << ": ";
                cout << binvar.get(GRB_StringAttr_VarName) << " == " << binval << " -> ";

                for (int j = 0; j < expr.size(); j++)
                {
                    GRBVar var = expr.getVar(j);
                    double coeff = expr.getCoeff(j);
                    if (coeff > maxGen)
                    {
                        maxGen = coeff;
                    }
                    if (coeff < minGen)
                    {
                        minGen = coeff;
                    }
                    cout << coeff << "*" << var.get(GRB_StringAttr_VarName);
                    if (j < expr.size() - 1)
                    {
                        cout << " + ";
                    }
                }
                cout << " ";

                switch (sense)
                {
                case GRB_LESS_EQUAL:
                    cout << "<= ";
                    break;
                case GRB_EQUAL:
                    cout << "= ";
                    break;
                case GRB_GREATER_EQUAL:
                    cout << ">= ";
                    break;
                }

                cout << rhs << endl;
            }
            catch (GRBException &e)
            {
                cout << "Error code 1 = " << e.getErrorCode() << endl;
                cout << e.getMessage() << endl;
            }
            cout << "GENERATOR MAX =" << maxGen << "MIN = " << minGen << endl;
        }
    }
}

void printObjective(GRBModel &model)
{
    cout << "I am printing the objective" << endl;
    try
    {
        int sense = model.get(GRB_IntAttr_ModelSense);
        std::string senseStr = (sense == GRB_MINIMIZE) ? "Minimize" : "Maximize";

        // Print the sense of the objective
        std::cout << senseStr << " ";

        // Get the linear part of the objective
        GRBQuadExpr quadExpr = model.getObjective();
        GRBLinExpr linearObj = quadExpr.getLinExpr();

        int numVars = linearObj.size();
        bool isFirstTerm = true;

        // Print linear terms
        for (int i = 0; i < numVars; ++i)
        {
            GRBVar var = linearObj.getVar(i);
            double coeff = linearObj.getCoeff(i);

            if (coeff != 0.0)
            {
                if (!isFirstTerm && coeff > 0)
                {
                    std::cout << " + ";
                }
                else if (coeff < 0)
                {
                    std::cout << " - ";
                    coeff = -coeff;
                }

                if (coeff != 1.0)
                {
                    std::cout << coeff << "*";
                }

                std::cout << var.get(GRB_StringAttr_VarName);
                isFirstTerm = false;
            }
        }

        std::cout << std::endl;
    }
    catch (GRBException &e)
    {
        std::cout << "Error code 2 = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
}

void printExpression(const GRBLinExpr &expr)
{
    cout << "Expression: ";
    for (int j = 0; j < expr.size(); j++)
    {
        try
        {
            GRBVar var = expr.getVar(j);
            double coeff = expr.getCoeff(j);
            cout << coeff << "*" << var.get(GRB_StringAttr_VarName);
            if (j < expr.size() - 1)
            {
                cout << " + ";
            }
        }
        catch (GRBException &e)
        {
            std::cout << "Error code 200 = " << e.getErrorCode() << std::endl;
            std::cout << e.getMessage() << std::endl;
        }
    }
    cout << endl;
}

void printVariableNames(GRBModel &model)
{
    try
    {
        // Get all the variables in the model
        int numVars = model.get(GRB_IntAttr_NumVars);
        GRBVar *vars = model.getVars();

        // Print the names of each variable
        for (int i = 0; i < numVars; i++)
        {
            GRBVar var = vars[i];
            std::string varName = var.get(GRB_StringAttr_VarName);
            std::cout << "Variable " << i << ": " << varName << std::endl;
        }

        delete[] vars;
    }
    catch (GRBException &e)
    {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
    catch (...)
    {
        std::cout << "Exception during printing variable names" << std::endl;
    }
}

// input: sorted partition, boolean that determines if we need min/max summary, alpha value to consider
// G(alpha) = |alpha*Пz| elements to consider in order to find min/max
template <typename T1, typename T2>
T1 calculateSummary(std::vector<pair<T1, T2>> &vect, bool maxS, double alpha)
{
    int G_alpha = (int)ceil(alpha * (double)vect.size());
    T1 summary = vect[0].first;
    for (int i = 1; i < G_alpha; i++)
    {
        if (maxS)
        {
            if (vect[i].first > summary)
            {
                summary = vect[i].first;
            }
        }
        else
        {
            if (vect[i].first < summary)
            {
                summary = vect[i].first;
            }
        }
    }
    return summary;
}

// function that initializes given vector
template <typename T>
void initializeVector(std::vector<T> &v, int sz, T init)
{
    for (int i = 0; i < sz; i++)
    {
        v.push_back(init);
    }
}

// function that finds the indices of solution x that are non-zero
template <typename T>
void findNonzero(std::vector<int> &selectIds, std::vector<T> &x)
{
    for (int i = 0; i < x.size(); i++)
    {
        if (x[i] != 0)
        {
            selectIds.push_back(i + 1);
        }
    }
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
// for each of the Stochastic Constraints we need to check how many are satisfied
// get the string sign from each prob constraint and perform the right operation
// the value of v is retrieved from the spq
int countSatisfied(int numScenarios, std::vector<double> &innerConst, shared_ptr<ProbConstraint> probCon, shared_ptr<StochasticPackageQuery> spq)
{
    int satisfied = 0;
    double epsilon = 1e-7;
    // cout << numScenarios << endl;
    for (int j = 0; j < numScenarios; j++)
    {
        if (probCon->vsign == Inequality::gteq) // >=
        {
            if (innerConst[j] >= spq->getValue(probCon->v) - epsilon)
            {
                satisfied++;
            }
        }
        else
        {
            if (innerConst[j] <= spq->getValue(probCon->v) + epsilon) // <=
            {
                satisfied++;
            }
        }
    }

    return satisfied;
}

// If the constraint is of form P(sth)>=p rk = satisfied/numScenarios - p
// else we need to convert the constraint in the right format

double calculateRk(shared_ptr<ProbConstraint> probCon, int satisfied, int numScenarios, shared_ptr<StochasticPackageQuery> spq)
{
    double rk;
    cout << "Satisfied: " << satisfied << endl;
    if (probCon->psign == Inequality::gteq)
    {
        rk = (double)satisfied / numScenarios - spq->getValue(probCon->p);
        cout << "rk-gteq" << " " << rk << endl;
    }
    else
    {
        rk = (1 - (double)satisfied / numScenarios) - (1 - spq->getValue(probCon->p));
    }

    return rk;
}

template <typename T>
double SummarySearch::calculateCntObj(std::vector<T> &x, std::vector<int> &selectIds, shared_ptr<CountObjective> cntObj, string DB_optim)
{
    cout << "Calculating Count Objective" << endl;
    double sum = 0;
    for (int i = 0; i < x.size(); i++)
    {
        sum += x[i];
    }
    return sum;
}

template <typename T>
double SummarySearch::calculateSumObj(std::vector<T> &x, std::vector<int> &selectIds, shared_ptr<AttrObjective> attrObj, string DB_optim)
{
    cout << "Calculating Sum Objective" << endl;
    double sum = 0;

    string selectcols = fmt::format("{},{}", "id", attrObj->obj);
    string sql = fmt::format("SELECT {} FROM \"{}\" WHERE id = ANY(ARRAY[{}])", selectcols, DB_optim, fmt::join(selectIds, ","));

    SingleRow sr = SingleRow(sql);
    int i = 0;
    while (sr.fetchRow())
    {
        double si = sr.getNumeric(1);
        i = sr.getBigInt(0) - 1;
        sum += si * x[i];
    }
    return sum;
}

template <typename T>
double SummarySearch::calculateExpSumObj(std::vector<T> &x, std::vector<int> &selectIds, shared_ptr<AttrObjective> attrObj, string DB_optim)
{
    cout << "Calculating Expected Sum Objective" << endl;
    double sum = 0;

    string table = fmt::format("{}_{}", DB_optim, "summary");
    string objMean = fmt::format("{}_{}", attrObj->obj, "mean");
    string selectcols = fmt::format("{},{}", "id", objMean);
    string sql = fmt::format("SELECT {} FROM \"{}\" WHERE id = ANY(ARRAY[{}])", selectcols, table, fmt::join(selectIds, ","));

    SingleRow sr = SingleRow(sql);
    int i = 0;
    while (sr.fetchRow())
    {
        double si = sr.getNumeric(1);
        i = sr.getBigInt(0) - 1;
        sum += si * x[i];
    }
    return sum;
}

template <typename T>
double SummarySearch::calculateObj(std::vector<T> &x, std::vector<int> &selectIds, shared_ptr<Objective> obj, string DB_optim)
{
    double w = 0;
    shared_ptr<CountObjective> cntObj = getCount(obj);
    if (cntObj)
    {
        w = calculateCntObj(x, selectIds, cntObj, DB_optim);
        return w;
    }

    shared_ptr<AttrObjective> attrObj;
    bool isDet = isDeterministic(obj, attrObj);
    if (isDet && attrObj->objType == numeric_type)
    {
        w = calculateSumObj(x, selectIds, attrObj, DB_optim);
        return w;
    }

    shared_ptr<AttrObjective> attrObj2;
    bool isDet2 = isDeterministic(obj, attrObj2);
    if (isDet2 && attrObj2->objType == array_type)
    {
        w = calculateExpSumObj(x, selectIds, attrObj2, DB_optim);
        return w;
    }
    return -1;
}

template <typename T>
void SummarySearch::validate(std::vector<T> &x, shared_ptr<StochasticPackageQuery> spq, int M_hat, string DB_optim)
{
    // empty previous surplus vector
    cout << "Validating!" << endl;
    // clear r from previous iteration -> holds all of the rk values for each iteration
    this->r.clear();
    // we need to clear the innerConstraints in order to update them with new values
    this->innerConstraints.clear();
    // because we get row by row we need a way to compute the innerConst for all scenarios in the same time
    int cons_num = spq->cons.size();

    std::vector<int> selectIds;
    findNonzero(selectIds, x);
    // for all of the constrants, check if they are probConstraints
    // query the optim table
    // for each tuple update innerConst
    for (int i = 0; i < cons_num; i++)
    {
        shared_ptr<ProbConstraint> probCon;
        shared_ptr<AttrConstraint> attrCon;
        bool isStoch = isStochastic(spq->cons[i], probCon, attrCon);
        if (isStoch)
        {
            std::vector<double> innerConst;
            int row;
            int numScenarios = 0;
            bool initialized;
            // initialize innerConst with 0 in order to computer the innerConst for that particular probCons
            // string sqlInit = fmt::format("SELECT {} FROM \"{}\" LIMIT 1",attrCon->attr,DB_optim);
            // SingleRow SrIn = SingleRow(sqlInit);
            cout << "Attribute: " << attrCon->attr << " " << attrCon->attrType << endl;
            cout << "Expression: " << spq->getValue(probCon->p) << " " << probCon->psign << " " << spq->getValue(probCon->v) << " " << probCon->vsign << endl;

            // select only the rows where x is != 0, also select the array
            string selectcols = fmt::format("{},{}", "id", attrCon->attr);
            string sql;
            try
            {
                sql = fmt::format("SELECT {} FROM \"{}\" WHERE id = ANY(ARRAY[{}])", selectcols, DB_optim, fmt::join(selectIds, ","));
            }
            catch (...)
            {
                cout << "Query Exception Happened" << endl;
            }
            SingleRow sr = SingleRow(sql);

            while (sr.fetchRow())
            {
                std::vector<double> realization;
                sr.getArray(1, realization);
                int row = sr.getBigInt(0) - 1;
                numScenarios = realization.size();
                if (!initialized)
                {
                    cout << "Initializing InnerConst" << endl;
                    initializeVector(innerConst, numScenarios, 0.0);
                    initialized = true;
                }
                for (int j = 0; j < numScenarios; j++)
                {
                    double value = (double)x[row] * realization[j];
                    innerConst[j] += value;
                }
            }
            // push innerConst to innerConstraints
            this->innerConstraints.push_back(innerConst);
            int satisfied = countSatisfied(numScenarios, innerConst, probCon, spq);
            double rk = calculateRk(probCon, satisfied, numScenarios, spq);
            this->r.push_back(rk);
        }
    }
    this->W_q = calculateObj(x, selectIds, spq->obj, DB_optim);
    cout << "Objective Value = " << W_q << endl;
    cout << "Validation Completed" << endl;
}

// get the vector pairs with realizations and scores
// shuffle them and then find what is partitionSize
// start from beginning, count partitionSize numbers and add them into a partition
// once >= partitionSize reset cnt and start new partition add the prev into partitions
void partition(int M, int Z, std::vector<pair<double, double>> &realizationScore, std::vector<int> &shuffler, std::vector<std::vector<pair<double, double>>> &partitions)
{
    // cout << "Partitioning" << endl;
    int totalLength = realizationScore.size();
    double partLength = ceil(totalLength / (double)Z);
    int partitionSize = static_cast<int>(partLength);
    // create a partition vector of pairs that contains the realization and scores
    std::vector<pair<double, double>> partition;
    int cnt = 0;
    for (int i = 0; i < totalLength; i++)
    {
        // use the shuffler array which is already shuffled in order to retrieve partitionSize number of RANDOM elements into the partition
        partition.push_back(realizationScore[shuffler[i]]);
        cnt++;
        if (cnt >= partitionSize)
        {
            //  push partition to partition
            //  clear the current partition to calculate the next one
            //  reset counter to 0 to count for the next partition
            cnt = 0;
            partitions.push_back(partition);
            partition.clear();
        }
    }
}

template <typename T1, typename T2>
void printVectPair(std::vector<pair<T1, T2>> &v)
{
    for (int i = 0; i < v.size(); i++)
    {
        cout << "(" << v[i].first << ", " << v[i].second << ") ";
    }
    cout << endl;
}

// probConst contains information such as vsign, psign that helps determine min/max summary and sorting order asc/desc
// attrCon helps with SQL query
// conOrder helps to retreive the values of the innerConstraints for each constraint/scenario
// Idea: we have already computed the innerConstraints --> we can form vector<pair<scenario value, score=innerconst>> helps for sorting based on score
// then we do the random partitions
// we go partition by partition and sort each of them depending on the vsign?
// then find summary for each row
// push into a vector
// return the summary-vector
std::vector<std::vector<double>> SummarySearch::summarize(std::vector<int> &x, int M, int Z, double alpha, string DB_optim, shared_ptr<ProbConstraint> probCon, shared_ptr<AttrConstraint> attrCon, int conOrder, vector<int> &reducedIds, bool reduced)
{
    cout << "summarizing" << endl;
    string selectcols = fmt::format("{},{}", "id", attrCon->attr);
    // I have put limit 1 to test easier REMOVE LIMIT 1 for future use
    string sql;
    if (reduced)
    {
        sql = fmt::format("SELECT {} FROM {} WHERE id = ANY(ARRAY[{}]) ORDER BY id", selectcols, DB_optim, fmt::join(reducedIds, ","));
    }
    else
    {
        sql = fmt::format("SELECT {} FROM \"{}\"", selectcols, DB_optim);
    }
    SingleRow sr = SingleRow(sql);
    // we have a vector of doubles, each representing the summary for each partition
    std::vector<std::vector<double>> summaries;
    std::vector<double> summary;
    bool init;

    while (sr.fetchRow())
    {
        std::vector<double> realization;
        std::vector<pair<double, double>> realizationScore;
        sr.getArray(1, realization);
        int id = sr.getBigInt(0) - 1;

        int numScenarios = realization.size();
        // make pairs of realization of tuple read by fetchRow(), and the corresponding score of x for scenario i
        for (int i = 0; i < numScenarios; i++)
        {
            realizationScore.push_back(make_pair(realization[i], this->innerConstraints[conOrder][i]));
        }
        std::vector<std::vector<pair<double, double>>> partitions;
        partition(M, Z, realizationScore, this->shuffler, partitions);

        // summaries is initialized with a number of partitions, summary vectors
        if (!init)
        {
            if (!reduced)
            {
                initializeVector(summary, NTuples, 0.0);
            }
            initializeVector(summaries, partitions.size(), summary);
            init = true;
        }
        //  partition by partition, sort
        //  summarize
        for (int i = 0; i < partitions.size(); i++)
        {
            double summary;
            // decide when ASC when DESC
            if (probCon->vsign == Inequality::gteq) //>= min summary
            {
                sort(partitions[i].begin(), partitions[i].end(), sortbysecDESC);
                // deb(partitions[i]);
                summary = calculateSummary(partitions[i], false, alpha);
            }
            else
            {
                sort(partitions[i].begin(), partitions[i].end(), sortbysecASC);
                summary = calculateSummary(partitions[i], true, alpha);
            }
            if (!reduced)
            {
                summaries[i][id] = summary;
            }
            else
            {
                summaries[i].push_back(summary);
            }
        }
    }
    return summaries;
}

void SummarySearch::reshuffleShuffler(std::vector<int> &shuffle)
{
    cout << "Shuffling" << endl;
    srand(unsigned(time(0)));
    random_shuffle(shuffle.begin(), shuffle.end());
}
/*fetch the attr from the DB and add the corresponding value to the coeff[i]
then formulate the linearExpression and based on the objSense add Maximize or Minimize
*/
void SummarySearch::formSumObj(GRBModel &model, shared_ptr<Objective> obj, GRBVar *xx, vector<int> &reducedIds, bool reduced)
{
    shared_ptr<AttrObjective> attrObj;
    bool isDet = isDeterministic(obj, attrObj);
    if (!isDet || attrObj->objType != numeric_type)
    {
        return;
    }

    cout << "Sum Objective" << endl;

    string selectcols = fmt::format("{},{}", "id", attrObj->obj);
    // string selectcols = attrObj->obj;
    string sql;
    if (reduced)
    {
        string sql = fmt::format("SELECT {} FROM {} WHERE id = ANY(ARRAY[{}]) ORDER BY id", selectcols, DB_optim, fmt::join(reducedIds, ","));
    }
    else
    {
        string sql = fmt::format("SELECT {} FROM {}", selectcols, DB_optim);
    }

    SingleRow sr = SingleRow(sql);

    GRBLinExpr sumObjExpr;
    while (sr.fetchRow())
    {
        int id = sr.getBigInt(0) - 1;
        double coeffVal = sr.getNumeric(1);
        sumObjExpr += coeffVal * xx[id];
        // sumObjExpr.addTerm(coeffVal, xx[id]);
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
void SummarySearch::formExpSumObj(GRBModel &model, shared_ptr<Objective> obj, GRBVar *xx, vector<int> &reducedIds, bool reduced)
{
    shared_ptr<AttrObjective> attrObj;
    bool isDet = isDeterministic(obj, attrObj);
    if (!isDet || attrObj->objType != array_type)
    {
        return;
    }
    cout << "Expected Sum Objective" << endl;

    string table = fmt::format("{}_{}", DB_optim, "summary");
    string columnObj = fmt::format("{}_{}", attrObj->obj, "mean");
    string selectcols = fmt::format("{},{}", "id", columnObj);

    string sql;
    if (reduced)
    {
        sql = fmt::format("SELECT {} FROM {} WHERE id = ANY(ARRAY[{}]) ORDER BY id", selectcols, table, fmt::join(reducedIds, ","));
    }
    else
    {
        sql = fmt::format("SELECT {} FROM {}", selectcols, table);
    }
    SingleRow sr = SingleRow(sql);
    GRBLinExpr expSumObjExpr;

    while (sr.fetchRow())
    {
        int id = sr.getBigInt(0) - 1;
        double coeffVal = sr.getNumeric(1);
        expSumObjExpr += coeffVal * xx[id];
        // expSumObjExpr.addTerm(coeffVal,xx[id]);
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
void SummarySearch::formCntObj(GRBModel &model, shared_ptr<Objective> obj, GRBVar *xx, vector<int> &reducedIds, bool reduced)
{
    shared_ptr<CountObjective> cntObj = getCount(obj);
    if (!cntObj)
    {
        return;
    }

    cout << "Count Objective" << endl;

    GRBLinExpr cntObjExpr;
    if (reduced)
    {
        double coeffval = 1;
        for (int i = 0; i < reducedIds.size(); i++)
        {
            int id = reducedIds[i] - 1;
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

/*
Count Constraint is basically Sum (xi) with bounds
*/
void SummarySearch::formCountCons(GRBModel &model, shared_ptr<Constraint> cons, GRBVar *xx, vector<int> &reducedIds, bool reduced, map<string, Option> &options)
{
    shared_ptr<CountConstraint> cntCons = getCount(cons);
    double ub_eps = 1e30;
    double lb_eps = -1e30;

    bool omitCnt = boost::get<bool>(options.at("omit count constraint"));
    bool scale = boost::get<bool>(options.at("scale down"));
    double scaleFactor = boost::get<double>(options.at("scale factor"));

    if (!cntCons || omitCnt)
    {
        return;
    }

    GRBLinExpr cntConsExpr;
    if (reduced)
    {
        double coeffval = 1;
        for (int i = 0; i < reducedIds.size(); i++)
        {
            int id = reducedIds[i] - 1;
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
    if (scale)
    {
        lb = lb * scaleFactor;
        ub = ub * scaleFactor;
    }

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
}
// select the attr from the DB, and then add the values of each tuple to the corresponding coeff[i]
// formulate GRBLinExpr with addTerms.(coeff,xx,NTuples)
// based on the bounds formulate the constraints ignore the bounds that are not set
void SummarySearch::formSumCons(GRBModel &model, shared_ptr<Constraint> cons, GRBVar *xx, vector<int> &reducedIds, bool reduced)
{
    shared_ptr<BoundConstraint> boundCon = dynamic_pointer_cast<BoundConstraint>(cons);
    shared_ptr<AttrConstraint> attrCon = dynamic_pointer_cast<AttrConstraint>(cons);

    double ub_eps = 1e30;
    double lb_eps = -1e30;

    if (!boundCon || !attrCon || attrCon->attrType != numeric_type)
    {
        return;
    }
    cout << "compare attrCon->attrType  =  " << (attrCon->attrType) << endl;

    string selectcols = fmt::format("{},{}", "id", attrCon->attr);
    // string selectcols = attrCon->attr;
    string sql;
    if (reduced)
    {
        sql = fmt::format("SELECT {} FROM {} WHERE id = ANY(ARRAY[{}]) ORDER BY id", selectcols, DB_optim, fmt::join(reducedIds, ","));
    }
    else
    {
        sql = fmt::format("SELECT {} FROM \"{}\"", selectcols, DB_optim);
    }

    SingleRow sr = SingleRow(sql);
    GRBLinExpr sumConsExpr;

    while (sr.fetchRow())
    {
        int id = sr.getBigInt(0) - 1;
        double coeffVal = sr.getNumeric(1);
        sumConsExpr += coeffVal * xx[id];
        // sumConsExpr.addTerm(coeffVal, xx[id]);
    }

    double lb = spq->getValue(boundCon->lb);
    double ub = spq->getValue(boundCon->ub);
    cout << lb << " " << ub << endl;
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
    cout << "formulated sum constraint" << endl;
}

void clearConstr(GRBModel &model)
{
    cout << "deleting constraints" << endl;
    GRBConstr *constraints = model.getConstrs();

    for (int i = 0; i < model.get(GRB_IntAttr_NumConstrs); i++)
    {
        model.remove(constraints[i]);
    }

    GRBGenConstr *genConstraints = model.getGenConstrs();
    for (int i = 0; i < model.get(GRB_IntAttr_NumGenConstrs); i++)
    {
        model.remove(genConstraints[i]);
    }
    model.update();
}

// go through the constraints and delete them
void removeProbConstr(GRBModel &model, std::vector<GRBVar> &yy, std::vector<GRBGenConstr> &genCon, std::vector<GRBConstr> &sumyCon)
{
    for (int i = 0; i < genCon.size(); i++)
    {
        // cout << "deleting Indicator" << endl;
        model.remove(genCon[i]);
    }
    model.update();
    for (int i = 0; i < sumyCon.size(); i++)
    {
        // cout << "deleting SUM yy[i]" << endl;
        model.remove(sumyCon[i]);
    }
    model.update();
    for (int i = 0; i < yy.size(); i++)
    {
        // cout << "deleting a yy[i]" << endl;
        model.remove(yy[i]);
    }
    model.update();
    yy.clear();
    genCon.clear();
    sumyCon.clear();
}
// check if constraint is probabilistic
// find the number of summaries for that constraint
// remove the variables and constraints related to this prob constraint if any
// find pZ for the SumY
// add indicator constraints for the variables Y[z]
void SummarySearch::formProbCons(GRBModel &model, shared_ptr<Constraint> cons, GRBVar *xx,
                                 std::vector<std::vector<std::vector<double>>> &summaries, int &probConOrder,
                                 std::vector<GRBVar> &yy, std::vector<GRBGenConstr> &genCon, std::vector<GRBConstr> &sumyCon,
                                 vector<int> &reducedIds, bool reduced)
{
    shared_ptr<ProbConstraint> probCon;
    shared_ptr<AttrConstraint> attrCon;

    bool isstoch = isStochastic(cons, probCon, attrCon);
    if (!isstoch)
    {
        return;
    }
    cout << "removing prob constr" << endl;
    removeProbConstr(model, yy, genCon, sumyCon);
    model.update();
    cout << "removed prob constr" << endl;
    // calculate the number of yk indicator variables
    int Z = summaries[probConOrder].size();
    double v = spq->getValue(probCon->v);
    double p = spq->getValue(probCon->p);
    double pZ = ceil(p * Z);
    GRBVar y[Z];
    double coeff_Y[Z];
    for (int z = 0; z < Z; z++)
    {
        y[z] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y[" + to_string(z) + "]");
        yy.push_back(y[z]);
        coeff_Y[z] = 1;
    }
    model.update();
    GRBLinExpr sumYz;
    std::vector<std::vector<double>> S = summaries[probConOrder];
    cout << probConOrder << " " << summaries.size() << endl;
    for (int z = 0; z < S.size(); z++)
    {
        GRBLinExpr innerCons;
        int sz = S[z].size();
        for (int i = 0; i < sz; i++)
        {
            if (reduced)
            {
                double coeffVal = S[z][i];
                int id = reducedIds[i] - 1;
                innerCons += coeffVal * xx[id];
                // innerCons.addTerm(coeffVal,xx[id]);
            }
            else
            {
                double coeffVal = S[z][i];
                int id = i;
                innerCons += coeffVal * xx[id];
                // innerCons.addTerm(coeffVal,xx[id]);
            }
        }

        try
        {
            if (probCon->vsign == Inequality::gteq)
            {
                GRBGenConstr indicator = model.addGenConstrIndicator(y[z], 1, innerCons, GRB_GREATER_EQUAL, v);
                genCon.push_back(indicator);
            }
            else
            {
                GRBGenConstr indicator = model.addGenConstrIndicator(y[z], 1, innerCons, GRB_LESS_EQUAL, v);
                genCon.push_back(indicator);
            }
        }
        catch (GRBException &e)
        {
            cout << "Error code 8 = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        }
    }

    sumYz.addTerms(coeff_Y, y, Z);
    try
    {
        if (probCon->psign == Inequality::gteq)
        {
            GRBConstr constr = model.addConstr(sumYz, GRB_GREATER_EQUAL, pZ);
            sumyCon.push_back(constr);
        }
        else
        {
            GRBConstr constr = model.addConstr(sumYz, GRB_LESS_EQUAL, pZ);
            sumyCon.push_back(constr);
        }
    }
    catch (GRBException &e)
    {
        cout << "Error code 9 = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    probConOrder += 1;
}

int getQuantileIdx(double p, int n)
{
    int idx = n * p;
    return idx;
}

double getExpectedValueQtile(vector<double> &v, int qtileIdx)
{
    double qtile = v[qtileIdx];
    double sum = 0.0;
    int cnt = 0;
    for (int i = 0; i < v.size(); i++)
    {
        if (v[i] <= qtile)
        {
            sum += v[i];
            cnt += 1;
        }
    }

    double mean = sum / cnt;
    return mean;
}

// figure out which cvar corresponds to each var
void SummarySearch::formCVaR(GRBModel &model, shared_ptr<Constraint> cons, GRBVar *xx)
{
    shared_ptr<ProbConstraint> probCon;
    shared_ptr<AttrConstraint> attrCon;

    bool isstoch = isStochastic(cons, probCon, attrCon);
    if (!isstoch)
    {
        return;
    }
    cout << "Formulating CVAR" << endl;
    double v = spq->getValue(probCon->v);
    double p = spq->getValue(probCon->p);

    string selectcols = fmt::format("{},{}_{}", "id", attrCon->attr, "quantiles");
    string table = fmt::format("{}_{}", DB_optim, "summary");
    cout << "here" << endl;
    int qtileNumber = pg.getColumnLength(table, fmt::format("{}_{}", attrCon->attr, "quantiles"));
    cout << "qtileNumber" << qtileNumber << endl;
    int qtileIdx = getQuantileIdx(1 - p, qtileNumber);

    cout << "QTILEIDX = " << qtileIdx << endl;

    // string sql = fmt::format(
    //     "SELECT id, "
    //     "AVG(profit_quantiles[1:{}]) AS avg_quantiles "
    //     "FROM \"{}\"",
    //     qtileIdx, table);

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
        double coeff[NTuples];
        while (sr.fetchRow())
        {
            int id = sr.getBigInt(0) - 1;
            double coeffVal = sr.getNumeric(1);
            coeff[id] = coeffVal;
        }

        GRBLinExpr CVar_m;
        CVar_m.addTerms(coeff, xx, NTuples);
        model.addConstr(CVar_m, GRB_GREATER_EQUAL, v);
    }
    else
    {
        cout << "Currently There's no Implementation for this combination psign and vsign" << endl;
    }
    cout << "formulated CVaR" << endl;
}

/*
check if constraint has bounds, an attribute and is of array_type -> expCons
we can make use of the percomputed table for tuple-wise mean -> determine the table by concatenating "summary"
select cols = concatenate "mean"
assign coeff[NTuples] the tuple-wise means, and then formulate the linearExpr with the coeff and the xx
add constraint based on the bounds
*/
void SummarySearch::formExpCons(GRBModel &model, shared_ptr<Constraint> cons, GRBVar *xx, vector<int> &reducedIds, bool reduced)
{
    cout << "Formulating Exp Constraint" << endl;
    shared_ptr<BoundConstraint> boundCon = dynamic_pointer_cast<BoundConstraint>(cons);
    shared_ptr<AttrConstraint> attrCon = dynamic_pointer_cast<AttrConstraint>(cons);
    if (!boundCon || !attrCon || attrCon->attrType != array_type)
    {
        return;
    }
    cout << "Expectation sum constraint" << endl;
    double ub_eps = 1e30;
    double lb_eps = -1e30;
    string table = fmt::format("{}_{}", DB_optim, "summary");
    string selectAttr = fmt::format("{}_{}", attrCon->attr, "mean");
    string selectcols = fmt::format("{},{}", "id", selectAttr);

    string sql;
    if (reduced)
    {
        sql = fmt::format("SELECT {} FROM \"{}\" WHERE id = ANY(ARRAY[{}]) ORDER BY id", selectcols, table, fmt::join(reducedIds, ","));
    }
    else
    {
        sql = fmt::format("SELECT {} FROM  \"{}\"", selectcols, table);
    }

    SingleRow sr = SingleRow(sql);
    GRBLinExpr expSumConsExpr;

    while (sr.fetchRow())
    {
        int id = sr.getBigInt(0) - 1;
        double coeffVal = sr.getNumeric(1);
        expSumConsExpr += coeffVal * xx[id];
        // expSumConsExpr.addTerm(coeffVal,xx[id]);
    }

    double lb = spq->getValue(boundCon->lb);
    double ub = spq->getValue(boundCon->ub);
    cout << lb << " " << ub << endl;
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
    cout << "Formulated Expectation Constraint" << endl;
}

// naive
// formulateSAA(summary[M],q=0);
void SummarySearch::formulateSAA(GRBModel &model, std::vector<std::vector<std::vector<double>>> &summaries, int q, GRBVar *xx, vector<int> &reducedIds, bool reduced, map<string, Option> &cntoptions)
{
    cout << "formulating" << endl;
    int numCons = spq->cons.size();
    int probConOrder = 0;
    for (int i = 0; i < numCons; i++)
    {
        if (q == 0)
        {
            cout << "tuka" << endl;
            formCountCons(model, spq->cons[i], xx, reducedIds, reduced, cntoptions);
            formSumCons(model, spq->cons[i], xx, reducedIds, reduced);
            formExpCons(model, spq->cons[i], xx, reducedIds, reduced);
        }
        else
        {
            formProbCons(model, spq->cons[i], xx, summaries, probConOrder, yy, genCon, sumyCon, reducedIds, reduced);
        }
    }
    if (q == 0)
    {
        formSumObj(model, spq->obj, xx, reducedIds, reduced);
        formExpSumObj(model, spq->obj, xx, reducedIds, reduced);
        formCntObj(model, spq->obj, xx, reducedIds, reduced);
    }
    model.update();
}

void SummarySearch::formulateLP(GRBModel &model, GRBVar *xx, bool reduced, map<string, Option> &cntoptions)
{
    cout << "Formulating LP" << endl;
    vector<int> reducedIds;
    int numCons = spq->cons.size();
    int probConOrder = 0;
    for (int i = 0; i < numCons; i++)
    {
        formCountCons(model, spq->cons[i], xx, reducedIds, reduced, cntoptions);
        formSumCons(model, spq->cons[i], xx, reducedIds, reduced);
        formExpCons(model, spq->cons[i], xx, reducedIds, reduced);
        formCVaR(model, spq->cons[i], xx);
    }
    formSumObj(model, spq->obj, xx, reducedIds, reduced);
    formExpSumObj(model, spq->obj, xx, reducedIds, reduced);
    formCntObj(model, spq->obj, xx, reducedIds, reduced);
    model.update();
}

void SummarySearch::guessOptimalConservativeness(std::vector<std::vector<pair<double, double>>> &history, std::vector<double> &alpha)
{
    fitter.fit_and_predict(history, alpha);
}

bool solutionInH(vector<int> &x, vector<double> &alpha, History &H)
{
    for (auto &tuple : H.bestSolMetadata)
    {
        bool xIN = true;
        bool alphaIN = true;

        // traverse x from tuple
        for (int j = 0; j < get<0>(tuple).size(); j++)
        {
            if (x[j] != get<0>(tuple)[j])
            {
                xIN = false;
                break;
            }
        }

        for (int j = 0; j < get<1>(tuple).size(); j++)
        {
            if (alpha[j] != get<1>(tuple)[j])
            {
                alphaIN = false;
                break;
            }
        }

        if (xIN && alphaIN)
        {
            return true;
        }
    }
    return false;
}

bool isFeasible(vector<double> &r)
{
    for (int k = 0; k < r.size(); k++)
    {
        if (r[k] < 0)
        {
            return false;
        }
    }
    return true;
}

double calculateEpsilonQ(shared_ptr<StochasticPackageQuery> spq, double W_q, double W0)
{
    double epsilonQ;
    if (spq->obj->objSense == maximize)
    {
        epsilonQ = W0 / W_q - 1;
    }
    else
    {
        epsilonQ = (W_q / W0) - 1;
    }

    return epsilonQ;
}

void SummarySearch::Best(shared_ptr<Objective> obj, std::vector<int> &x, History &H)
{
    if (!H.foundFeasible)
    {
        cout << "Current Best is NOT Feasible" << endl;
        if (isFeasible(this->r))
        {
            cout << "New Solution is Feasible -> BEST" << endl;
            H.foundFeasible = true;
            H.wBest = W_q;
            H.epsilonBest = calculateEpsilonQ(this->spq, this->W_q, this->W0);
            for (int i = 0; i < x.size(); i++)
            {
                H.xBest.push_back(x[i]);
            }
        }
        return;
    }
    cout << "There's already a feasible solution" << endl;
    if (obj->objSense == minimize)
    {
        if (isFeasible(this->r) && W_q < H.wBest)
        {
            cout << "New Solution has better Objective value MINIMIZE" << endl;
            cout << "OLD BEST = " << H.wBest << " NEW BEST = " << W_q << endl;
            H.xBest.clear();
            H.wBest = W_q;
            H.epsilonBest = calculateEpsilonQ(this->spq, this->W_q, this->W0);
            for (int i = 0; i < x.size(); i++)
            {
                H.xBest.push_back(x[i]);
            }
        }
    }
    else
    {
        if (isFeasible(this->r) && W_q > H.wBest)
        {
            cout << "New Solution has better Objective value MAXIMIZE" << endl;
            cout << "OLD BEST = " << H.wBest << " NEW BEST = " << W_q << endl;
            H.xBest.clear();
            H.wBest = W_q;
            H.epsilonBest = calculateEpsilonQ(this->spq, this->W_q, this->W0);
            for (int i = 0; i < x.size(); i++)
            {
                H.xBest.push_back(x[i]);
            }
        }
    }
}

// solve() -> model.optimize(); then translate the x[i] GRBVar vector into a x vector so you can manipulate it.
template <typename T>
void SummarySearch::solve(GRBModel &model, std::vector<T> &x, GRBVar *xx, vector<int> &reducedIds, bool reduced)
{
    // this->model.reset();
    cout << "Optimizing Model" << endl;
    model.optimize();
    if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
    {
        std::cout << "Initial model is infeasible." << std::endl;
        x.clear();
        return;
    }
    else
    {
        std::cout << "Model optimization ended with status: "
                  << model.get(GRB_IntAttr_Status) << std::endl;
    }
    for (int i = 0; i < this->NTuples; i++)
    {
        x[i] = static_cast<T>(xx[i].get(GRB_DoubleAttr_X));
    }
}

SolutionMetadata SummarySearch::CSASolve(std::vector<int> &x, int M, int Z, vector<int> &reducedIds, bool reduced, map<string, Option> &cntoptions)
{
    int q = 0;
    H.bestSolMetadata.clear();
    H.curveFitMetadata.clear();
    std::vector<double> alpha;
    std::vector<pair<double, double>> alphaK_rK;
    initializeVector(H.curveFitMetadata, probConstCnt, alphaK_rK);
    initializeVector(alpha, probConstCnt, 0.0);
    while (true)
    {
        cout << "Iteration =" << q << " Z = " << Z << endl;
        if (solutionInH(x, alpha, this->H))
        {
            cout << "SOLUTION IS IN THE SET OF SOLUTIONS" << endl;
            SolutionMetadata sol(H.xBest, H.wBest, H.epsilonBest, H.foundFeasible);
            return sol;
        }

        // validate() -> this will set the rk values, and calculate the Wq
        validate(x, this->spq, 1e6, this->DB_optim);
        Best(spq->obj, x, H);

        // update the metadata for finding the best folution and for curve fit
        bool isFeas = isFeasible(r);
        H.bestSolMetadata.emplace_back(x, alpha, isFeas, W_q);

        for (int k = 0; k < probConstCnt; k++)
        {
            H.curveFitMetadata[k].push_back(make_pair(alpha[k], r[k]));
        }

        if (q == 0)
        {
            W0 = W_q;
        }

        // loop through the rk values and check if they are all >= 0 -> if true
        // calculate the e^Q -> if <= epsilon
        double epsilonQ = calculateEpsilonQ(this->spq, this->W_q, this->W0);
        cout << "EPSILON = " << epsilonQ << endl;
        cout << (isFeas ? "Feasible" : "Infeasible") << endl;
        cout << epsilon << " " << (epsilonQ <= this->epsilon ? "Good Bound" : "Bad Bound") << endl;

        //-> return this solution -> this means update xBest to be equal to X, and set all metadata
        if (isFeas && epsilonQ <= this->epsilon)
        {
            SolutionMetadata sol(x, W_q, epsilonQ, true);
            return sol;
        }

        q = q + 1;
        guessOptimalConservativeness(H.curveFitMetadata, alpha);
        deb(H.curveFitMetadata);
        deb(alpha);

        int cnt = 0;
        std::vector<std::vector<std::vector<double>>> summaries;
        for (int i = 0; i < spq->cons.size(); i++)
        {
            shared_ptr<ProbConstraint> probCon;
            shared_ptr<AttrConstraint> attrCon;
            bool isstoch = isStochastic(spq->cons[i], probCon, attrCon);
            if (isstoch)
            {
                reshuffleShuffler(shuffler);
                std::vector<std::vector<double>> summariesCons = summarize(x, M, Z, alpha[cnt], DB_optim, probCon, attrCon, cnt, reducedIds, reduced);
                summaries.push_back(summariesCons);
                cnt++;
            }
        }

        formulateSAA(this->modelILP, summaries, q, this->xxILP.get(), reducedIds, reduced, cntoptions);
        solve(this->modelILP, x, this->xxILP.get(), reducedIds, reduced);
        if (x.empty())
        {
            cout << "FORMULATION IS INFEASIBLE" << endl;
            SolutionMetadata sol(x, 0, 0, false);
            return sol;
        }
    }
}

BinarySearchMetadata SummarySearch::guessOptimalConservativenessBinarySearch(double low, double high, double rk)
{
    double alpha = low + (high - low) / 2;
    BinarySearchMetadata metadata(0.0, 0.0, 0.0);
    if (rk < 0)
    {
        // the solution is infeasible -> use more conservative summary
        cout << "UPPPPPPP" << endl;
        metadata.low = alpha;
        metadata.high = high;
        metadata.alpha = alpha;
    }
    else
    {
        cout << "DOWNNN" << endl;
        // the solution is feasible but suboptimal or the system is infeasible -> use less conservative summary
        metadata.low = low;
        metadata.high = alpha;
        metadata.alpha = alpha;
    }

    return metadata;
}

SolutionMetadata SummarySearch::CSASolveBinSearch(std::vector<int> &x, int M, int Z, vector<int> &reducedIds, bool reduced, map<string, Option> &cntoptions)
{
    cout << "SOLVING WITH BINARY SEARCH FITTER" << endl;
    BinarySearchMetadata alpha_KMetadata(1e-5, 1.0, 0.0);
    std::vector<BinarySearchMetadata> history;
    initializeVector(history, probConstCnt, alpha_KMetadata);
    std::vector<double> alpha;
    initializeVector(alpha, probConstCnt, 0.0);
    int q = 0;
    int qAfterZequalsM = 0;

    SolutionMetadata bestSol(x, 0, 0, false);

    while (true)
    {
        cout << "Iteration: " << q << " Z = " << Z << endl;
        // validate() -> this will set the rk values, and calculate the Wq
        if (x.size() > 0)
        {

            validate(x, this->spq, 1e6, this->DB_optim);

            // we need both isFeas and rk for the binary search

            bool isFeas = isFeasible(r);

            if (q == 0)
            {
                W0 = W_q;
            }

            // loop through rk values and check if they are all >= 0 -> if true
            // calculate the e^Q -> if <= epsilon
            double epsilonQ = calculateEpsilonQ(this->spq, this->W_q, this->W0);
            cout << "EPSILON = " << epsilonQ << endl;
            cout << (isFeas ? "Feasible" : "Infeasible") << endl;
            cout << epsilon << " " << (epsilonQ <= this->epsilon ? "Good Bound" : "Bad Bound") << endl;

            //-> return this solution -> this means update xBest to be equal to X, and set all metadata
            if (isFeas && epsilonQ <= this->epsilon)
            {
                SolutionMetadata sol(x, W_q, epsilonQ, true);
                bestSol.Z = Z;
                bestSol.x = x;
                bestSol.w = W_q;
                bestSol.epsilon = epsilonQ;
                bestSol.isFeasible = true;
                bestSol.isOptimal = true;
                return bestSol;
            }
            else if (isFeas)
            {
                bestSol.Z = Z;
                bestSol.x = x;
                bestSol.w = W_q;
                bestSol.epsilon = epsilonQ;
                bestSol.isFeasible = true;
                bestSol.isOptimal = false;
            }
        }
        else
        {
            SolutionMetadata sol(x, 0, 0, false);
            return sol;
        }

        q = q + 1;
        if (Z == cntScenarios)
        {
            qAfterZequalsM += 1;
        }

        // call the binary search function with for each constraint

        for (int i = 0; i < probConstCnt; i++)
        {
            int ScenariosLow = (int)ceil(history[i].low * (this->cntScenarios / Z));
            int ScenariosRight = (int)ceil(history[i].high * (this->cntScenarios / Z));
            if (ScenariosLow != ScenariosRight || (Z == cntScenarios && qAfterZequalsM <= 1))
            {
                if (x.size() == 0)
                {
                    cout << "Before call of Binary Search INFEASIBLE"<<endl;
                }
                else
                {
                    cout << "Before call of Binary Search: " << history[i].low << " " << history[i].high << " " << r[i] << endl;
                    BinarySearchMetadata metadata = guessOptimalConservativenessBinarySearch(history[i].low, history[i].high, r[i]);
                    history[i].low = metadata.low;
                    history[i].high = metadata.high;
                    history[i].alpha = metadata.alpha;
                    alpha[i] = metadata.alpha;
                    cout << "After call of Binary Search: " << history[i].low << " " << history[i].high << " " << r[i] << endl;
                }
            }
            else
            {
                cout << "BINARY SEARCH END CONDITION MET" << endl;
                return bestSol;
            }
        }

        deb(alpha);

        int cnt = 0;
        std::vector<std::vector<std::vector<double>>> summaries;
        for (int i = 0; i < spq->cons.size(); i++)
        {
            shared_ptr<ProbConstraint> probCon;
            shared_ptr<AttrConstraint> attrCon;
            bool isstoch = isStochastic(spq->cons[i], probCon, attrCon);
            if (isstoch)
            {
                reshuffleShuffler(shuffler);
                std::vector<std::vector<double>> summariesCons = summarize(x, M, Z, alpha[cnt], DB_optim, probCon, attrCon, cnt, reducedIds, reduced);
                summaries.push_back(summariesCons);
                cnt++;
            }
        }
        formulateSAA(this->modelILP, summaries, q, this->xxILP.get(), reducedIds, reduced, cntoptions);
        solve(this->modelILP, x, this->xxILP.get(), reducedIds, reduced);
    }
}

SolutionMetadata SummarySearch::summarySearch(shared_ptr<StochasticPackageQuery> spq, int M_hat, int M, int m, int z, vector<int> &reducedIds, bool reduced, map<string, Option> &cntoptions, map<string, Option> &curveFitOptions)
{
    std::vector<std::vector<std::vector<double>>> summaries;
    // modelILP is global so we need to delete the constraints for each time we run new Summary Search
    clearConstr(this->modelILP);
    yy.clear();
    genCon.clear();
    sumyCon.clear();
    // formulate Deterministic ILP
    formulateSAA(this->modelILP, summaries, 0, this->xxILP.get(), reducedIds, reduced, cntoptions);
    //printConstraints(this->modelILP);
    vector<int> x0;
    initializeVector(x0, NTuples, 0);
    solve(this->modelILP, x0, this->xxILP.get(), reducedIds, reduced);

    bool binSearch = boost::get<bool>(curveFitOptions.at("binarySearch"));
    bool curveFit = boost::get<bool>(curveFitOptions.at("arctan"));

    int Z = 1;
    while (true)
    {
        vector<int> x;
        for (int i = 0; i < x0.size(); i++)
        {
            x.push_back(x0[i]);
        }

        SolutionMetadata sol;
        if (binSearch)
        {
            sol = CSASolveBinSearch(x, M, Z, reducedIds, reduced, cntoptions);
        }
        else
        {
            sol = CSASolve(x, M, Z, reducedIds, reduced, cntoptions);
        }
        if (sol.isFeasible && sol.epsilon <= this->epsilon)
        {
            return sol;
        }
        else
        {
            if (cntScenarios == Z)
            {
                vector<int> xxx;
                SolutionMetadata sol(xxx, 0, 0, false);
                sol.Z = cntScenarios;
                return sol;
            }
            Z = Z + min(z, cntScenarios - Z);
            z = z * 2;
        }
    }
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
    for(int i = 0; i < solDet.size(); i++)
    {
        if(solCVaR[i] > 0 || solDet[i] > 0)
        {
            reducedIds.insert(i + 1);
        }
    }
}

void findUnionDetCVaR2(vector<double> &solDet, vector<double> &solCVaR, vector<double> &solLP1)
{
    for(int i = 0; i < solDet.size(); i++)
    {
        double value = (solDet[i] + solCVaR[i])/2;
        solLP1[i] = value;
    }
}

void populateReducedIds(vector<double> &solLP, set<int> &reducedIds)
{
    for(int i = 0; i < solLP.size(); i++)
    {
        reducedIds.insert(i + 1);
    }
}

void SummarySearch::solveLP2(GRBModel &model, vector<double> &sol, GRBVar *xx, double ub, vector<int> dummyVect)
{
    for (int i = 0; i < NTuples; i++)
    {
        xx[i].set(GRB_DoubleAttr_UB, ub);
    }

    solve(model, sol, xx, dummyVect, false);
}

void SummarySearch::formulateDeterministicLP(GRBModel &model, GRBVar *xx, bool reduced, map<string, Option> &cntoptions)
{
    cout << "Formulating Deterministic LP" << endl;
    vector<int> reducedIds;
    int numCons = spq->cons.size();
    int probConOrder = 0;
    for (int i = 0; i < numCons; i++)
    {
        formCountCons(model, spq->cons[i], xx, reducedIds, reduced, cntoptions);
        formSumCons(model, spq->cons[i], xx, reducedIds, reduced);
    }
    formSumObj(model, spq->obj, xx, reducedIds, reduced);
    formCntObj(model, spq->obj, xx, reducedIds, reduced);
    formExpSumObj(model, spq->obj, xx, reducedIds, reduced);
    model.update();
}

void SummarySearch::formulateCVaRLP(GRBModel &model, GRBVar *xx, bool reduced, map<string, Option> &cntoptions)
{
    cout << "Formulating CVaR LP" << endl;
    vector<int> reducedIds;
    int numCons = spq->cons.size();
    int probConOrder = 0;
    for (int i = 0; i < numCons; i++)
    {
        formCVaR(model, spq->cons[i], xx);
    }
    formSumObj(model, spq->obj, xx, reducedIds, reduced);
    formCntObj(model, spq->obj, xx, reducedIds, reduced);
    formExpSumObj(model, spq->obj, xx, reducedIds, reduced);
    model.update();
}

// '''The idea here is that if the LP1 system is infeasible
//      --> formulate separate LPDet1 and LPCVaR1, so that we get tuples that will satisfy both constraints
//          i.e. narrow the number of tuples down only to tuples that are meaningful for the problem itself'''
SolutionMetadata SummarySearch::stochDualReducer(shared_ptr<StochasticPackageQuery> spq, int qSz, map<string, Option> &cntoptions, map<string, Option> &curveFitOptions)
{
    bool stochLP1feas = true;

    // GRBModel &model, GRBVar *xx, bool reduced
    xxLP1 = std::make_unique<GRBVar[]>(NTuples);
    xxDetLP1 = std::make_unique<GRBVar[]>(NTuples);
    xxCVaRLP1 = std::make_unique<GRBVar[]>(NTuples);

    for (int i = 0; i < NTuples; i++)
    {
        xxLP1[i] = modelLP1.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "xxLP1[" + to_string(i) + "]");
        xxDetLP1[i] = modelDetLP1.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "xxDetLP1[" + to_string(i) + "]");
        xxCVaRLP1[i] = modelCVaRLP1.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "xxCVaRLP1[" + to_string(i) + "]");
    }

    vector<int> dummyVect; // empty vector passing to solve
    vector<double> solLP1;
    vector<double> solDetLP1;
    vector<double> solCVaRLP1;

    initializeVector(solLP1, NTuples, 0.0);
    formulateLP(modelLP1, xxLP1.get(), false, cntoptions);
    solve(modelLP1, solLP1, xxLP1.get(), dummyVect, false);

    if (solLP1.size() == 0)
    {
        stochLP1feas = false;
        initializeVector(solLP1, NTuples, 0.0);
        initializeVector(solDetLP1, NTuples, 0.0);
        initializeVector(solCVaRLP1, NTuples, 0.0);
        // formulate and solve the deterministic and CVaR LP
        formulateDeterministicLP(modelDetLP1, xxDetLP1.get(), false, cntoptions);
        formulateCVaRLP(modelCVaRLP1, xxCVaRLP1.get(), false, cntoptions);
        solve(modelDetLP1, solDetLP1, xxDetLP1.get(), dummyVect, false);
        solve(modelCVaRLP1, solCVaRLP1, xxCVaRLP1.get(), dummyVect, false);
        findUnionDetCVaR2(solDetLP1, solCVaRLP1, solLP1);
    }
    double E = calculateE(solLP1);

    while (true)
    {
        double ub = E / qSz;
        vector<double> solLP2;
        initializeVector(solLP2, NTuples, 0.0);
        set<int> reducedIds2;
        if (ub > 0)
        {
            if (stochLP1feas)
            {
                solveLP2(modelLP1, solLP2, xxLP1.get(), ub, dummyVect);
            }
            else
            {
                vector<double> solDetLP2;
                vector<double> solCVaRLP2;
                initializeVector(solDetLP2, NTuples, 0.0);
                initializeVector(solCVaRLP2, NTuples, 0.0);
                // find the union of both solutions and store it in solLP1 in solve2 we ubdate the upper bound
                solveLP2(modelDetLP1, solDetLP2, xxDetLP1.get(), ub, dummyVect);
                solveLP2(modelCVaRLP1, solCVaRLP2, xxCVaRLP1.get(), ub, dummyVect);
                findUnionDetCVaR2(solDetLP2, solCVaRLP2, solLP2);
            }
        }
        else
        {
            // if the solution is empty even on deterministic LP then it's infeasible
            vector<int> v;
            SolutionMetadata sol(v, -1.0, -1.0, 0);
            return sol;
        }

        vector<int> reducedIds;
        findUnion(solLP1, solLP2, reducedIds);

        int steps = 0;
        if (reducedIds.size() < qSz)
        {
            double low = ub;
            double high = 1.0;
            double eps = 1e-6;
            while (high - low > eps)
            {
                steps += 1;
                reducedIds.clear();
                reducedIds2.clear();
                double mid = low + (high - low) / 2;
                cout << "Low = " << low << endl;
                cout << "High = " << high << endl;
                cout << "Mid = " << mid << endl;
                if (stochLP1feas)
                {
                    //solveLP2 updates the upper bound as well
                    solveLP2(modelLP1, solLP2, xxLP1.get(), mid, dummyVect);
                }
                else
                {
                    vector<double> solDetLP2;
                    vector<double> solCVaRLP2;
                    initializeVector(solDetLP2, NTuples, 0.0);
                    initializeVector(solCVaRLP2, NTuples, 0.0);
                    // find the union of both solutions and store it in solLP1 in solve2 we ubdate the upper bound
                    solveLP2(modelDetLP1, solDetLP2, xxDetLP1.get(), ub, dummyVect);
                    solveLP2(modelCVaRLP1, solCVaRLP2, xxCVaRLP1.get(), ub, dummyVect);
                    findUnionDetCVaR2(solDetLP2, solCVaRLP2, solLP2);
                }

                findUnion(solLP1, solLP2, reducedIds);

                if (reducedIds.size() >= qSz)
                {
                    break;
                }

                if (reducedIds.size() < qSz && solLP2.size() > 0)
                {
                    high = mid;
                }

                if (solLP2.size() == 0)
                {
                    low = mid;
                }
            }
        }

        cntoptions["omit count constraint"] = false;
        cntoptions["scale down"] = false;

        SolutionMetadata sol = summarySearch(this->spq, M_hat, M, 5, 1, reducedIds, true, cntoptions, curveFitOptions);
        if (sol.isFeasible || qSz == NTuples)
        {
            sol.binarySearchSteps = steps;
            sol.qSz = qSz;
            return sol;
        }
        else
        {
            qSz = min(2 * qSz, NTuples);
        }
    }
}



// SolutionMetadata SummarySearch::stochDualReducer(shared_ptr<StochasticPackageQuery> spq, int qSz, map<string, Option> &cntoptions, map<string, Option> &curveFitOptions)
// {
//     bool stochLP1feas = true;

//     // GRBModel &model, GRBVar *xx, bool reduced
//     xxLP1 = std::make_unique<GRBVar[]>(NTuples);
//     xxDetLP1 = std::make_unique<GRBVar[]>(NTuples);
//     xxCVaRLP1 = std::make_unique<GRBVar[]>(NTuples);

//     for (int i = 0; i < NTuples; i++)
//     {
//         xxLP1[i] = modelLP1.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "xxLP1[" + to_string(i) + "]");
//         xxDetLP1[i] = modelDetLP1.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "xxDetLP1[" + to_string(i) + "]");
//         xxCVaRLP1[i] = modelCVaRLP1.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "xxCVaRLP1[" + to_string(i) + "]");
//     }

//     vector<int> dummyVect; // empty vector passing to solve
//     vector<double> solLP1;
//     vector<double> solDetLP1;
//     vector<double> solCVaRLP1;
//     set<int> reducedIds1;

//     initializeVector(solLP1, NTuples, 0.0);
//     formulateLP(modelLP1, xxLP1.get(), false, cntoptions);
//     solve(modelLP1, solLP1, xxLP1.get(), dummyVect, false);
//     populateReducedIds(solLP1, reducedIds1);


//     double E = calculateE(solLP1);


//     if (solLP1.size() == 0)
//     {
//         stochLP1feas = false;
//         initializeVector(solLP1, NTuples, 0.0);
//         initializeVector(solDetLP1, NTuples, 0.0);
//         initializeVector(solCVaRLP1, NTuples, 0.0);
//         // formulate and solve the deterministic and CVaR LP
//         formulateDeterministicLP(modelDetLP1, xxDetLP1.get(), false, cntoptions);
//         formulateCVaRLP(modelCVaRLP1, xxCVaRLP1.get(), false, cntoptions);
//         solve(modelDetLP1, solDetLP1, xxDetLP1.get(), dummyVect, false);
//         solve(modelCVaRLP1, solCVaRLP1, xxCVaRLP1.get(), dummyVect, false);
//         double E1 = calculateE(solDetLP1);
//         double E2 = calculateE(solCVaRLP1);
//         // find the union of both solutions and store it in solLP1
//         //findUnion between solDetLP1, solCVaRLP1 --> for each index if either lp has a positive value store it in reducedIds vector
//         findUnionDetCVaR(solDetLP1, solCVaRLP1, reducedIds1);
//         E = (E1 + E2) / 2;
//         cout << "Size = " << E << endl;
//     }

//     while (true)
//     {
//         double ub = E / qSz;
//         vector<double> solLP2;
//         initializeVector(solLP2, NTuples, 0.0);
//         set<int> reducedIds2;
//         if (ub > 0)
//         {
//             if (stochLP1feas)
//             {
//                 solveLP2(modelLP1, solLP2, xxLP1.get(), ub, dummyVect);
//                 populateReducedIds(solLP2, reducedIds2);
//             }
//             else
//             {
//                 vector<double> solDetLP2;
//                 vector<double> solCVaRLP2;
//                 initializeVector(solDetLP2, NTuples, 0.0);
//                 initializeVector(solCVaRLP2, NTuples, 0.0);
//                 // find the union of both solutions and store it in solLP1 in solve2 we ubdate the upper bound
//                 solveLP2(modelDetLP1, solDetLP2, xxDetLP1.get(), ub, dummyVect);
//                 solveLP2(modelCVaRLP1, solCVaRLP2, xxCVaRLP1.get(), ub, dummyVect);
//                 findUnionDetCVaR(solDetLP2, solCVaRLP2, reducedIds2);
//             }
//         }
//         else
//         {
//             // if the solution is empty even on deterministic LP then it's infeasible
//             vector<int> v;
//             SolutionMetadata sol(v, -1.0, -1.0, 0);
//             return sol;
//         }

//         vector<int> reducedIds;
//         std::set_union(reducedIds2.begin(), reducedIds2.end(),reducedIds1.begin(), reducedIds1.end(),std::back_inserter(reducedIds));

//         int steps = 0;
//         if (reducedIds.size() < qSz)
//         {
//             double low = ub;
//             double high = 1.0;
//             double eps = 1e-6;
//             while (high - low > eps)
//             {
//                 steps += 1;
//                 reducedIds.clear();
//                 reducedIds2.clear();
//                 double mid = low + (high - low) / 2;
//                 cout << "Low = " << low << endl;
//                 cout << "High = " << high << endl;
//                 cout << "Mid = " << mid << endl;
//                 if (stochLP1feas)
//                 {
//                     //solveLP2 updates the upper bound as well
//                     solveLP2(modelLP1, solLP2, xxLP1.get(), mid, dummyVect);
//                     populateReducedIds(solLP2, reducedIds2);
//                 }
//                 else
//                 {
//                     vector<double> solDetLP2;
//                     vector<double> solCVaRLP2;
//                     initializeVector(solDetLP2, NTuples, 0.0);
//                     initializeVector(solCVaRLP2, NTuples, 0.0);
//                     // find the union of both solutions and store it in solLP1 in solve2 we ubdate the upper bound
//                     solveLP2(modelDetLP1, solDetLP2, xxDetLP1.get(), ub, dummyVect);
//                     solveLP2(modelCVaRLP1, solCVaRLP2, xxCVaRLP1.get(), ub, dummyVect);
//                     findUnionDetCVaR(solDetLP2, solCVaRLP2, reducedIds2);
//                 }
//                 std::set_union(reducedIds2.begin(), reducedIds2.end(),reducedIds1.begin(), reducedIds1.end(),std::back_inserter(reducedIds));
//                 cout << "Count of positive indices = " << reducedIds.size() << endl;
//                 if (reducedIds.size() >= qSz)
//                 {
//                     break;
//                 }

//                 if (reducedIds.size() < qSz && solLP2.size() > 0)
//                 {
//                     high = mid;
//                 }

//                 if (solLP2.size() == 0)
//                 {
//                     low = mid;
//                 }
//             }
//         }

//         cntoptions["omit count constraint"] = false;
//         cntoptions["scale down"] = false;

//         SolutionMetadata sol = summarySearch(this->spq, M_hat, M, 5, 1, reducedIds, true, cntoptions, curveFitOptions);
//         if (sol.isFeasible || qSz == NTuples)
//         {
//             sol.binarySearchSteps = steps;
//             sol.qSz = qSz;
//             return sol;
//         }
//         else
//         {
//             qSz = min(2 * qSz, NTuples);
//         }
//     }
// }



// The commented implementation contains an implementation that only solves a deterministic LP, if the LP1 is system infeasible
//  SolutionMetadata SummarySearch::stochDualReducer(shared_ptr<StochasticPackageQuery> spq, int qSz, map<string, Option> &cntoptions, map<string, Option> &curveFitOptions)
//  {
//      bool stochLP1feas = true;

//     // GRBModel &model, GRBVar *xx, bool reduced
//     xxLP1 = std::make_unique<GRBVar[]>(NTuples);
//     xxDetLP1 = std::make_unique<GRBVar[]>(NTuples);

//     for (int i = 0; i < NTuples; i++)
//     {
//         xxLP1[i] = modelLP1.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "xxLP1[" + to_string(i) + "]");
//         xxDetLP1[i] = modelDetLP1.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "xxDetLP1[" + to_string(i) + "]");
//     }

//     vector<int> dummyVect; // empty vector passing to solve
//     vector<double> solLP1;
//     initializeVector(solLP1, NTuples, 0.0);
//     formulateLP(modelLP1, xxLP1.get(), false, cntoptions);
//     solve(modelLP1, solLP1, xxLP1.get(), dummyVect, false);

//     double E = calculateE(solLP1);

//     deb(solLP1);

//     if (solLP1.size() == 0)
//     {
//         stochLP1feas = false;
//         initializeVector(solLP1, NTuples, 0.0);
//         formulateDeterministicLP(modelDetLP1, xxDetLP1.get(), false, cntoptions);
//         solve(modelDetLP1, solLP1, xxDetLP1.get(), dummyVect, false);
//         E = calculateE(solLP1);
//         deb(solLP1);
//         cout<<"Size = "<<E<<endl;
//     }

//     while (true)
//     {
//         double ub = E / qSz;
//         vector<double> solLP2;
//         initializeVector(solLP2, NTuples, 0.0);
//         if (ub > 0)
//         {
//             if(stochLP1feas)
//             {
//                 solveLP2(modelLP1, solLP2, xxLP1.get(), ub, dummyVect);
//             }else
//             {
//                 solveLP2(modelDetLP1, solLP2, xxDetLP1.get(), ub, dummyVect);
//             }
//         }
//         else
//         {
//             // if the solution is empty even on deterministic LP then it's infeasible
//             vector<int> v;
//             SolutionMetadata sol(v, -1.0, -1.0, 0);
//             return sol;
//         }

//         vector<int> reducedIds;
//         findUnion(solLP1, solLP2, reducedIds);

//         int steps = 0;
//         if (reducedIds.size() < qSz)
//         {
//             double low = ub;
//             double high = 1.0;
//             double eps = 1e-6;
//             while (high - low > eps)
//             {
//                 steps += 1;
//                 reducedIds.clear();
//                 double mid = low + (high - low) / 2;
//                 cout << "Low = "<< low <<endl;
//                 cout<< "High = " << high << endl;
//                 cout << "Mid = " << mid << endl;
//                 if(stochLP1feas)
//                 {
//                     solveLP2(modelLP1, solLP2, xxLP1.get(), mid, dummyVect);
//                 }else
//                 {
//                     solveLP2(modelDetLP1, solLP2, xxDetLP1.get(), mid, dummyVect);
//                 }
//                 findUnion(solLP1, solLP2, reducedIds);
//                 cout << "Count orf positive indices = " << reducedIds.size() << endl;
//                 if (reducedIds.size() >= qSz)
//                 {
//                     break;
//                 }

//                 if (reducedIds.size() < qSz && solLP2.size() > 0)
//                 {
//                     high = mid;
//                 }

//                 if (solLP2.size() == 0)
//                 {
//                     low = mid;
//                 }
//             }
//         }

//         deb(reducedIds);
//         cout << "Number of Reduced Ids = " << reducedIds.size() << endl;

//         cntoptions["omit count constraint"] = false;
//         cntoptions["scale down"] = false;

//         SolutionMetadata sol = summarySearch(this->spq, M_hat, M, 5, 1, reducedIds, true, cntoptions, curveFitOptions);
//         if (sol.isFeasible || qSz == NTuples)
//         {
//             sol.binarySearchSteps = steps;
//             sol.qSz = qSz;
//             return sol;
//         }
//         else
//         {
//             qSz = min(2 * qSz, NTuples);
//         }
//     }
// }

void SummarySearch::populateShuffler(vector<int> &v)
{
    for (int i = 0; i < this->spq->cons.size(); i++)
    {
        shared_ptr<ProbConstraint> probCon;
        shared_ptr<AttrConstraint> attrCon;
        bool isstoch = isStochastic(spq->cons[i], probCon, attrCon);
        if (isstoch)
        {
            string selectcols = attrCon->attr;
            string sql = fmt::format("SELECT {} FROM \"{}\" LIMIT 1", selectcols, DB_optim);
            SingleRow sr = SingleRow(sql);
            int row;
            int numScenarios = 0;
            while (sr.fetchRow())
            {
                vector<double> realization;
                sr.getArray(0, realization);
                numScenarios = realization.size();
                cout << "POPULATING SHUFFLER WITH " << numScenarios << " ELEMENTS" << endl;
                for (int j = 0; j < numScenarios; j++)
                {
                    v.push_back(j);
                }
            }
            break;
        }
    }
}

SummarySearch::SummarySearch(int M, int M_hat, shared_ptr<StochasticPackageQuery> spq, double epsilon)
{
    this->M = M;
    this->M_hat = M_hat;
    this->spq = spq;
    // get NTuples from spq->tableName and populate X with NTuples zeroes
    this->DB_optim = spq->tableName;
    // test
    this->DB_valid = fmt::format("{}_{}", DB_optim, "validate");
    this->NTuples = pg.getTableSize(spq->tableName);
    this->cntScenarios = pg.getColumnLength(spq->tableName, "profit");
    this->probConstCnt = countProbConst(spq);
    for (int i = 0; i < probConstCnt; i++)
    {
        std::vector<pair<double, double>> alphaK_rK;
        H.curveFitMetadata.push_back(alphaK_rK);
    }
    // change this
    populateShuffler(shuffler);
    xxILP = std::make_unique<GRBVar[]>(NTuples);
    for (int i = 0; i < NTuples; i++)
    {
        // up to repeat, and fix GRB_BINARY GRB_INTEGER
        xxILP[i] = modelILP.addVar(0.0, 1.0, 0.0, GRB_BINARY, "xx[" + to_string(i) + "]");
    }

    this->epsilon = epsilon;
    cout << "Success Constructor" << endl;
};