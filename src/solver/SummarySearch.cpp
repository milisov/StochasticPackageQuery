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

    int numGenConstrs = model.get(GRB_IntAttr_NumGenConstrs);
    cout << "NUMBER OF GENERAL CONSTRAINTS ADDED: " << numGenConstrs << endl;

    GRBGenConstr *genConstrs = model.getGenConstrs();
    for (int i = 0; i < numGenConstrs; i++)
    {
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
        }
    }
}

void printObjective(GRBModel &model)
{
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
T1 calculateSummary(vector<pair<T1, T2>> &vect, bool maxS, double alpha)
{
    int G_alpha = (int)ceil(alpha * (double)vect.size());
    T1 summary = vect[0].first;
    // cout << "G(alpha)=" << G_alpha << endl;
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
void initializeVector(vector<T> &v, int sz, T init)
{
    for (int i = 0; i < sz; i++)
    {
        v.push_back(init);
    }
}

// function that finds the indices of solution x that are non-zero
void findNonzero(vector<int> &selectIds, vector<int> &x)
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
int countSatisfied(int numScenarios, vector<double> &innerConst, shared_ptr<ProbConstraint> probCon, shared_ptr<StochasticPackageQuery> spq)
{
    int satisfied = 0;
    for (int j = 0; j < numScenarios; j++)
    {
        if (probCon->vsign == Inequality::gteq) // >=
        {
            if (innerConst[j] >= spq->getValue(probCon->v))
            {
                satisfied++;
            }
        }
        else
        {
            if (innerConst[j] <= spq->getValue(probCon->v)) // <=
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

void SummarySearch::validate(vector<int> &x, shared_ptr<StochasticPackageQuery> spq, int M_hat, string DB_optim)
{
    // empty previous surplus vector
    cout << "Validating!" << endl;
    // clear r from previous iteration -> holds all of the rk values for each iteration
    this->r.clear();
    // we need to clear the innerConstraints in order to update them with new values
    this->innerConstraints.clear();
    // because we get row by row we need a way to compute the innerConst for all scenarios in the same time
    int cons_num = spq->cons.size();

    vector<int> selectIds;
    findNonzero(selectIds, x);
    // WARNING: NEED TO DELETE THE NEXT 2 LINES: THIS IS JUST FOR TESTING THE SELECTION QUERY
    selectIds.push_back(1);
    selectIds.push_back(2);
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
            // initialize innerConst with 0 in order to computer the innerConst for that particular probCons
            vector<double> innerConst;
            initializeVector(innerConst, M, (double)0);
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
            int row;
            int numScenarios = 0;
            while (sr.fetchRow())
            {
                vector<double> realization;
                sr.getArray(1, realization);
                int row = sr.getBigInt(0) - 1;
                numScenarios = realization.size();
                for (int j = 0; j < numScenarios; j++)
                {
                    innerConst[j] += x[row] * realization[j];
                }
            }
            // push innerConst to innerConstraints
            this->innerConstraints.push_back(innerConst);
            int satisfied = countSatisfied(numScenarios, innerConst, probCon, spq);
            double rk = calculateRk(probCon, satisfied, numScenarios, spq);
            r.push_back(rk);
            // compute v
        }
    }

    cout << "Validation Completed" << endl;
};

// get the vector pairs with realizations and scores
// shuffle them and then find what is partitionSize
// start from beginning, count partitionSize numbers and add them into a partition
// once >= partitionSize reset cnt and start new partition add the prev into partitions
void partition(int Z, vector<pair<double, double>> &realizationScore, vector<int> shuffler, vector<vector<pair<double, double>>> &partitions)
{
    int totalLength = realizationScore.size();
    double partLength = ceil(totalLength / (double)Z);

    int partitionSize = static_cast<int>(partLength);
    // create a partition vector of pairs that contains the realization and scores
    vector<pair<double, double>> partition;
    int cnt = 0;
    for (int i = 0; i < totalLength; i++)
    {
        // use the shuffler array which is already shuffled in order to retrieve partitionSize number of RANDOM elements into the partition
        partition.push_back(realizationScore[shuffler[i]]);
        cnt++;
        if (cnt >= partitionSize)
        {
            // push partition to partition
            // clear the current partition to calculate the next one
            // reset counter to 0 to count for the next partition
            cnt = 0;
            partitions.push_back(partition);
            partition.clear();
        }
    }
}

template <typename T1, typename T2>
void printVectPair(vector<pair<T1, T2>> &v)
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
vector<vector<double>> SummarySearch::summarize(vector<int> &x, int Z, double alpha, string DB_optim, shared_ptr<ProbConstraint> probCon, shared_ptr<AttrConstraint> attrCon, int conOrder)
{
    string selectcols = attrCon->attr;
    // I have put limit 1 to test easier REMOVE LIMIT 1 for future use
    string sql = fmt::format("SELECT {} FROM \"{}\"", selectcols, DB_optim);
    SingleRow sr = SingleRow(sql);
    // we have a vector of doubles, each representing the summary for each partition
    vector<vector<double>> summaries;
    vector<double> summary;
    bool init;

    while (sr.fetchRow())
    {
        vector<double> realization;
        vector<pair<double, double>> realizationScore;
        sr.getArray(0, realization);

        int numScenarios = realization.size();
        // cout << "NUMBER OF SCENARIOS:" << numScenarios << endl;
        //  make pairs of realization of tuple read by fetchRow(), and the corresponding score of x for scenario i
        for (int i = 0; i < numScenarios; i++)
        {
            realizationScore.push_back(make_pair(realization[i], this->innerConstraints[conOrder][i]));
            // cout << "(" << realization[i] << ", " << innerConstraints[conOrder][i] << ") ";
        }
        // cout << endl;

        // partition
        vector<vector<pair<double, double>>> partitions;
        partition(Z, realizationScore, this->shuffler, partitions);
        // summaries is initialized with a number of partitions, summary vectors
        if (!init)
        {
            initializeVector(summaries, partitions.size(), summary);
            init = true;
        }

        // partition by partition, sort
        // summarize
        for (int i = 0; i < partitions.size(); i++)
        {
            // cout << "Before Sort - (Realization, Score)" << endl;
            // printVectPair(partitions[i]);
            double summary;
            // decide when ASC when DESC
            if (probCon->vsign == Inequality::gteq) //>= min summary
            {
                sort(partitions[i].begin(), partitions[i].end(), sortbysecDESC);
                summary = calculateSummary(partitions[i], false, alpha);
            }
            else
            {
                sort(partitions[i].begin(), partitions[i].end(), sortbysecASC);
                summary = calculateSummary(partitions[i], true, alpha);
            }

            // cout << "After Sort - (Realization, Score)" << endl;
            // printVectPair(partitions[i]);
            // cout << "The summary is: " << summary << endl;
            summaries[i].push_back(summary);
            // I need to define a vector of vector for summaries and push the summary value into the appropriate place
            // return the vector of vectors
        }
    }
    return summaries;
}

void SummarySearch::reshuffleShuffler(vector<int> &shuffle)
{
    srand(unsigned(time(0)));
    random_shuffle(shuffle.begin(), shuffle.end());
}
/*fetch the attr from the DB and add the corresponding value to the coeff[i]
then formulate the linearExpression and based on the objSense add Maximize or Minimize
*/
void SummarySearch::formSumObj(shared_ptr<Objective> obj, GRBVar *xx)
{
    shared_ptr<AttrObjective> attrObj;
    bool isDet = isDeterministic(obj, attrObj);
    if (!isDet || attrObj->objType != numeric_type)
    {
        return;
    }

    cout << "Sum Objective" << endl;

    string selectcols = attrObj->obj;
    string sql = fmt::format("SELECT {} FROM {}", selectcols, DB_optim);

    SingleRow sr = SingleRow(sql);
    double coeff[NTuples];
    int i = 0;

    while (sr.fetchRow())
    {
        coeff[i] = sr.getNumeric(0);
        i++;
    }

    GRBLinExpr sumObjExpr;
    sumObjExpr.addTerms(coeff, xx, NTuples);

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
void SummarySearch::formExpSumObj(shared_ptr<Objective> obj, GRBVar *xx)
{
    shared_ptr<AttrObjective> attrObj;
    bool isDet = isDeterministic(obj, attrObj);
    if (!isDet || attrObj->objType != array_type)
    {
        return;
    }
    cout << "Expected Sum Objective" << endl;

    string table = fmt::format("{}_{}", DB_optim, "summary");
    string selectcols = fmt::format("{}_{}", attrObj->obj, "mean");
    string sql = fmt::format("SELECT {} FROM {}", selectcols, table);
    SingleRow sr = SingleRow(sql);
    double coeff[NTuples];
    int i = 0;

    while (sr.fetchRow())
    {
        coeff[i] = sr.getNumeric(0);
        i++;
    }

    GRBLinExpr expSumObjExpr;
    expSumObjExpr.addTerms(coeff, xx, NTuples);

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
void SummarySearch::formCntObj(shared_ptr<Objective> obj, GRBVar *xx)
{
    shared_ptr<CountObjective> cntObj = getCount(obj);
    if (!cntObj)
    {
        return;
    }

    cout << "Count Objective" << endl;

    double coeff[NTuples];

    for (int i = 0; i < NTuples; i++)
    {
        coeff[i] = 1;
    }

    GRBLinExpr cntObjExpr;
    cntObjExpr.addTerms(coeff, xx, NTuples);

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
void SummarySearch::formCountCons(shared_ptr<Constraint> cons, GRBVar *xx)
{
    shared_ptr<CountConstraint> cntCons = getCount(cons);
    double ub_eps = 1e30;
    double lb_eps = -1e30;

    if (!cntCons)
    {
        return;
    }

    GRBLinExpr cntConsExpr;
    double coeff[NTuples];
    for (int i = 0; i < NTuples; i++)
    {
        coeff[i] = 1;
    }
    cntConsExpr.addTerms(coeff, xx, NTuples);

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
}
//select the attr from the DB, and then add the values of each tuple to the corresponding coeff[i]
//formulate GRBLinExpr with addTerms.(coeff,xx,NTuples)
//based on the bounds formulate the constraints ignore the bounds that are not set
void SummarySearch::formSumCons(shared_ptr<Constraint> cons, GRBVar *xx)
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
    GRBLinExpr sumConsExpr;

    string selectcols = attrCon->attr;
    string sql = fmt::format("SELECT {} FROM \"{}\"", selectcols, DB_optim);
    SingleRow sr = SingleRow(sql);
    double coeff[NTuples];
    int i = 0;
    while (sr.fetchRow())
    {
        coeff[i] = sr.getNumeric(0);
        i++;
    }

    sumConsExpr.addTerms(coeff, xx, NTuples);

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
}

//go through the constraints and delete them
void removeProbConstr(GRBModel &model, vector<GRBVar> &yy, vector<GRBGenConstr> &genCon, vector<GRBConstr> &sumyCon)
{
    for (int i = 0; i < genCon.size(); i++)
    {
        cout << "deleting Indicator" << endl;
        model.remove(genCon[i]);
    }
    model.update();
    for (int i = 0; i < sumyCon.size(); i++)
    {
        cout << "deleting SUM yy[i]" << endl;
        model.remove(sumyCon[i]);
    }
    model.update();
    for (int i = 0; i < yy.size(); i++)
    {
        cout << "deleting a yy[i]" << endl;
        model.remove(yy[i]);
    }
    model.update();
    yy.clear();
    genCon.clear();
    sumyCon.clear();
}
//check if constraint is probabilistic
//find the number of summaries for that constraint
//remove the variables and constraints related to this prob constraint if any
//find pZ for the SumY
//add indicator constraints for the variables Y[z]
void SummarySearch::formProbCons(shared_ptr<Constraint> cons, GRBVar *xx, vector<vector<vector<double>>> &summaries, int &probConOrder, vector<GRBVar> &yy, vector<GRBGenConstr> &genCon, vector<GRBConstr> &sumyCon)
{
    shared_ptr<ProbConstraint> probCon;
    shared_ptr<AttrConstraint> attrCon;

    bool isstoch = isStochastic(cons, probCon, attrCon);
    if (!isstoch)
    {
        return;
    }
    removeProbConstr(this->model, yy, genCon, sumyCon);
    model.update();
    // calculate the number of yk indicator variables
    int Z = summaries[probConOrder].size();
    double v = spq->getValue(probCon->v);
    double p = spq->getValue(probCon->p);
    double pZ = ceil(p * Z);
    GRBVar y[Z];
    double coeff_Y[Z];
    for (int z = 0; z < Z; z++)
    {
        y[z] = model.addVar(0.0, 0.1, 0.0, GRB_BINARY, "y[" + to_string(z) + "]");
        yy.push_back(y[z]);
        coeff_Y[z] = 1;
    }
    model.update();
    GRBLinExpr sumYz;
    vector<vector<double>> S = summaries[probConOrder];
    cout << probConOrder << " " << summaries.size() << endl;
    for (int z = 0; z < S.size(); z++)
    {
        GRBLinExpr innerCons;
        int sz = S[z].size();
        double coeff[sz];
        cout << sz << endl;
        for (int i = 0; i < sz; i++)
        {
            coeff[i] = S[z][i];
        }
        innerCons.addTerms(coeff, xx, NTuples);
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
/*
check if constraint has bounds, an attribute and is of array_type -> expCons
we can make use of the percomputed table for tuple-wise mean -> determine the table by concatenating "summary"
select cols = concatenate "mean"
assign coeff[NTuples] the tuple-wise means, and then formulate the linearExpr with the coeff and the xx
add constraint based on the bounds
*/
void SummarySearch::formExpCons(shared_ptr<Constraint> cons, GRBVar *xx)
{
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
    string selectcols = fmt::format("{}_{}", attrCon->attr, "mean");
    string sql = fmt::format("SELECT {} FROM  \"{}\"", selectcols, table);

    SingleRow sr = SingleRow(sql);
    double coeff[NTuples];
    int i = 0;
    while (sr.fetchRow())
    {
        coeff[i] = sr.getNumeric(0);
        i++;
    }

    GRBLinExpr expSumConsExpr;
    expSumConsExpr.addTerms(coeff, xx, NTuples);

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
}

//naive
//formulateSAA(summary[M],q=0);
void SummarySearch::formulateSAA(vector<vector<vector<double>>> &summaries, int q)
{
    cout << "formulating" << endl;
    int numCons = spq->cons.size();
    int probConOrder = 0;
    for (int i = 0; i < numCons; i++)
    {
        if (q == 0)
        {
            formCountCons(spq->cons[i], xx.get());
            formSumCons(spq->cons[i], xx.get());
            formProbCons(spq->cons[i], xx.get(), summaries, probConOrder, yy, genCon, sumyCon);
            formExpCons(spq->cons[i], xx.get());
        }
        else
        {
            formProbCons(spq->cons[i], xx.get(), summaries, probConOrder, yy, genCon, sumyCon);
        }
    }
    if (q == 0)
    {
        formSumObj(spq->obj, xx.get());
        formExpSumObj(spq->obj, xx.get());
        formCntObj(spq->obj, xx.get());
    }
    model.update();
    // model.write("debug.lp");
    printConstraints(this->model);
    // printObjective(this->model);
}

SummarySearch::SummarySearch(int M, int M_hat, shared_ptr<StochasticPackageQuery> spq)
{
    this->M = M;
    this->M_hat = M_hat;
    this->spq = spq;
    // get NTuples from spq->tableName and populate X with NTuples zeroes
    this->DB_optim = spq->tableName;
    // test
    this->DB_valid = fmt::format("{}_{}", DB_optim, "validate");
    this->NTuples = pg.getTableSize(spq->tableName);
    this->probConstCnt = countProbConst(spq);
    for (int i = 0; i < M; i++)
    {
        shuffler.push_back(i);
    }
    // populate solution vector
    initializeVector(x, NTuples, 0);
    xx = std::make_unique<GRBVar[]>(NTuples);
    for (int i = 0; i < NTuples; i++)
    {
        //up to repeat, and fix GRB_BINARY GRB_INTEGER
        xx[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "xx[" + to_string(i) + "]");
    }
    cout << "Success Constructor" << endl;
};