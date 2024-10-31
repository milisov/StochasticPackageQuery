#include "Naive.hpp"
#include <boost/algorithm/string/join.hpp>
#include <fmt/ranges.h>

void printNaiveConstraints(GRBModel &model)
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
    for (int i = 0; i < 3; i++)
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

void printNaiveObjective(GRBModel &model)
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




// function that initializes given vector
template <typename T>
void initializeVector(std::vector<T> &v, int sz, T init)
{
    for (int i = 0; i < sz; i++)
    {
        v.push_back(init);
    }
}


/*fetch the attr from the DB and add the corresponding value to the coeff[i]
then formulate the linearExpression and based on the objSense add Maximize or Minimize
*/
void Naive::formSumObj(shared_ptr<Objective> obj, GRBVar *xx)
{
    shared_ptr<AttrObjective> attrObj;
    bool isDet = isDeterministic(obj, attrObj);
    if (!isDet || attrObj->objType != numeric_type)
    {
        return;
    }

    cout << "Sum Objective" << endl;

    string selectcols = fmt::format("{},{}", "id", attrObj->obj);
    string sql = fmt::format("SELECT {} FROM {}", selectcols, DB_optim);

    SingleRow sr = SingleRow(sql);
    double coeff[NTuples];

    while (sr.fetchRow())
    {
        int id = sr.getBigInt(0) - 1; 
        double coeffVal = sr.getNumeric(1);
        coeff[id] = coeffVal;
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
void Naive::formExpSumObj(shared_ptr<Objective> obj, GRBVar *xx)
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
    string sql = fmt::format("SELECT {} FROM {}", selectcols, table);
    SingleRow sr = SingleRow(sql);
    double coeff[NTuples];

    while (sr.fetchRow())
    {
        int id = sr.getBigInt(0) - 1;
        double coeffVal = sr.getNumeric(1);
        coeff[id] = coeffVal;
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
void Naive::formCntObj(shared_ptr<Objective> obj, GRBVar *xx)
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
void Naive::formCountCons(shared_ptr<Constraint> cons, GRBVar *xx)
{
    shared_ptr<CountConstraint> cntCons = getCount(cons);
    double ub_eps = 1e30;
    double lb_eps = -1e30;

    if (!cntCons)
    {
        return;
    }
    cout<<"Count Constraint"<<endl;
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
// select the attr from the DB, and then add the values of each tuple to the corresponding coeff[i]
// formulate GRBLinExpr with addTerms.(coeff,xx,NTuples)
// based on the bounds formulate the constraints ignore the bounds that are not set
void Naive::formSumCons(shared_ptr<Constraint> cons, GRBVar *xx)
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

    string selectcols = fmt::format("{},{}", "id", attrCon->attr);
    // string selectcols = attrCon->attr;
    string sql = fmt::format("SELECT {} FROM \"{}\"", selectcols, DB_optim);
    SingleRow sr = SingleRow(sql);
    double coeff[NTuples];
    while (sr.fetchRow())
    {
        int id = sr.getBigInt(0) - 1;
        double coeffVal = sr.getNumeric(1);
        coeff[id] = coeffVal;
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


// check if constraint is probabilistic
// find the number of summaries for that constraint
// remove the variables and constraints related to this prob constraint if any
// find pZ for the SumY
// add indicator constraints for the variables Y[z]
void Naive::formProbCons(shared_ptr<Constraint> cons, GRBVar *xx)
{
    shared_ptr<ProbConstraint> probCon;
    shared_ptr<AttrConstraint> attrCon;

    bool isstoch = isStochastic(cons, probCon, attrCon);
    if (!isstoch)
    {
        return;
    }
    cout<<"FORMULATING PROB CONSTRAINT"<<endl;
    // calculate the number of yk indicator variables
    double v = spq->getValue(probCon->v);
    double p = spq->getValue(probCon->p);

    string selectcols = fmt::format("{},{}", "id", attrCon->attr);
    // string selectcols = attrCon->attr;
    string sql = fmt::format("SELECT {} FROM \"{}\" LIMIT 1", selectcols, DB_optim);
    SingleRow sr2 = SingleRow(sql);


    while(sr2.fetchRow())
    {
        vector<double>realization;
        sr2.getArray(1,realization);
        cntScenarios = realization.size();
    }


    double pM = ceil(p * cntScenarios);
    GRBVar y[cntScenarios];
    double coeff_Y[cntScenarios];
    for (int j = 0; j < cntScenarios; j++)
    {
        y[j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y[" + to_string(j) + "]");
        //yy.push_back(y[j]);
        coeff_Y[j] = 1;
    }
    model.update();

    GRBLinExpr sumYz;
    vector<double>vect;
    initializeVector(vect, NTuples, (double) 0.0);
    vector<vector<double>>coeffs;
    initializeVector(coeffs, cntScenarios, vect);
    sql = fmt::format("SELECT {} FROM \"{}\"", selectcols, DB_optim);
    SingleRow sr = SingleRow(sql);
    while(sr.fetchRow())
    {
        int id = sr.getBigInt(0) - 1;
        vector<double>realization;
        sr.getArray(1,realization);
        for(int i = 0; i < realization.size(); i++)
        {
            coeffs[i][id] = realization[i];
        }
    }


    for(int i = 0; i < coeffs.size(); i++)
    {
        GRBLinExpr innerCons;
        double coeff[NTuples];
        for(int j = 0; j < coeffs[i].size(); j++)
        {
            coeff[j] = coeffs[i][j];
        }
        innerCons.addTerms(coeff, xx, NTuples);
        try
        {
            if (probCon->vsign == Inequality::gteq)
            {
                GRBGenConstr indicator = model.addGenConstrIndicator(y[i], 1, innerCons, GRB_GREATER_EQUAL, v);
            }
            else
            {
                GRBGenConstr indicator = model.addGenConstrIndicator(y[i], 1, innerCons, GRB_LESS_EQUAL, v);

            }
        }
        catch (GRBException &e)
        {
            cout << "Error code 8 = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        }
    }

    sumYz.addTerms(coeff_Y, y, cntScenarios);
    try
    {
        if (probCon->psign == Inequality::gteq)
        {
            GRBConstr constr = model.addConstr(sumYz, GRB_GREATER_EQUAL, pM);
        }
        else
        {
            GRBConstr constr = model.addConstr(sumYz, GRB_LESS_EQUAL, pM);
        }
    }
    catch (GRBException &e)
    {
        cout << "Error code 9 = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
}
/*
check if constraint has bounds, an attribute and is of array_type -> expCons
we can make use of the percomputed table for tuple-wise mean -> determine the table by concatenating "summary"
select cols = concatenate "mean"
assign coeff[NTuples] the tuple-wise means, and then formulate the linearExpr with the coeff and the xx
add constraint based on the bounds
*/
void Naive::formExpCons(shared_ptr<Constraint> cons, GRBVar *xx)
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
    string selectAttr = fmt::format("{}_{}", attrCon->attr, "mean");
    string selectcols = fmt::format("{},{}", "id", selectAttr);
    string sql = fmt::format("SELECT {} FROM  \"{}\"", selectcols, table);

    SingleRow sr = SingleRow(sql);
    double coeff[NTuples];
    while (sr.fetchRow())
    {
        int id = sr.getBigInt(0) - 1;
        double coeffVal = sr.getNumeric(1);
        coeff[id] = coeffVal;
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

void Naive::formulateSAA()
{
    cout << "formulating" << endl;
    int numCons = spq->cons.size();
    for (int i = 0; i < numCons; i++)
    {
        formCountCons(spq->cons[i], xx.get());
        formSumCons(spq->cons[i], xx.get());
        formProbCons(spq->cons[i], xx.get());
        formExpCons(spq->cons[i], xx.get());
    }
    formSumObj(spq->obj, xx.get());
    formExpSumObj(spq->obj, xx.get());
    formCntObj(spq->obj, xx.get());
    model.update();
}

inline void findNonzero(std::vector<int> &selectIds, std::vector<int> &x)
{
    for (int i = 0; i < x.size(); i++)
    {
        if (x[i] != 0)
        {
            selectIds.push_back(i + 1);
        }
    }
}

// for each of the Stochastic Constraints we need to check how many are satisfied
// get the string sign from each prob constraint and perform the right operation
// the value of v is retrieved from the spq
inline int countSatisfied(int numScenarios, std::vector<double> &innerConst, shared_ptr<ProbConstraint> probCon, shared_ptr<StochasticPackageQuery> spq)
{
    int satisfied = 0;
    double epsilon = 1e-7;
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

inline double calculateRk(shared_ptr<ProbConstraint> probCon, int satisfied, int numScenarios, shared_ptr<StochasticPackageQuery> spq)
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


inline bool isFeasible(vector<double>&r)
{
    for(int i = 0; i < r.size(); i++)
    {
        if(r[i] < 0)
        {
            return false;
        }
    }
    return true;
}


double Naive::calculateCntObj(std::vector<int> &x, std::vector<int> &selectIds, shared_ptr<CountObjective> cntObj, string DB_optim)
{
    cout << "Calculating Count Objective" << endl;
    double sum = 0;
    for (int i = 0; i < x.size(); i++)
    {
        sum += x[i];
    }
    return sum;
}

double Naive::calculateSumObj(std::vector<int> &x, std::vector<int> &selectIds, shared_ptr<AttrObjective> attrObj, string DB_optim)
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

double Naive::calculateExpSumObj(std::vector<int> &x, std::vector<int> &selectIds, shared_ptr<AttrObjective> attrObj, string DB_optim)
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

double Naive::calculateObj(std::vector<int> &x, std::vector<int> &selectIds, shared_ptr<Objective> obj, string DB_optim)
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


void Naive::validate(std::vector<int> &x, shared_ptr<StochasticPackageQuery> spq, int M_hat, string DB_optim)
{
    cout << "Validating!" << endl;
    // clear r from previous iteration -> holds all of the rk values for each iteration
    this->r.clear();
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
            // because we get row by row we need a way to compute the innerConst for all scenarios in the same time
            std::vector<double> innerConst;
            int row;
            int numScenarios = 0;
            bool initialized;
            // initialize innerConst with 0 in order to computer the innerConst for that particular probCons
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
            int satisfied = countSatisfied(numScenarios, innerConst, probCon, spq);
            double rk = calculateRk(probCon, satisfied, numScenarios, spq);
            this->r.push_back(rk);
        }
    }
    this->W_q = calculateObj(x, selectIds, spq->obj, DB_optim);
    cout << "Validation Completed" << endl;
}


void Naive::solve(std::vector<int> &x, GRBVar *xx)
{
    cout << "Optimizing Model" << endl;
    this->model.optimize();
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
        x[i] = static_cast<int>(xx[i].get(GRB_DoubleAttr_X));
    }

    cout<<"Solution Found"<<endl;
}

Solution Naive::solveNaive(shared_ptr<StochasticPackageQuery> spq, int m, int M, int M_hat)
{
    formulateSAA();
    vector<int>x;
    initializeVector(x,NTuples,0);
    solve(x,this->xx.get());
    validate(x, spq, M_hat, DB_optim);
    deb(r);
    if(isFeasible(r))
    {
        Solution sol(x,this->W_q,true);
        return sol;
    }
    Solution sol(x,-1,false);
    return sol;
}


Naive::Naive(int M, int M_hat, shared_ptr<StochasticPackageQuery> spq)
{
    this->M = M;
    this->M_hat = M_hat;
    this->spq = spq;
    // get NTuples from spq->tableName and populate X with NTuples zeroes
    this->DB_optim = spq->tableName;
    // test
    this->DB_valid = fmt::format("{}_{}", DB_optim, "validate");
    this->NTuples = pg.getTableSize(spq->tableName);

    xx = std::make_unique<GRBVar[]>(NTuples);
    for (int i = 0; i < NTuples; i++)
    {
        // up to repeat, and fix GRB_BINARY GRB_INTEGER
        xx[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "xx[" + to_string(i) + "]");
    }
    cout << "Success Constructor" << endl;
}