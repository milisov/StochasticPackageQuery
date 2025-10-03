#include <iostream>
#include "solver/solversCore.hpp"
#include "CurveFit.hpp"
#include "spq/ssformulator.hpp"
#include "spq/rsformulator.hpp"
#pragma once
using namespace std;

inline void printObjective(GRBModel &model)
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

inline void printExpression(const GRBLinExpr &expr)
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

// If the constraint is of form P(sth)>=p rk = satisfied/numScenarios - p
// else we need to convert the constraint in the right format

template <typename T1, typename T2>
inline void printVectPair(std::vector<pair<T1, T2>> &v)
{
    for (int i = 0; i < v.size(); i++)
    {
        cout << "(" << v[i].first << ", " << v[i].second << ") ";
    }
    cout << endl;
}

inline double getValueQtile(vector<double> &v, int qtileIdx)
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

class History
{
public:
    std::vector<int> xBest;
    double wBest;
    double epsilonBest;
    bool foundFeasible = false;
    // alpha_k, r_k pairs for each probCons K
    std::vector<std::vector<pair<double, double>>> curveFitMetadata;

    History() {};
};

class SummarySearch : public Solver
{
private:
    std::vector<std::vector<std::pair<double, double>>> history;
    std::string fitFunction = "atan";
    double W0;
    double epsilonQ;

public:
    Data &data;
    History H;
    CurveFitter fitter;
    double epsilon;
    std::string SPQPath;

    // Constructor
    SummarySearch(int M = 1e4,
                  std::shared_ptr<StochasticPackageQuery> spq = nullptr,
                  double epsilon = 1e-5) : data(Data::getInstance())
    {
        this->M = M;
        this->spq = spq;
        this->DB_optim = spq->tableName;
        this->DB_valid = fmt::format("{}_{}", DB_optim, "validate");
        this->NTuples = pg.getTableSize(spq->tableName);
        this->cntScenarios = pg.getColumnLength(spq->tableName, "profit");
        this->probConstCnt = countProbConst(spq);
        this->epsilon = epsilon;
    }

    // Extra methods specific to SummarySearch
    void guessOptimalConservativeness(std::vector<std::vector<std::pair<double, double>>> &history,
                                      std::vector<double> &alpha);

    void Best(std::shared_ptr<Objective> obj,
              std::vector<int> &x,
              History &H);

    template <typename T>
    SolutionMetadata<T> CSASolveBinSearch(GRBModel &model,
                                          std::vector<T> &x,
                                          Formulator &formulator,
                                          FormulateOptions &formOptions);

    template <typename T>
    SolutionMetadata<T> CSASolveBinSearchRS(Formulator &formulator,
                                            FormulateOptions &formOptions);

    void solveLP2(GRBModel &model,
                  std::vector<double> &sol,
                  GRBVar *xx,
                  double ub,
                  std::vector<int> dummyVect);

    template <typename T>
    SolutionMetadata<T> summarySearch(shared_ptr<StochasticPackageQuery> spq, Formulator &formulator, FormulateOptions &formOptions, map<string, Option> &curveFitOptions, int z);

    template <typename T>
    SolutionMetadata<T> summarySearchRS(shared_ptr<StochasticPackageQuery> spq, Formulator &formulator,
                                        FormulateOptions &formOptions, map<string, Option> &curveFitOptions, int z);

    SolutionMetadata<int> stochDualReducer(std::shared_ptr<StochasticPackageQuery> spq,
                                           FormulateOptions &formOptions,
                                           std::map<std::string, Option> &curveFitOptions);
};

// template <typename T>
// SolutionMetadata<T> SummarySearch::CSASolveBinSearch(GRBModel &model, std::vector<T> &x, Formulator &formulator, FormulateOptions &formOptions)
// {
//     SolveOptions options;
//     BinarySearchMetadata alpha_KMetadata(1e-5, 1.0, 0.0);
//     std::vector<BinarySearchMetadata> history;
//     initializeVector(history, probConstCnt, alpha_KMetadata);
//     std::vector<double> alpha;
//     initializeVector(alpha, probConstCnt, -1.0);
//     int q = 0;
//     int qAfterZequalsM = 0;

//     SolutionMetadata bestSol(x, 0, 0, false);

//     formulator.partitions.clear();
//     formulator.reshuffleShuffler(formulator.shuffler);
//     for(int i = 0; i < formulator.spq->cons.size(); i++)
//     {
//         int conOrder = 0;
//         shared_ptr<ProbConstraint> probCon;
//         shared_ptr<AttrConstraint> attrCon;
    
//         bool isstoch = isStochastic(formulator.spq->cons[i], probCon, attrCon);
//         if (isstoch)
//         {
//             formulator.partition(formOptions.Z, formOptions.innerConstraints[conOrder], formulator.shuffler, formulator.partitions);
//             conOrder++;
//         }
//     }

//     while (true)
//     {
//         int Z = formOptions.Z;
//         // cout << "Iteration: " << q << " Z = " << Z << endl;
//         //  validate() -> this will set the rk values, and calculate the Wq
//         if (x.size() > 0)
//         {
//             validate(model, x, this->spq, options);
//             bool isFeas = isFeasible(r);
//             if (q == 0)
//             {
//                 W0 = W_q;
//             }
//             // loop through rk values and check if they are all >= 0 -> if true
//             // calculate the e^Q -> if <= epsilon
//             double epsilonQ = calculateEpsilonQ(this->spq, this->W_q, this->W0);
//             // cout << "EPSILON = " << epsilonQ << endl;
//             // cout << (isFeas ? "Feasible" : "Infeasible") << endl;
//             // cout << epsilon << " " << (epsilonQ <= this->epsilon ? "Good Bound" : "Bad Bound") << endl;

//             //-> return this solution -> this means update xBest to be equal to X, and set all metadata
//             if (isFeas && epsilonQ <= this->epsilon)
//             {
//                 bestSol.setSolution(x, W_q, epsilonQ, true, true, Z);
//                 return bestSol;
//             }
//             else if (isFeas)
//             {
//                 bestSol.setSolution(x, W_q, epsilonQ, false, false, Z);
//             }
//         }

//         q = q + 1;
//         if (Z == cntScenarios)
//         {
//             qAfterZequalsM += 1;
//         }

//         // call the binary search function with for each constraint
//         for (int i = 0; i < probConstCnt; i++)
//         {
//             int ScenariosLeft = (int)ceil(history[i].low * (this->cntScenarios / Z));
//             int ScenariosRight = (int)ceil(history[i].high * (this->cntScenarios / Z));
//             double eps = 1e-3;
//             if (history[i].high - history[i].low > eps && (Z <= cntScenarios && qAfterZequalsM <= 1))
//             {
//                 // cout << "Before call of Binary Search: " << history[i].low << " " << history[i].high << " " << alpha[i] << endl;
//                 if (alpha[i] != -1.0)
//                 {
//                     if (x.size() == 0)
//                     {
//                         // the solution is feasible but suboptimal or the system is infeasible -> use less conservative summary
//                         // cout << "USE LESS CONSERVATIVE" << endl;
//                         history[i].high = alpha[i];
//                     }
//                     else
//                     {
//                         if (r[i] < 0)
//                         {
//                             // cout << "USE MORE CONSERVATIVE" << endl;
//                             //  the solution is infeasible -> use more conservative summary
//                             history[i].low = alpha[i];
//                         }
//                         else
//                         {
//                             // the solution is feasible but suboptimal or the system is infeasible -> use less conservative summary
//                             // cout << "USE LESS CONSERVATIVE" << endl;
//                             history[i].high = alpha[i];
//                         }
//                     }
//                 }
//                 alpha[i] = history[i].low + (history[i].high - history[i].low) / 2;
//                 // cout << "After call of Binary Search: " << history[i].low << " " << history[i].high << " " << alpha[i] << endl;
//             }
//             else
//             {
//                 // cout << "BINARY SEARCH END CONDITION MET" << endl;
//                 return bestSol;
//             }
//         }
//         formOptions.alpha = alpha;
//         formOptions.iteration = q;
//         formOptions.innerConstraints = this->innerConstraints;
//         GRBModel model = formulator.formulate(spq, formOptions);
//         model.update();
//         initializeVector(x, NTuples, T(0));
//         options.reduced = formOptions.reduced;
//         options.reducedIds = formOptions.reducedIds;
//         options.computeActiveness = false;
//         solve(model, x, options);
//     }
// }

template <typename T>
SolutionMetadata<T> SummarySearch::CSASolveBinSearchRS(Formulator &formulator, FormulateOptions &formOptions)
{
    std::vector<T> x;
    BinarySearchMetadata alpha_KMetadata(1e-5, 1.0, 0.0);
    std::vector<BinarySearchMetadata> history;
    initializeVector(history, probConstCnt, alpha_KMetadata);
    std::vector<double> alpha;
    initializeVector(alpha, probConstCnt, -1.0);

    int q = 0;
    int qAfterZequalsM = 0;
    std::vector<int> vbasis = formOptions.vbasis;
    std::vector<int> cbasis = formOptions.cbasis;
    SolutionMetadata bestSol(x, 0, 0);
    int Z = formOptions.Z;


    formulator.partitions.clear();
    formulator.reshuffleShuffler(formulator.shuffler);
    for(int i = 0; i < formulator.spq->cons.size(); i++)
    {
        int conOrder = 0;
        shared_ptr<ProbConstraint> probCon;
        shared_ptr<AttrConstraint> attrCon;
    
        bool isstoch = isStochastic(formulator.spq->cons[i], probCon, attrCon);
        if (isstoch)
        {
            formulator.partition(Z, formOptions.innerConstraints[conOrder], formulator.shuffler, formulator.partitions);
            conOrder++;
        }
    }
    while (true)
    {
        cout << "Iteration: " << q << " Z = " << Z << " "<<qAfterZequalsM << " "<<cntScenarios<<endl;

        if (Z == cntScenarios)
        {
            qAfterZequalsM += 1;
        }

        for (int i = 0; i < probConstCnt; i++)
        {
            alpha[i] = history[i].low + (history[i].high - history[i].low) / 2;
        }
        formOptions.alpha = alpha;
        deb(alpha);
        // deb(formOptions.alpha);
        formOptions.iteration = q;
        // formulate the I/LP
        GRBModel model = formulator.formulate(spq, formOptions);
        // optimize and store solution in x
        model.update();
        //model.write("/home/fm2288/StochasticPackageQuery/src/solver/try+to_str(q).lp");
        initializeVector(x, NTuples, T(0));
        SolveOptions options;
        options.reduced = formOptions.reduced;
        options.reducedIds = formOptions.reducedIds;
        options.computeActiveness = formOptions.computeActiveness;
        if (!vbasis.empty() && !cbasis.empty())
        {
            GRBVar *xx = model.getVars();
            for (size_t i = 0; i < vbasis.size(); ++i)
            {
                xx[i].set(GRB_IntAttr_VBasis, vbasis[i]);
            }
            
            GRBConstr *constrs = model.getConstrs();
            for (size_t i = 0; i < cbasis.size(); ++i)
            {
                constrs[i].set(GRB_IntAttr_CBasis, cbasis[i]);
            }
        }
        solve(model, x, options);
        if (x.size() > 0)
        {
            validate(model, x, this->spq, options);
            formOptions.innerConstraints = this->innerConstraints;
            bool isFeas = isFeasible(r);
            if (q == 0)
            {
                W0 = W_q;
            }
            // loop through rk values and check if they are all >= 0 -> if true
            cout << (isFeas ? "Feasible" : "Infeasible") << endl;

            //-> return this solution -> this means update xBest to be equal to X, and set all metadata
            if (isFeas)
            {
                bestSol.setSolution(x, W_q, true, true, Z);
                bestSol.bestRk = r[0];
                bestSol.bestPosActivenessRS = posActivenessRS;
                bestSol.bestNegActivenessRS = negActivenessRS;
                return bestSol;
            }
            else if(r[0] > bestSol.bestRk)
            {
                //We abuse the fact that we know it's only 1 VaR constr;
                bestSol.setSolution(x,W_q, false, true, Z);
                bestSol.bestRk = r[0];
                bestSol.bestPosActivenessRS = posActivenessRS;
                bestSol.bestNegActivenessRS = negActivenessRS;
            }

            vbasis.resize(model.get(GRB_IntAttr_NumVars));
            cbasis.resize(model.get(GRB_IntAttr_NumConstrs));

            GRBVar *vars = model.getVars();
            int numVars = model.get(GRB_IntAttr_NumVars);
            vbasis.resize(numVars);
            for (int i = 0; i < numVars; ++i)
            {
                vbasis[i] = vars[i].get(GRB_IntAttr_VBasis);
            }

            GRBConstr *constrs = model.getConstrs();
            int numConstrs = model.get(GRB_IntAttr_NumConstrs);
            cbasis.resize(numConstrs);
            for (int i = 0; i < numConstrs; ++i)
            {
                cbasis[i] = constrs[i].get(GRB_IntAttr_CBasis);
            }
            formOptions.vbasis = vbasis;
            formOptions.cbasis = cbasis;
        }

        for (int i = 0; i < probConstCnt; i++)
        {
            int ScenariosLeft = (int)ceil(history[i].low * (this->cntScenarios / Z));
            int ScenariosRight = (int)ceil(history[i].high * (this->cntScenarios / Z));
            double eps = 1e-3;
            if (history[i].high - history[i].low > eps && (Z <= cntScenarios && qAfterZequalsM < 1))
            {
                if (alpha[i] != -1.0)
                {
                    if (x.size() == 0)
                    {
                        // if No Solution try smaller alpha
                        history[i].high = alpha[i];
                    }
                    else
                    {
                        if (r[i] < 0)
                        {
                            // if infeasible, go for more conservative alpha
                            history[i].low = alpha[i];
                        }
                        else
                        {
                            // if feasible but suboptimal try less conservative alpha
                            history[i].high = alpha[i];
                        }
                    }
                }
            }
            else
            {
                return bestSol;
            }
        }
        q = q + 1;
    }
}

template <typename T>
SolutionMetadata<T> SummarySearch::summarySearchRS(shared_ptr<StochasticPackageQuery> spq, Formulator &formulator,
                                                   FormulateOptions &formOptions, map<string, Option> &curveFitOptions, int z)
{
    std::vector<std::vector<std::vector<double>>> summaries;
    formOptions.M = this->M;
    bool binSearch = boost::get<bool>(curveFitOptions.at("binarySearch"));
    bool curveFit = boost::get<bool>(curveFitOptions.at("arctan"));
    std::vector<std::vector<std::pair<int, double>>> innerConstraintsDet = formOptions.innerConstraints;
    SolutionMetadata<T> sol;
    formOptions.innerConstraints = innerConstraintsDet;
    sol = CSASolveBinSearchRS<T>(formulator, formOptions);
    deb(sol.x);
    return sol;
}

// template <typename T>
// SolutionMetadata<T> SummarySearch::summarySearch(shared_ptr<StochasticPackageQuery> spq, Formulator &formulator, FormulateOptions &formOptions, map<string, Option> &curveFitOptions, int z)
// {
//     std::vector<std::vector<std::vector<double>>> summaries;
//     formOptions.M = this->M;
//     // formulate Deterministic ILP
//     GRBModel model = formulator.formulate(spq, formOptions); // need to add the right values of formOptions from wherever you are calling summary search
//     int numConstraints = model.get(GRB_IntAttr_NumConstrs);
//     // printVariableNames(model);
//     std::vector<T> x0;
//     initializeVector(x0, NTuples, T(0));
//     SolveOptions options;
//     options.reduced = formOptions.reduced;
//     options.reducedIds = formOptions.reducedIds;
//     solve(model, x0, options);

//     bool binSearch = boost::get<bool>(curveFitOptions.at("binarySearch"));
//     bool curveFit = boost::get<bool>(curveFitOptions.at("arctan"));

//     formOptions.Z = 1;
//     while (true)
//     {
//         std::vector<T> x = x0; // copy solution to deterministic
//         SolutionMetadata<T> sol;
//         if (binSearch)
//         {
//             sol = CSASolveBinSearch(model, x, formulator, formOptions);
//         }
//         else
//         {
//             // sol = CSASolve(x, M, Z, reducedIds, reduced, cntoptions);
//         }
//         if (sol.isFeasible && sol.epsilon <= this->epsilon)
//         {
//             return sol;
//         }
//         else
//         {
//             if (cntScenarios == formOptions.Z)
//             {
//                 std::vector<T> xxx;
//                 SolutionMetadata<T> sol(xxx, 0, 0, false);
//                 sol.Z = cntScenarios;
//                 return sol;
//             }
//             formOptions.Z = formOptions.Z + min(z, cntScenarios - formOptions.Z);
//             z = z * 2;
//         }
//     }
// }