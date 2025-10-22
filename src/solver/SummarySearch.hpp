#include <iostream>
#include "solver/solversCore.hpp"
#include "CurveFit.hpp"
#include "spq/ssformulator.hpp"
#include "spq/rsformulator.hpp"
#include <chrono>
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

    // solution vector x, all alpha_k, indicator if it's feasible or not, objective value Wq
    std::vector<std::tuple<std::vector<int>, std::vector<double>, bool, double>> bestSolMetadata;
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
                                          FormulateOptions &formOptions,
                                          const std::chrono::steady_clock::time_point& start_time,
                                          double timeout_seconds);

    template <typename T>
    SolutionMetadata<T> CSASolveBinSearchRS(Formulator &formulator,
                                            FormulateOptions &formOptions);

    void solveLP2(GRBModel &model,
                  std::vector<double> &sol,
                  GRBVar *xx,
                  double ub,
                  std::vector<int> dummyVect);

    template <typename T>
    SolutionMetadata<T> summarySearch(shared_ptr<StochasticPackageQuery> spq, Formulator &formulator, 
        FormulateOptions &formOptions, map<string, Option> &curveFitOptions, int z,
        const std::chrono::steady_clock::time_point& start_time,
        double timeout_seconds);

    template <typename T>
    SolutionMetadata<T> summarySearchRS(shared_ptr<StochasticPackageQuery> spq, Formulator &formulator,
                                        FormulateOptions &formOptions, map<string, Option> &curveFitOptions, int z);

    SolutionMetadata<int> stochDualReducer(std::shared_ptr<StochasticPackageQuery> spq,
                                           FormulateOptions &formOptions,
                                           std::map<std::string, Option> &curveFitOptions);
};

template <typename T>
SolutionMetadata<T> SummarySearch::CSASolveBinSearch(GRBModel &model, std::vector<T> &x, Formulator &formulator, FormulateOptions &formOptions,
                                                     const std::chrono::steady_clock::time_point& start_time,
                                                     double timeout_seconds)
{
    deb("solve csa")
    SolveOptions options;
    BinarySearchMetadata alpha_KMetadata(1e-5, 1.0, 0.0);
    std::vector<BinarySearchMetadata> history;
    initializeVector(history, probConstCnt, alpha_KMetadata);
    std::vector<double> alpha;
    initializeVector(alpha, probConstCnt, -1.0);
    int q = 0;
    int qAfterZequalsM = 0;

    SolutionMetadata<T> bestSol(x, 0, 0, false);
    int Z = formOptions.Z;
    bool partitioned = false;
    while (true)
    {
        cout << "Iteration: " << q << " Z = " << Z << endl;
        //  validate() -> this will set the rk values, and calculate the Wq
        if (x.size() > 0)
        {
            validate(model, x, this->spq, options);
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
                bestSol.setSolution(x, W_q, epsilonQ, true, true, Z);
                return bestSol;
            }else if(isFeas && bestSol.isFeasible)
            {
                cout<<"Best Solve Before:"<<"rk = "<<bestSol.bestRk<<" objective = "<<bestSol.w<<endl;
                if(this->W_q > bestSol.w) //abusing maximization
                {
                    bestSol.x = x;
                    bestSol.isFeasible = true;
                    bestSol.bestRk = r[0];
                    bestSol.w = this->W_q;
                    bestSol.Z = Z;
                    cout<<"Best Solve After:"<<"rk = "<<bestSol.bestRk<<" objective = "<<bestSol.w<<endl;
                }
            }else if(isFeas && !bestSol.isFeasible)
            {
                cout<<"Best Solve Before:"<<"rk = "<<bestSol.bestRk<<" objective = "<<bestSol.w<<endl;
                bestSol.x = x;
                bestSol.isFeasible = true;
                bestSol.bestRk = r[0];
                bestSol.w = this->W_q;
                bestSol.Z = Z;
                cout<<"Best Solve After:"<<"rk = "<<bestSol.bestRk<<" objective = "<<bestSol.w<<endl;
            }else if(!isFeas && !bestSol.isFeasible)
            {
                cout<<"Best Solve Before:"<<"rk = "<<bestSol.bestRk<<" objective = "<<bestSol.w<<endl;
                if(r[0] > bestSol.bestRk)
                {
                    bestSol.x = x;
                    bestSol.isFeasible = false;
                    bestSol.bestRk = r[0];
                    bestSol.w = this->W_q;
                    bestSol.Z = Z;
                    cout<<"Best Solve After:"<<"rk = "<<bestSol.bestRk<<" objective = "<<bestSol.w<<endl;
                }
            }
        }

        auto current_time = std::chrono::steady_clock::now();
        double elapsed_seconds = std::chrono::duration<double>(current_time - start_time).count();
        if (elapsed_seconds > timeout_seconds)
        {
            cout<<"TIME OUT HAPPENED"<<endl;
            return bestSol;
        }

        q = q + 1;
        if (Z == cntScenarios)
        {
            qAfterZequalsM += 1;
        }

        // call the binary search function with for each constraint
        for (int i = 0; i < probConstCnt; i++)
        {
            int ScenariosLeft = (int)ceil(history[i].low * (this->cntScenarios / Z));
            int ScenariosRight = (int)ceil(history[i].high * (this->cntScenarios / Z));
            double eps = 1e-3;
            if (history[i].high - history[i].low > eps && (Z <= cntScenarios && qAfterZequalsM <= 1))
            {
                // cout << "Before call of Binary Search: " << history[i].low << " " << history[i].high << " " << alpha[i] << endl;
                if (alpha[i] != -1.0)
                {
                    if (x.size() == 0)
                    {
                        // the solution is feasible but suboptimal or the system is infeasible -> use less conservative summary
                        // cout << "USE LESS CONSERVATIVE" << endl;
                        history[i].high = alpha[i];
                    }
                    else
                    {
                        if (r[i] < 0)
                        {
                            // cout << "USE MORE CONSERVATIVE" << endl;
                            //  the solution is infeasible -> use more conservative summary
                            history[i].low = alpha[i];
                        }
                        else
                        {
                            // the solution is feasible but suboptimal or the system is infeasible -> use less conservative summary
                            // cout << "USE LESS CONSERVATIVE" << endl;
                            history[i].high = alpha[i];
                        }
                    }
                }
                alpha[i] = history[i].low + (history[i].high - history[i].low) / 2;
                // cout << "After call of Binary Search: " << history[i].low << " " << history[i].high << " " << alpha[i] << endl;
            }
            else
            {
                // cout << "BINARY SEARCH END CONDITION MET" << endl;
                return bestSol;
            }
        }
        formOptions.alpha = alpha;
        formOptions.iteration = q;
        formOptions.innerConstraints = this->innerConstraints;
        if(partitioned == false)
        {
            partitioned = true;
            formulator.createPartitions(formOptions);
        }
        GRBModel model = formulator.formulate(spq, formOptions);
        model.update();
        initializeVector(x, NTuples, T(0));
        options.reduced = formOptions.reduced;
        options.reducedIds = formOptions.reducedIds;
        options.computeActiveness = false;
        solve(model, x, options);
    }
}

template <typename T>
SolutionMetadata<T> SummarySearch::CSASolveBinSearchRS(Formulator &formulator, FormulateOptions &formOptions)
{
    double totalTimeSolving = 0.0;
    std::vector<T> x;
    BinarySearchMetadata alpha_KMetadata(1e-5, 1.0, 0.0);
    std::vector<BinarySearchMetadata> history;
    initializeVector(history, probConstCnt, alpha_KMetadata);
    std::vector<double> alpha;
    initializeVector(alpha, probConstCnt, -1.0);
    int q = 0;
    int qAfterZequalsM = 0;
    std::vector<int> vbasis;
    std::vector<int> cbasis;
    SolutionMetadata bestSol(x, 0, 0, false);
    int Z = formOptions.Z;
    formulator.createPartitions(formOptions);
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
        // deb(formOptions.alpha);
        formOptions.iteration = q;
        // formulate the I/LP
        GRBModel model = formulator.formulate(spq, formOptions);
        // optimize and store solution in x
        model.update();
        //model.write("/home/fm2288/StochasticPackageQuery/src/solver/try" + to_string(q) +".lp");
        initializeVector(x, NTuples, T(0));
        // model.write("mod"+to_string(Z)+"_"+to_string(formOptions.iteration)+".lp");
        SolveOptions options;
        options.reduced = formOptions.reduced;
        options.reducedIds = formOptions.reducedIds;
        options.computeActiveness = formOptions.computeActiveness;

        cout<<vbasis.size()<<" "<<cbasis.size()<<endl;
        if (!vbasis.empty() && !cbasis.empty())
        {
            GRBVar *xx = model.getVars();
            cout<<model.get(GRB_IntAttr_NumVars)<<endl;
            for (size_t i = 0; i < vbasis.size(); ++i)
            {
                xx[i].set(GRB_IntAttr_VBasis, vbasis[i]);
            }
            cout<<model.get(GRB_IntAttr_NumConstrs)<<endl;
            GRBConstr *constrs = model.getConstrs();
            for (size_t i = 0; i < cbasis.size(); ++i)
            {
                constrs[i].set(GRB_IntAttr_CBasis, cbasis[i]);
            }
        }
        cout<<"solving"<<endl;
        solve(model, x, options);

        if (x.size() > 0)
        {
            validate(model, x, this->spq, options);
            formOptions.innerConstraints = this->innerConstraints;
            bool isFeas = isFeasible(r);
            cout << (isFeas ? "Feasible" : "Infeasible") << endl;

            //-> return this solution -> this means update xBest to be equal to X, and set all metadata
            if (isFeas)
            {
                bestSol.x = x;
                bestSol.isFeasible = true;
                bestSol.bestRk = r[0];
                bestSol.bestPosActivenessRS = posActivenessRS;
                bestSol.bestNegActivenessRS = negActivenessRS;
                return bestSol;
            }
            else if(r[0] > bestSol.bestRk)
            {
                bestSol.x = x;
                bestSol.isFeasible = false;
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
            // formOptions.vbasis = vbasis;
            // formOptions.cbasis = cbasis;
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
    std::vector<std::vector<std::pair<int, double>>> innerConstraintsDet = formOptions.innerConstraints;
    SolutionMetadata<T> sol;
    formOptions.innerConstraints = innerConstraintsDet;
    sol = CSASolveBinSearchRS<T>(formulator, formOptions);
    return sol;
}

template <typename T>
SolutionMetadata<T> SummarySearch::summarySearch(shared_ptr<StochasticPackageQuery> spq, Formulator &formulator, 
        FormulateOptions &formOptions, map<string, Option> &curveFitOptions, int z,
        const std::chrono::steady_clock::time_point& start_time,
        double timeout_seconds)
{
    std::vector<std::vector<std::vector<double>>> summaries;
    formOptions.M = this->M;
    // formulate Deterministic ILP
    GRBModel model = formulator.formulate(spq, formOptions); // need to add the right values of formOptions from wherever you are calling summary search
    int numConstraints = model.get(GRB_IntAttr_NumConstrs);
    // printVariableNames(model);
    std::vector<T> x0;
    initializeVector(x0, NTuples, T(0));
    SolveOptions options;
    options.reduced = formOptions.reduced;
    options.reducedIds = formOptions.reducedIds;
    solve(model, x0, options);
    
    SolutionMetadata bestSol(x0, 0, 0, false);
    bool binSearch = boost::get<bool>(curveFitOptions.at("binarySearch"));
    bool curveFit = boost::get<bool>(curveFitOptions.at("arctan"));

    formOptions.Z = 1;
    while (true)
    {
        std::vector<T> x = x0; // copy solution to deterministic
        SolutionMetadata<T> sol;
        if (binSearch)
        {
            sol = CSASolveBinSearch(model, x, formulator, formOptions, start_time, timeout_seconds);
        }
        else
        {
            // sol = CSASolve(x, M, Z, reducedIds, reduced, cntoptions);
        }
        if (sol.isFeasible && sol.epsilon <= this->epsilon)
        {
            return sol;
        }
        else
        {
            if(sol.isFeasible && bestSol.isFeasible)
            {
                cout<<"SUMMARY SEARCH Best Solve Before:"<<"rk = "<<bestSol.bestRk<<" objective = "<<bestSol.w<<endl;
                if(this->W_q > bestSol.w) //abusing maximization
                {
                    bestSol.x = sol.x;
                    bestSol.isFeasible = true;
                    bestSol.bestRk = sol.isFeasible;
                    bestSol.w = sol.w;
                    bestSol.Z = sol.Z;
                    cout<<"SUMMARY SEARCH Best Solve After:"<<"rk = "<<bestSol.bestRk<<" objective = "<<bestSol.w<<endl;
                }
            }else if(sol.isFeasible && !bestSol.isFeasible)
            {
                cout<<"SUMMARY SEARCH Best Solve Before:"<<"rk = "<<bestSol.bestRk<<" objective = "<<bestSol.w<<endl;
                bestSol.x = sol.x;
                bestSol.isFeasible = true;
                bestSol.bestRk = sol.bestRk;
                bestSol.w = sol.w;
                bestSol.Z = sol.Z;
                cout<<"SUMMARY SEEARCH Best Solve After:"<<"rk = "<<bestSol.bestRk<<" objective = "<<bestSol.w<<endl;
            }else if(!sol.isFeasible && !bestSol.isFeasible)
            {
                cout<<"SUMMARY SEARCH Best Solve Before:"<<"rk = "<<bestSol.bestRk<<" objective = "<<bestSol.w<<endl;
                if(r[0] > bestSol.bestRk)
                {
                    bestSol.x = sol.x;
                    bestSol.isFeasible = false;
                    bestSol.bestRk = sol.bestRk;
                    bestSol.w = sol.w;
                    bestSol.Z = sol.Z;
                    cout<<"SUMMARY SEARCH Best Solve After:"<<"rk = "<<bestSol.bestRk<<" objective = "<<bestSol.w<<endl;
                }
            }
            auto current_time = std::chrono::steady_clock::now();
            double elapsed_seconds = std::chrono::duration<double>(current_time - start_time).count();
            cout<<"ELAPSED TIME= "<<elapsed_seconds<<" "<<"LIMIT = "<<timeout_seconds<<endl;
            if (elapsed_seconds > timeout_seconds)
            {
                cout<<"TIMEOUT HAPPENED"<<endl;
                return bestSol;
            }else if(cntScenarios == formOptions.Z)
            {
                return bestSol;
            }
            formOptions.Z = formOptions.Z + min(z, cntScenarios - formOptions.Z);
            z = z * 2;
        }
    }
}