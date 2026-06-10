#include <iostream>
#include "solver/solversCore.hpp"
#include "CurveFit.hpp"
#include "spq/ssformulator.hpp"
#include "spq/rsformulator.hpp"
#include "core/checker.hpp"
#include <chrono>
#pragma once
using namespace std;

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

template <typename T>
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
    SolutionMetadata<T> bestSol;

    // Constructor
    SummarySearch(int M = 1e4,
                  std::shared_ptr<StochasticPackageQuery> spq = nullptr,
                  double epsilon = 1e-5, std::string validateTableNameBase = "") : data(Data::getInstance())
    {
        this->M = M;
        this->spq = spq;
        this->DB_optim = spq->tableName;
        this->DB_valid = validateTableNameBase;
        this->NTuples = pg.getTableSize(spq->tableName);
        this->cntScenarios = pg.getColumnLength(spq->tableName, "profit");
        this->probConstCnt = countProbConst(spq);
        this->epsilon = epsilon;

        spq->setTableName(DB_valid);
        this->checker = std::make_unique<SPQChecker>(this->spq);
        spq->setTableName(DB_optim);
    }

    // Extra methods specific to SummarySearch
    void guessOptimalConservativeness(std::vector<std::vector<std::pair<double, double>>> &history,
                                      std::vector<double> &alpha);

    void Best(std::shared_ptr<Objective> obj,
              std::vector<int> &x,
              History &H);

    SolutionMetadata<T> CSASolveBinSearch(GRBModel &model,
                                          std::vector<T> &x,
                                          Formulator &formulator,
                                          FormulateOptions &formOptions,
                                          SolveOptions &solveOptions);

    SolutionMetadata<T> CSASolveBinSearchStage3(Formulator &formulator,
                                                FormulateOptions &formOptions,
                                                SolveOptions &solveOptions);

    SolutionMetadata<T> CSASolveBinSearchRS(Formulator &formulator,
                                            FormulateOptions &formOptions, SolveOptions &solveOptions);

    SolutionMetadata<T> CSASolveBinSearchRSStage1(Formulator &formulator,
                                                  FormulateOptions &formOptions, SolveOptions &solveOptions);

    SolutionMetadata<T> summarySearch(shared_ptr<StochasticPackageQuery> spq, Formulator &formulator,
                                      FormulateOptions &formOptions,  SolveOptions &solveOptions, int z);

    SolutionMetadata<T> summarySearchRS(shared_ptr<StochasticPackageQuery> spq, Formulator &formulator,
                                        FormulateOptions &formOptions, SolveOptions &solveOptions);

    SolutionMetadata<T> summarySearchStage3(shared_ptr<StochasticPackageQuery> spq, Formulator &formulator, 
                                            FormulateOptions &formOptions, SolveOptions &solveOptions);
};

template <typename T>
SolutionMetadata<T> SummarySearch<T>::CSASolveBinSearchStage3(Formulator &formulator, FormulateOptions &formOptions,
                                                              SolveOptions &solveOptions)
{
    BinarySearchMetadata alpha_KMetadata(1e-5, 1.0, 0.0);
    std::vector<BinarySearchMetadata> history;
    initializeVector(history, probConstCnt, alpha_KMetadata);
    std::vector<double> alpha;
    initializeVector(alpha, probConstCnt, -1.0);
    int q = 0;
    int qAfterZequalsM = 0;

    int Z = formOptions.Z;
    bool partitioned = false;
    vector<T> x;

    while (true)
    {
        formOptions.iteration = q;
        // call the binary search function with for each constraint
        for (int i = 0; i < probConstCnt; i++)
        {
            double eps = 1e-3;
            if (formOptions.iteration > 0 && history[i].high - history[i].low > eps)
            {
                if (alpha[i] != -1.0)
                {
                    if (x.size() == 0)
                    {
                        //no solution make the system easier
                        history[i].high = alpha[i];
                    }
                    else
                    {
                        if (r[i] < 0)
                        {
                            //  the solution is infeasible -> make the system harder
                            history[i].low = alpha[i];
                        }
                        else
                        {
                            // the solution is feasible -> use less conservative summary
                            history[i].high = alpha[i];
                            // if(solveOptions.returnFirstFeasible)
                            // {
                            //     formOptions.alpha = alpha;
                            //     formOptions.innerConstraints = this->innerConstraints;
                            //     return bestSol;
                            // }
                        }
                    }
                }
                alpha[i] = history[i].low + (history[i].high - history[i].low) / 2;
                formOptions.alpha = alpha;
                deb(alpha);
                formOptions.innerConstraints = this->innerConstraints;
            }
            else if(formOptions.iteration != 0)
            {
                // cout << "BINARY SEARCH END CONDITION MET" << endl;
                return bestSol;
            }
        }
        if (partitioned == false && q > 0)
        {
            deb("Creating partitions");
            partitioned = true;
            formulator.createPartitions(spq, formOptions);
        }
        GRBModel model = formulator.formulate(spq, formOptions);
        if (x.size() > 0)
        {
            GRBVar *vars = model.getVars();
            int numVars = model.get(GRB_IntAttr_NumVars);
            for (int i = 0; i < formOptions.reducedIds.size(); i++)
            {
                int id = formOptions.reducedIds[i] - 1;
                try
                {
                    vars[i].set(GRB_DoubleAttr_Start, static_cast<double>(x[id]));
                }
                catch (const std::exception &e)
                {
                    std::cerr << e.what() << '\n';
                }
            }
            delete[] vars;
        }
        model.update();
        if (x.size() == 0)
        {
            initializeVector(x, NTuples, T(0));
        }
        // solveOptions.reduced = formOptions.reduced;
        // solveOptions.reducedIds = formOptions.reducedIds;
        // options.computeActiveness = false;
        solve(model, x, solveOptions);
        cout << "Iteration: " << q << " Z = " << Z << endl;
        //  validate() -> this will set the rk values, and calculate the Wq
        if (x.size() > 0)
        {
            validate(model, x, this->spq, solveOptions);
            bool isFeas = isFeasible(r);
            cout << (isFeas ? "Feasible" : "Infeasible") << endl;
            gpro.stop("effectiveRuntime");
            vector<double> feasibilities;
            vector<double> surpluses;
            SolType res;
            for (int i = 0; i < NTuples; i++)
            {
                if(x[i] > 0)
				{
					res[i + 1] = x[i];
				}
            }
            bool deterFeasible, probFeasible;
            bool feas = checker->feasible(res, feasibilities, surpluses, deterFeasible, probFeasible);
            double validObj = checker->getObjective(res);
            bestSol.solutions.push_back(make_tuple(gpro.getTime("effectiveRuntime"), satisfiedScenarios, r, this->W_q, deterFeasible, probFeasible, feasibilities, surpluses, validObj, formOptions.reducedIds.size(), formOptions.Z));
            gpro.clock("effectiveRuntime");
            int sense;
            if(formOptions.includeObjectiveFunction)
            {
                sense = model.get(GRB_IntAttr_ModelSense);
            }else
            {
                sense = GRB_MAXIMIZE;
            }
            bool currentIsBetter = isCurrentBetterThanBest(r, this->W_q, bestSol, sense);
            if (currentIsBetter || 1 == 1)
            {
                bestSol.x = x;
                bool isFeas = isFeasible(r);
                deb(isFeasible(r));
                bestSol.isFeasible = isFeas;
                bestSol.bestRk = r;
                bestSol.w = this->W_q;
                bestSol.Z = Z;
                bestSol.meanPerVaR = this->meanPerVaR;
                bestSol.variancePerVaR = this->variancePerVaR;
                bestSol.minInnerConstPerVaR = this->minInnerConstPerVaR;
            }
        }
        gpro.stop("effectiveRuntime");
        double elapsed_seconds = gpro.getTime("effectiveRuntime")/1000;
        gpro.clock("effectiveRuntime");
        deb(elapsed_seconds);
        if (elapsed_seconds > solveOptions.timeout_seconds)
        {
            cout << "TIME OUT HAPPENED" << endl;
            return bestSol;
        }

        if(qAfterZequalsM)
        {
            cout << "I SHOULD BE RUNNING NAIVE" <<endl;
            return bestSol;
        }

        q = q + 1;
        if (Z == cntScenarios)
        {
            qAfterZequalsM += 1;
        }
    }
}

template <typename T>
SolutionMetadata<T> SummarySearch<T>::CSASolveBinSearch(GRBModel &model, std::vector<T> &x, Formulator &formulator, FormulateOptions &formOptions,
                                                        SolveOptions &solveOptions)
{
    deb("solve csa");
    BinarySearchMetadata alpha_KMetadata(1e-5, 1.0, 0.0);
    std::vector<BinarySearchMetadata> history;
    initializeVector(history, probConstCnt, alpha_KMetadata);
    std::vector<double> alpha;
    initializeVector(alpha, probConstCnt, -1.0);
    int q = 0;
    int qAfterZequalsM = 0;

    int Z = formOptions.Z;
    bool partitioned = false;
    while (true)
    {
        cout << "Iteration: " << q << " Z = " << Z << endl;
        //  validate() -> this will set the rk values, and calculate the Wq
        if (x.size() > 0)
        {
            validate(model, x, this->spq, solveOptions);
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
            int sense;
            if(formOptions.includeObjectiveFunction)
            {
                sense = model.get(GRB_IntAttr_ModelSense);
            }else
            {
                sense = GRB_MAXIMIZE;
            }
            bool isCurrentBetter = isCurrentBetterThanBest(r, this->W_q, bestSol, sense);
            if (isCurrentBetter)
            {
                bestSol.x = x;
                bestSol.isFeasible = isFeasible(r);
                bestSol.bestRk = r;
                bestSol.w = this->W_q;
                bestSol.Z = Z;
            }
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
        if (partitioned == false)
        {
            deb("Creating partitions");
            partitioned = true;
            formulator.createPartitions(spq, formOptions);
        }
        GRBModel model = formulator.formulate(spq, formOptions);
        model.update();
        if (x.size() == 0)
        {
            initializeVector(x, NTuples, T(0));
        }
        solve(model, x, solveOptions);
    }
}

template <typename T>
SolutionMetadata<T> SummarySearch<T>::CSASolveBinSearchRS(Formulator &formulator, FormulateOptions &formOptions, SolveOptions &solveOptions)
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
    int Z = formOptions.Z;
    formulator.createPartitions(spq, formOptions);
    while (true)
    {
        cout << "Iteration: " << q << " Z = " << Z << " " << qAfterZequalsM << " " << cntScenarios << endl;

        if (Z == cntScenarios)
        {
            qAfterZequalsM += 1;
        }

        for (int i = 0; i < probConstCnt; i++)
        {
            alpha[i] = history[i].low + (history[i].high - history[i].low) / 2;
        }
        cout << "current alpha = " << alpha[0] << endl;
        formOptions.alpha = alpha;
        // deb(formOptions.alpha);
        formOptions.iteration = q;
        // formulate the I/LP
        GRBModel model = formulator.formulate(spq, formOptions);
        // optimize and store solution in x
        model.update();
        if (x.size() == 0)
        {
            initializeVector(x, NTuples, T(0));
        }

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
        solve(model, x, solveOptions);

        if (x.size() > 0)
        {
            validate(model, x, this->spq, solveOptions);
            formOptions.innerConstraints = this->innerConstraints;
            bool isFeas = isFeasible(r);
            int sense;
            if(formOptions.includeObjectiveFunction)
            {
                sense = model.get(GRB_IntAttr_ModelSense);
            }else
            {
                sense = GRB_MAXIMIZE;
            }
            bool isCurrentBetter = isCurrentBetterThanBest(r, this->W_q, bestSol, sense);
            if (isCurrentBetter)
            {
                bestSol.x = x;
                bestSol.isFeasible = isFeasible(r);
                deb(isFeasible(r));
                bestSol.bestRk = r;
                bestSol.bestActiveness = activeness;
                bestSol.bestActivenessNoAbsValue = activenessNoAbsValue;
                bestSol.w = this->W_q;
            }

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
SolutionMetadata<T> SummarySearch<T>::summarySearchRS(shared_ptr<StochasticPackageQuery> spq, Formulator &formulator,
                                                      FormulateOptions &formOptions, SolveOptions &solveOptions)
{
    std::vector<std::vector<std::vector<double>>> summaries;
    formOptions.M = this->M;
    std::vector<std::vector<std::pair<int, double>>> innerConstraintsDet = formOptions.innerConstraints;
    SolutionMetadata<T> sol;
    formOptions.innerConstraints = innerConstraintsDet;
    sol = CSASolveBinSearchRS(formulator, formOptions, solveOptions);

    return sol;
}

template <typename T>
SolutionMetadata<T> SummarySearch<T>::summarySearchStage3(shared_ptr<StochasticPackageQuery> spq, Formulator &formulator,
                                                          FormulateOptions &formOptions,SolveOptions &solveOptions)
{
    deb("Solving Stage 3 Summary Search");
    formOptions.M = this->M;
    SolutionMetadata<T> sol;
    sol = CSASolveBinSearchStage3(formulator, formOptions, solveOptions);
    deb(bestSol.bestRk);
    return sol;
}

template <typename T>
SolutionMetadata<T> SummarySearch<T>::summarySearch(shared_ptr<StochasticPackageQuery> spq, Formulator &formulator,
                                                    FormulateOptions &formOptions, SolveOptions &solveOptions, int z)
{
    formOptions.M = this->M;
    // formulate Deterministic ILP
    GRBModel model = formulator.formulate(spq, formOptions); // need to add the right values of formOptions from wherever you are calling summary search
    int numConstraints = model.get(GRB_IntAttr_NumConstrs);
    // printVariableNames(model);
    std::vector<T> x0;
    initializeVector(x0, NTuples, T(0));
    solve(model, x0, solveOptions);

    formOptions.Z = 1;
    while (true)
    {
        std::vector<T> x = x0; // copy solution to deterministic
        SolutionMetadata<T> sol;
        sol = CSASolveBinSearch(model, x, formulator, formOptions, solveOptions);
        if (sol.isFeasible && sol.epsilon <= this->epsilon)
        {
            return sol;
        }
        else
        {
            if (cntScenarios == formOptions.Z)
            {
                return bestSol;
            }
            formOptions.Z = formOptions.Z + min(z, cntScenarios - formOptions.Z);
            z = z * 2;
        }
    }
}