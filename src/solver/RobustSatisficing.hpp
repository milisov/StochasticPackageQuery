#include "spq/rsformulator.hpp"
#include "solver/solversCore.hpp"
#include "solver/SummarySearch.hpp"
#include "solver/Naive.hpp"
#include "util/data.hpp"
#include <iostream>

class RobustSatisficing : public Solver
{
public:
    double epsilon;
    SolutionMetadata<int> bestSolGlobal;
    std::vector<std::vector<pair<int, double>>> bestPosActiveness;
    std::vector<std::vector<pair<int, double>>> bestNegActiveness;
    vector<pair<int, double>> rankingNonZero;
    vector<std::tuple<int, double, double>> rankingZero;
    RobustSatisficing(int M = 1e4,
                      std::shared_ptr<StochasticPackageQuery> spq = nullptr,
                      double epsilon = 1e-5, string validateTableNameBase = "")
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

    SolutionMetadata<int> solveDeterministic(std::shared_ptr<StochasticPackageQuery> spq, SolveOptions &solveOptions);

    SolutionMetadata<int> stochasticDualReducer(std::shared_ptr<StochasticPackageQuery> spq, SolveOptions &solveOptions, map<string, bool> ablate = {{"tuples", false}, {"scenarios", false}, {"ablate", false}});

    template <typename T>
    SolutionMetadata<T> solveSAA(GRBModel &model, FormulateOptions &formOptions, SolveOptions &solveOptions)
    {
        vector<T> x;
        initializeVector(x, NTuples, T(0));
        SolveOptions options;
        options.reduced = formOptions.reduced;
        options.reducedIds = formOptions.reducedIds;
        // options.computeActiveness = formOptions.computeActiveness;
        solve(model, x, options);
        SolutionMetadata<T> sol;
        if (x.size() > 0)
        {
            validate(model, x, spq, options);
            if (isFeasible(r))
            {
                sol.x = x;
                sol.isFeasible = true;
                return sol;
            }
        }
        sol.x = x;
        sol.isFeasible = false;
        return sol;
    }

    void populateMapNonZero(map<int, double> &reducedIds, const vector<double> &solDet);
    void populateMapFromVector(map<int, double> &reducedIdsMap, const vector<int> &reduced);
    void getReduced(vector<int> &reducedIds, vector<double> &solDet, vector<double> &solStage1, int qTarget);
    void getReducedFromRS(vector<int> &reducedIds, vector<double> &solStage1, int qTarget);

    void randomScenarioSelection(int Z, int conOrder, std::set<int> &selectedScenarios, unsigned int seed = 42);

    void selectTuplesFromRanking(vector<int> &reducedIds, int qTarget,
                                  vector<pair<int, double>> &rankingNonzero,
                                  vector<std::tuple<int, double, double>> &rankingZero);

    // Warmstart helper functions
    void applyWarmstart(GRBModel &model, const std::vector<int> &vbasis, const std::vector<int> &cbasis);
    void saveWarmstart(GRBModel &model, std::vector<int> &vbasis, std::vector<int> &cbasis);

    vector<int> reduceTuplesStageNoObjCons(std::shared_ptr<StochasticPackageQuery> spq,
                                           FormulateOptions &formOptions,
                                           SolveOptions &solveOptions, vector<double> &solDet, int qTarget);

    void reduceTuplesStageNoObjConsNoObj(std::shared_ptr<StochasticPackageQuery> spq,
                                        FormulateOptions &formOptions,
                                        SolveOptions &solveOptions, vector<double> &solDet, int qTarget);

    void reduceTuples(std::shared_ptr<StochasticPackageQuery> spq,
                                  FormulateOptions &formOptions,
                                  SolveOptions &solveOptions, vector<double> &solDet, int qTarget);

    void reduceScenarios(std::shared_ptr<StochasticPackageQuery> spq,
                                    FormulateOptions &formOptions,
                                    SolveOptions &solveOptions, int qTarget);

    void reduceTuplesLCVaR(std::shared_ptr<StochasticPackageQuery> spq,
                           FormulateOptions &formOptions,
                           SolveOptions &solveOptions, vector<double> &solDet, int qTarget);

    void reduceTuplesRobustSatisficingRelaxLP(std::shared_ptr<StochasticPackageQuery> spq,
                                   FormulateOptions &formOptions,
                                   SolveOptions &solveOptions, vector<double> &solDet, int qTarget);

    vector<int> reduceTuplesStageNoObjConsUpdatingBounds(std::shared_ptr<StochasticPackageQuery> spq,
                                                         FormulateOptions &formOptions,
                                                         std::map<std::string, Option> &curveFitOptions, int q);

    vector<int> reduceTuplesStage(std::shared_ptr<StochasticPackageQuery> spq,
                                  FormulateOptions &formOptions,
                                  std::map<std::string, Option> &curveFitOptions, double Z0);

    vector<int> finalReduce(map<int, double> &reducedIdsMap, vector<int> &reducedIds, int q, SolveOptions &solveOptions);

    vector<vector<pair<int, double>>> findBestObjectiveStage(std::shared_ptr<StochasticPackageQuery> spq,
                                  FormulateOptions &formOptions,
                                  SolveOptions &solveOptions, double Z0);

    SolutionMetadata<int> runFinalStage(std::shared_ptr<StochasticPackageQuery> spq,
                                        FormulateOptions &formOptions,
                                        SolveOptions &solveOptions,
                                        vector<int> &finalReducedIds,
                                        int qTarget);

    
    SolutionMetadata<int> runFinalStageAlternateVaRBound(std::shared_ptr<StochasticPackageQuery> spq,
                                        FormulateOptions &formOptions,
                                        SolveOptions &solveOptions,
                                        vector<int> &finalReducedIds,
                                        int qTarget);

    SolutionMetadata<int> runFinalStageVarianceControl(std::shared_ptr<StochasticPackageQuery> spq,
                                        FormulateOptions &formOptions,
                                        SolveOptions &solveOptions,
                                        vector<int> &finalReducedIds,
                                        int qTarget);

    bool finalStageTerminate(SolutionMetadata<int> &sol1, SolutionMetadata<int> &sol2, double bound, double threshold, int sense);

    SolutionMetadata<int> runNaiveVarianceControl(std::shared_ptr<StochasticPackageQuery> spq,
                                        FormulateOptions &formOptions,
                                        SolveOptions &solveOptions,
                                        int qTarget);

    
    SolutionMetadata<int> runNaiveFinalStage(std::shared_ptr<StochasticPackageQuery> spq,
                                        FormulateOptions &formOptions,
                                        SolveOptions &solveOptions,
                                        int qTarget);

    double findBestObjectiveStageSAA(GRBModel &model, std::shared_ptr<StochasticPackageQuery> spq,
                                     RSFormulator &formulator,
                                     FormulateOptions &formOptions,
                                     std::map<std::string, Option> &curveFitOptions, double Z0);

    SolutionMetadata<int> finalStageILP(GRBModel &model, std::shared_ptr<StochasticPackageQuery> spq,
                                        RSFormulator &formulator,
                                        FormulateOptions &formOptions, double bestEps, double Z0);

    // Compute variance control bound: x^T * mu - sqrt(p/(1-p)) * ||x .* sigma||_2
    template <typename T>
    double computeVarianceControlBound(const std::vector<T> &x, double p)
    {
        const std::vector<double> &mu = data.stockExpectedProfit;
        const std::vector<double> &sigma = data.stockProfitStdDev;

        double k = std::sqrt(p / (1.0 - p));

        double expectedProfit = 0.0;
        double normSquared = 0.0;

        for (size_t i = 0; i < x.size(); i++)
        {
            double xi = static_cast<double>(x[i]);
            expectedProfit += mu[i] * xi;
            double yi = sigma[i] * xi;
            normSquared += yi * yi;
        }

        deb(expectedProfit, k, normSquared);

        double norm = std::sqrt(normSquared);
        deb(norm);
        return expectedProfit - k * norm;
    }
};