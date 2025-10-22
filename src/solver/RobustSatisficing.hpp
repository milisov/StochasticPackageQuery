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
    std::vector<pair<int, double>> bestPosActivenessRS;
    std::vector<pair<int, double>> bestNegActivenessRS;
    RobustSatisficing(int M = 1e4,
                      std::shared_ptr<StochasticPackageQuery> spq = nullptr,
                      double epsilon = 1e-5)
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
    SolutionMetadata<int> stochasticDualReducer(std::shared_ptr<StochasticPackageQuery> spq, std::map<std::string, Option> &curveFitOptions);

    template <typename T>
    SolutionMetadata<T> solveSAA(GRBModel &model, FormulateOptions &formOptions)
    {
        vector<T>x;
        initializeVector(x,NTuples,T(0));
        SolveOptions options; 
        options.reduced = formOptions.reduced;
        options.reducedIds = formOptions.reducedIds;
        options.computeActiveness = formOptions.computeActiveness;
        solve(model,x, options);
        SolutionMetadata<T> sol; 
        if(x.size()>0)
        {
            validate(model, x, spq, options);
            if(isFeasible(r))
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

    vector<int> reduceTuplesStageNoObjCons(std::shared_ptr<StochasticPackageQuery> spq,
                                                         FormulateOptions &formOptions,
                                                         std::map<std::string, Option> &curveFitOptions);

    vector<int> reduceTuplesStage(std::shared_ptr<StochasticPackageQuery> spq,
                                  FormulateOptions &formOptions,
                                  std::map<std::string, Option> &curveFitOptions,double Z0);

    vector<int> finalReduce(map<int, double> &reducedIdsMap, vector<int> &reducedIds,int q);

    double findBestObjectiveStage(std::shared_ptr<StochasticPackageQuery> spq,
                                  FormulateOptions &formOptions,
                                  std::map<std::string, Option> &curveFitOptions, double Z0);


    double findBestObjectiveStageSAA(GRBModel &model,std::shared_ptr<StochasticPackageQuery> spq,
                                  RSFormulator &formulator,
                                  FormulateOptions &formOptions,
                                  std::map<std::string, Option> &curveFitOptions, double Z0);

    
    SolutionMetadata<int> finalStageILP(GRBModel &model,std::shared_ptr<StochasticPackageQuery> spq,
                                  RSFormulator &formulator,
                                  FormulateOptions &formOptions, double bestEps,double Z0);

                                
};