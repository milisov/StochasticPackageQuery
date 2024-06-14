#include "starter.hpp"

using std::make_unique;

map<string, Option> starterOptions = {
    {"number_of_cores", static_cast<int>(Config::getInstance()->nPhysicalCores)}
};

Starter::Starter(shared_ptr<StochasticPackageQuery> spq, const vector<long long>& ids, const map<string, Option>& options): spq(spq), options(options){
    stat = make_unique<Stat>();
    for (const auto& p : starterOptions){
        if (!options.count(p.first)) this->options[p.first] = p.second;
    }
    nCores = boost::get<int>(this->options["number_of_cores"]);
    auto tableSize = stat->pg->getTableSize(spq->tableName);
    for (const auto& id : ids){
        if (id > tableSize || id < 1){
            cerr << fmt::format("Cannot process id '{}' in table '{}'\n", id, spq->tableName);
            exit(1);
        }
    }
    idx = make_unique<UniqueIndexer>(nCores, tableSize, ids);
}

void Starter::solveCurrentSystem(SolIndType& currentSol){
    size_t nInds = idx->size();
    size_t numvars = nInds + spq->getCvarCount() + spq->isStochasticObjective();
    size_t m = spq->cons.size() + spq->isStochasticObjective();
    size_t nzN = sol.size();
    vector<int> cbeg, cind;
    vector<char> sense; sense.reserve(2*m);
    vector<double> cval, rhs;
    cbeg.reserve(2*m+1); cind.reserve(2*m*numvars); cval.reserve(2*m*numvars); rhs.reserve(2*m);
    cbeg.push_back(0);

    shared_ptr<AttrConstraint> attrCon;
    shared_ptr<BoundConstraint> boundCon;
    shared_ptr<ProbConstraint> probCon;
    for (const auto& con : spq->cons){
        if (isStochastic(con, probCon, attrCon) && attrCon){
            double v = spq->getValue(probCon->v);
            double p = spq->getValue(probCon->p);
            vector<double> stoCon (nInds);
            double stoBound;
            Inequality stoIneq = Inequality::lteq;
            if (getVar(con)){
                
            }
        }
    }
}

void Starter::solve(){
    map<string, pair<double, double>> intervals;
    shared_ptr<ProbConstraint> probCon;
    for (const auto& con : spq->cons){
        if (isStochastic(con, probCon)){
            double p1 = spq->getValue(probCon->p);
            double p2 = 0;
            if (getVar(con)){
                if (probCon->psign == Inequality::lteq) p2 = 1;
            }
            if (getCvar(con)){
                // To be added
            }
            intervals[boost::get<string>(probCon->p)] = {p1, p2};
        }
    }
    SolIndType currentSol;
    solveCurrentSystem(currentSol);
}

