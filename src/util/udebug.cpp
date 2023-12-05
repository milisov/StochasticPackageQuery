#include "udebug.hpp"

#include <omp.h>
#include <fmt/core.h>
#include <cassert>

const char* RED = "\033[31m";
const char* GREEN = "\033[32m";
const char* RESET = "\033[0m";

Profiler::Profiler(){
}

void Profiler::clock(const string& label){
    timePoints[label] = Clock::now();
}

void Profiler::stop(const string& label){
    if (timePoints.count(label)){
        auto now = Clock::now();
        auto duration = std::chrono::duration<double, std::milli>(now - timePoints[label]);
        if (!clocks.count(label)) clocks[label] = {0.0, 0};
        clocks[label] = clocks[label] + make_pair(duration.count(), 1);
    }
}

void Profiler::add(const Profiler& pro){
    for (const auto& cl : pro.clocks){
        auto label = cl.first;
        if (!clocks.count(label)) clocks[label] = {0.0, 0};
        clocks[label] = clocks[label] + cl.second;
    }
}

void Profiler::print() const{
    for (const auto& cl : clocks){
        auto label = cl.first;
        if (!label.size()) label = "Ø";
        fmt::println("{}[count={} avg={}ms]", label, cl.second.second, cl.second.first/cl.second.second);
    }
}

void checkGurobi(const bool& error, GRBenv* env, GRBmodel* model, const char* file, const int& line){
	if (error){
		cerr << RED << "File " << file  \
			<< ", Line " << line << RESET << "\n";
		cerr << GRBgeterrormsg(env) << '\n';
		if (model) GRBfreemodel(model);
		if (env) GRBfreeenv(env);
		exit(1);
	}
}

string getGurobiStatus(const int& status){
    switch (status) {
        case GRB_OPTIMAL:
            return "Optimal solution found";
        case GRB_INF_OR_UNBD:
            return "Model is infeasible or unbounded";
        case GRB_INFEASIBLE:
            return "Model is infeasible";
        case GRB_UNBOUNDED:
            return "Model is unbounded";
        case GRB_CUTOFF:
            return "Model is solved to a cutoff";
        default:
            return "Unknown or undefined solution status";
    }
}

void getBasicVariables(GRBmodel* model, const int& numvars, vector<size_t>& basics){
    basics.reserve(numvars);
    vector<int> basisStatus (numvars, 1);
    assert(!GRBgetintattrarray(model, GRB_INT_ATTR_VBASIS, 0, numvars, basisStatus.data()));
    for (size_t i = 0; i < basisStatus.size(); ++i){
        if (basisStatus[i] == GRB_BASIC) basics.push_back(i);
    }
}