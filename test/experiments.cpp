#include <iostream>
#include <fstream>
#include "util/udebug.hpp"
#include "util/fileio.hpp"
#include "spq/spq.hpp"
#include "core/stat.hpp"
#include "spq/bounder.hpp"
#include "core/checker.hpp"
#include "solver/starter.hpp"
#include "solver/SummarySearch.hpp"
#include "solver/Naive.hpp"
#include "gurobi_c++.h"
#include <gurobi_c.h>
#include "spq/formulator.hpp"

using std::map;
using std::vector;

#include "solver/taylor.hpp"


void testNaive(string path, int M_input)
{
	int M = M_input;
	int M_hat = 1e6;
	double epsilon = 0.46;

	string filePath = path;
	auto spq = parseSpaqlFromFile(filePath);

	if (spq)
	{
		spq->validate();
		unique_ptr<Stat> stat = std::make_unique<Stat>();
		stat->analyze(spq->tableName);
		size_t N = 10000;
		double E = 50;
		// N is number of scenarios in order to approx, E - expected package size in sol
		// set values that are variables in the query
		Bounder bounder(spq, N, E);
		deb(spq->executable(), spq);
		std::vector<std::string> headers = {"Hardness", "Naive-objective", "Naive-Feas", "RuntimeNaive"};
		string output = "/home/fm2288/StochasticPackageQuery/test/Experiments/Naive" + path + ".csv";
		DataWriter writer(output, headers);

		Profiler stopwatchNaive;

		for (int h = -4; h <= -4; h = h + 1)
		{
			double NaiveObjective;
			bool NaiveFeas;
			bounder.set(h);
			deb(spq);
			Naive Nai(M, M_hat, spq);

			string labelNaive = "Naive with Hardness" + std::to_string(h);
			stopwatchNaive.clock(labelNaive);

			FormulateOptions formOptions;
			formOptions.reduced = false; 

			DecisionVarOptions decVarOptions;
			
			decVarOptions.lb = 0.0;
			decVarOptions.ub = 1.0;
			decVarOptions.obj = 0.0;
			decVarOptions.varType = GrbVarType::Binary;
			
			Solution solu = Nai.solveNaive(Nai.spq, 10, M, M_hat, formOptions, decVarOptions);
			stopwatchNaive.stop(labelNaive);
			double totalTimeNaive = stopwatchNaive.getTime(labelNaive);

			SPQChecker Check(Nai.spq);
			SolType resNaive;

			if (solu.x.size() > 0)
			{
				for (int i = 0; i < Nai.NTuples; i++)
				{
					resNaive[i + 1] = solu.x[i];
				}
				NaiveObjective = Check.getObjective(resNaive);
				NaiveFeas = Check.feasible(resNaive);
			}
			else
			{
				NaiveObjective = -1;
				NaiveFeas = 0;
			}
			//{"Hardness","Naive-objective","Naive-Feas","TimeNaive"}
			writer.addRow(h, NaiveObjective, NaiveFeas, totalTimeNaive);
		}
	}
}


void testSummarySearch(string path, int M_input)
{
	int M = M_input;
	int M_hat = 1e6;
	double epsilon = 0.46;

	map<string, Option> countConstraintOptions = {
		{"omit count constraint", false},
		{"scale down", false},
		{"scale factor", -1.0},
	};

	map<string, Option> curveFitOptions = {
		{"arctan", false},
		{"binarySearch", true}};

	string filePath = path;
	auto spq = parseSpaqlFromFile(filePath);

	if (spq)
	{
		spq->validate();
		unique_ptr<Stat> stat = std::make_unique<Stat>();
		stat->analyze(spq->tableName);
		size_t N = 10000;
		double E = 50;
		// N is number of scenarios in order to approx, E - expected package size in sol
		// set values that are variables in the query
		Bounder bounder(spq, N, E);
		deb(spq->executable(), spq);
		std::vector<std::string> headers = {"Hardness", "SS-objective", "SS-Feas", "Z", "RuntimeSummarySearch"};
		string output = "/home/fm2288/StochasticPackageQuery/test/Experiments/SS" + path + ".csv";
		DataWriter writer(output, headers);

		Profiler stopwatchSS;

		vector<int> reducedIds;
		bool reduced = false;
		for (int h = -4; h <= -4; h = h + 1)
		{
			double SSObjective;
			bool SSFeas;
			int Z;
			bounder.set(h);
			deb(spq);
			SummarySearch SS(M, M_hat, spq, epsilon);

			string labelSS = "SummarySearch with Hardness" + std::to_string(h);
			stopwatchSS.clock(labelSS);
			SolutionMetadata sol = SS.summarySearch(SS.spq, SS.M_hat, SS.M, 5, 1, reducedIds, reduced, countConstraintOptions, curveFitOptions);
			stopwatchSS.stop(labelSS);
			double totalTimeSS = stopwatchSS.getTime(labelSS);

			SPQChecker Check(SS.spq);
			if (sol.x.size() > 0)
			{
				deb(sol.x);
				SolType res;
				for (int i = 0; i < SS.NTuples; i++)
				{
					res[i + 1] = sol.x[i];
				}
				SSObjective = Check.getObjective(res);
				SSFeas = Check.feasible(res);
				Z = sol.Z;
			}
			else
			{
				SSObjective = -1;
				SSFeas = 0;
				Z = sol.Z;
			}
			//{"Hardness","SS-objective","SS-Feas","TimeSummarySearch"}
			writer.addRow(h, SSObjective, SSFeas, Z, totalTimeSS);
		}
	}
}

void testStochDualRed(string path, int M_input)
{
	map<string, Option> countConstraintOptions = {
		{"omit count constraint", true},
		{"scale down", false},
		{"scale factor", -1.0},
	};

	map<string, Option> curveFitOptions = {
		{"arctan", false},
		{"binarySearch", true}};

	int M = M_input;
	int M_hat = 1e6;
	double epsilon = 0.46;
	int qSz = 500;

	string filePath = path;
	// convert query to class
	auto spq = parseSpaqlFromFile(filePath);
	if (spq)
	{
		spq->validate();
		unique_ptr<Stat> stat = std::make_unique<Stat>();
		stat->analyze(spq->tableName);
		size_t N = 10000;
		double E = 50;

		Bounder bounder(spq, N, E);
		deb(spq->executable(), spq);
		std::vector<std::string> headers = {"hardness", "SDR-Objective", "SDR-feas", "SDR-optimal", "Z", "qSz", "Binary Search Steps", "Runtime"};
		string output = "/home/fm2288/StochasticPackageQuery/test/Experiments/SDR" + path + ".csv";
		DataWriter writer(output, headers);

		Profiler stopwatchSDR;

		vector<int> reducedIds;
		bool reduced = false;
		for (int h = -4; h <= 4; h = h + 1)
		{
			bounder.set(h);
			deb(spq);
			SummarySearch SS(M, M_hat, spq, epsilon);

			string labelSS = "SummarySearch with Hardness" + std::to_string(h);
			stopwatchSDR.clock(labelSS);
			SolutionMetadata sol = SS.stochDualReducer(SS.spq, qSz, countConstraintOptions, curveFitOptions);
			stopwatchSDR.stop(labelSS);
			double totalTimeSS = stopwatchSDR.getTime(labelSS);

			if (sol.x.size() > 0)
			{
				deb(sol.x);

				SPQChecker Check(SS.spq);
				SolType res;
				for (int i = 0; i < SS.NTuples; i++)
				{
					res[i + 1] = sol.x[i];
				}

				bool feasible = Check.feasible(res);
				double obj = Check.getObjective(res);

				cout << (feasible ? "Validation Feasible Solution" : "Validation Infeasible Solution") << endl;
				cout << "Training Objective Summary Search = " << sol.w << endl;
				cout << "Validation Objective Summary Search = " << obj << endl;

				writer.addRow(h, obj, feasible, sol.isOptimal, sol.Z, sol.qSz, sol.binarySearchSteps, totalTimeSS);
			}
		}
	}
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <N> <M> <algorithm>" << std::endl;
        return 1;
    }

    int N = std::stoi(argv[1]);
    int M = std::stoi(argv[2]);
    std::string algorithm = argv[3];

    std::string dbPath = fmt::format("resource/sqls/_stocks_{}_{}.spaql", N, M);

    std::cout << "File Path: " << dbPath << std::endl;
    std::cout << "Algorithm: " << algorithm << std::endl;

	if(algorithm == "Naive")
	{
		testNaive(dbPath, M);
	}else
	if(algorithm == "SummarySearch")
	{
		testSummarySearch(dbPath, M);
	}else
	if(algorithm == "SDR")
	{
		testStochDualRed(dbPath, M);
	}else
	{
		cout<<"Please enter a correct algorithm name next time"<<endl;
	}

    return 0;
}