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

using std::map;
using std::vector;

void analyzeAll()
{
	std::vector<int> nStocks = {3, 4, 5, 6};
	std::vector<int> nPaths = {2, 4};
	for (auto nStock : nStocks)
	{
		for (auto nPath : nPaths)
		{
			INIT(pro);
			CLOCK("a");
			string filePath = fmt::format("resource/sqls/_stocks_{}_{}.spaql", nStock, nPath);
			auto spq = parseSpaqlFromFile(filePath);
			if (spq)
			{
				cout << "Success!\n"
					 << spq;
				deb(spq->validate());
				unique_ptr<Stat> stat = std::make_unique<Stat>();
				stat->analyze(spq->tableName);
			}
			STOP("a");
			PRINT(pro);
		}
	}
}

#include "solver/taylor.hpp"

// void experimentEpsilonvsObj(double h)
// {
// 	string filePath = fmt::format("resource/sqls/_stocks_3_2.spaql");
// 	// convert query to class
// 	auto spq = parseSpaqlFromFile(filePath);
// 	if (spq)
// 	{
// 		spq->validate();
// 		unique_ptr<Stat> stat = std::make_unique<Stat>();
// 		stat->analyze(spq->tableName);
// 		size_t N = 10000;
// 		double E = 50;
// 		// N is number of scenarios in order to approx, E - expected package size in sol
// 		// set values that are variables in the query
// 		Bounder bounder(spq, N, E);
// 		deb(spq->executable(), spq);
// 		std::vector<std::string> headers = {"Hardness", "Epsilon", "SummarySearchObj", "SS-Feas"};
// 		DataWriter writer("/home/fm2288/StochasticPackageQuery/test/epsilonVsObjective.csv", headers);

// 		bounder.set(h);

// 		for (double eps = 0.15; eps >= 0.025; eps = eps - 0.025)
// 		{
// 			cout << eps << endl;
// 			SummarySearch SS(100, 1e6, spq, eps);

// 			SolutionMetadata sol = SS.summarySearch(SS.spq, SS.M_hat, SS.M, 5, 1, eps);

// 			cout << "Summary Search Produces" << endl;
// 			cout << (sol.isFeasible ? "Feasible Solution" : "Infeasible Solution") << endl;
// 			deb(sol.x);

// 			SPQChecker Check(SS.spq);
// 			SolType res;
// 			for (int i = 0; i < SS.NTuples; i++)
// 			{
// 				res[i + 1] = sol.x[i];
// 			}

// 			cout << (Check.feasible(res) ? "Validation Feasible Solution" : "Validation Infeasible Solution") << endl;
// 			cout << "Training Objective Summary Search = " << sol.w << endl;
// 			cout << "Validation Objective Summary Search = " << Check.getObjective(res) << endl;

// 			writer.addRow(h, eps, Check.getObjective(res), Check.feasible(res));
// 		}
// 	}
// }

void testNaive()
{
	std::ofstream outputFile("/home/fm2288/StochasticPackageQuery/test/output.txt");
	string filePath = fmt::format("resource/sqls/_stocks_3_2.spaql");
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

		vector<pair<SolType, int>> solVectors;

		for (int h = 4; h >= -4; h = h - 1)
		{
			bounder.set(h);

			Naive Nai(100, 100, spq);
			Solution solu = Nai.solveNaive(Nai.spq, 100, 100, 100);
			SPQChecker Check(spq);
			SolType resNaive;

			for (int i = 0; i < Nai.NTuples; i++)
			{
				resNaive[i + 1] = solu.x[i];
			}

			if (solu.isFeasible)
			{
				deb(solu.x);
				cout << "Solution is Feasible" << endl;
				cout << "Training Objective Naive  = " << solu.W_q << endl;
			}
			else
			{
				cout << "Solution is Infeasible" << endl;
			}

			cout << (Check.feasible(resNaive) ? "Validation Feasible Solution" : "Validation Infeasible Solution") << endl;
			cout << "Training Objective Naive  = " << Check.getObjective(resNaive) << endl;

			solVectors.push_back(make_pair(resNaive, h));
			for (int i = 0; i < solVectors.size(); i++)
			{
				outputFile << "Solution for hardness = " << solVectors[i].second
						   << " returns " << (Check.feasible(solVectors[i].first) ? "Validation Feasible Solution" : "Validation Infeasible Solution")
						   << " for hardness level h = " << h << std::endl;
			}
		}
	}
	outputFile.close();
}

void debugNaive()
{
	string filePath = fmt::format("resource/sqls/_stocks_3_2.spaql");
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

		vector<pair<SolType, int>> solVectors;

		bounder.set(1);

		Naive Nai(100, 100, spq);
		Solution solu = Nai.solveNaive(Nai.spq, 100, 100, 100);
		SPQChecker Check(spq);
		SolType resNaive;

		for (int i = 0; i < Nai.NTuples; i++)
		{
			resNaive[i + 1] = solu.x[i];
		}

		if (solu.isFeasible)
		{
			deb(solu.x);
			cout << "Solution is Feasible" << endl;
			cout << "Training Objective Naive  = " << solu.W_q << endl;
		}
		else
		{
			cout << "Solution is Infeasible" << endl;
		}

		cout << (Check.feasible(resNaive) ? "Validation Feasible Solution" : "Validation Infeasible Solution") << endl;
		cout << "Training Objective Naive  = " << Check.getObjective(resNaive) << endl;
	}
}


void testSSvsNaive()
{
	int M = 10000;
	int M_hat = 1e6;
	double epsilon = 0.46;

	map<string, Option> countConstraintOptions = {
    {"omit count constraint", false}, 
    {"scale down", false}, 
    {"scale factor", -1.0},
	};

	map<string, Option> curveFitOptions = {
    {"arctan", false}, 
    {"binarySearch", true} 
	};

	string filePath = fmt::format("resource/sqls/_stocks_4_4.spaql");
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
		std::vector<std::string> headers = {"Hardness","SS-objective","SS-Feas","Z","RuntimeSummarySearch","Naive-objective","Naive-Feas","RuntimeNaive"};
		DataWriter writer("/home/fm2288/StochasticPackageQuery/test/BinarySearchCurveFit/SSvsNaive_4_4_eps046.csv", headers);

		Profiler stopwatchSS;
		Profiler stopwatchNaive;

		vector<int>reducedIds;
		bool reduced = false;
		for (int h = -4; h <= 4; h = h + 1)
		{
			double SSObjective;
			bool SSFeas;
			double NaiveObjective; 
			bool NaiveFeas;
			int Z;
			bounder.set(h);
			deb(spq);
			SummarySearch SS(M, M_hat, spq, epsilon);

			string labelSS = "SummarySearch with Hardness" + std::to_string(h);
			stopwatchSS.clock(labelSS);
			SolutionMetadata sol = SS.summarySearch(SS.spq, SS.M_hat, SS.M, 5, 1, reducedIds, reduced, countConstraintOptions,curveFitOptions);
			stopwatchSS.stop(labelSS);
			double totalTimeSS = stopwatchSS.getTime(labelSS);

			cout << "Summary Search Produces" << endl;
			cout << (sol.isFeasible ? "Feasible Solution" : "Infeasible Solution") << endl;
			deb(sol.x);

			SPQChecker Check(SS.spq);
			if(sol.x.size() > 0)
			{
				deb(sol.x);
				SolType res;
				for (int i = 0; i < SS.NTuples; i++)
				{
					res[i + 1] = sol.x[i];
				}

				cout << (Check.feasible(res) ? "Validation Feasible Solution" : "Validation Infeasible Solution") << endl;
				cout << "Training Objective Summary Search = " << sol.w << endl;
				cout << "Validation Objective Summary Search = " << Check.getObjective(res) << endl;

				SSObjective = Check.getObjective(res);
				SSFeas = Check.feasible(res);
				Z = sol.Z;

			}else
			{
				SSObjective = -1;
				SSFeas = 0;	
				Z = sol.Z;
			}

			Naive Nai(M, M_hat, spq);

			string labelNaive = "Naive with Hardness" + std::to_string(h);
			stopwatchNaive.clock(labelNaive);
			Solution solu = Nai.solveNaive(Nai.spq, 10, M, M_hat);
			stopwatchNaive.stop(labelNaive);
			double totalTimeNaive = stopwatchNaive.getTime(labelNaive);

			SolType resNaive;

			if(solu.x.size() > 0)
			{
				for (int i = 0; i < Nai.NTuples; i++)
				{
					resNaive[i + 1] = solu.x[i];
				}
				if (solu.isFeasible)
				{
					deb(solu.x);
					cout << "Solution is Feasible" << endl;
					cout << "Training Objective Naive  = " << solu.W_q << endl;
				}
				else
				{
					cout << "Solution is Infeasible" << endl;
				}

				cout << (Check.feasible(resNaive) ? "Validation Feasible Solution" : "Validation Infeasible Solution") << endl;
				cout << "Training Objective Naive  = " << Check.getObjective(resNaive) << endl;

				NaiveObjective = Check.getObjective(resNaive); 
				NaiveFeas = Check.feasible(resNaive);
			}else
			{
				NaiveObjective = -1; 
				NaiveFeas = 0;
			}
			//{"Hardness","SS-objective","SS-Feas","TimeSummarySearch","Naive-objective","Naive-Feas","TimeNaive"}
			writer.addRow(h, SSObjective, SSFeas, Z,totalTimeSS, NaiveObjective, NaiveFeas, totalTimeNaive);
		}
	}
}

void testStochDualRed()
{
	map<string, Option> countConstraintOptions = {
    {"omit count constraint", true}, 
    {"scale down", false}, 
    {"scale factor", -1.0},
	};

	map<string, Option> curveFitOptions = {
    {"arctan", false}, 
    {"binarySearch", true} 
	};

	int M = 100;
	int M_hat = 1e6;
	double epsilon = 0.46;
	int qSz = 500; 

	string filePath = fmt::format("resource/sqls/_stocks_5_2.spaql");
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
		std::vector<std::string> headers = {"hardness", "SDR-Objective", "SDR-feas","Z", "Binary Search Steps", "Runtime"};
		DataWriter writer("/home/fm2288/StochasticPackageQuery/test/StochDualReducerNoCountConstBinarySearchFit/stochdualred_5_2_eps046_qsz500.csv", headers);

		Profiler stopwatchSDR;


		vector<int>reducedIds;
		bool reduced = false;
		for (int h = -4; h <= 4; h = h + 1)
		{
			bounder.set(h);
			deb(spq);
			SummarySearch SS(M, M_hat, spq, epsilon);

			string labelSS = "SummarySearch with Hardness" + std::to_string(h);
			stopwatchSDR.clock(labelSS);
			SolutionMetadata sol = SS.stochDualReducer(SS.spq,qSz,countConstraintOptions,curveFitOptions);
			stopwatchSDR.stop(labelSS);
			double totalTimeSS = stopwatchSDR.getTime(labelSS);

			cout << "Summary Search Produces" << endl;
			cout << (sol.isFeasible ? "Feasible Solution" : "Infeasible Solution") << endl;
			
			if(sol.x.size() > 0)
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

				writer.addRow(h,obj,feasible,sol.Z,sol.binarySearchSteps, totalTimeSS);
			}else
			{
				writer.addRow(h, -1, 0, sol.Z,sol.binarySearchSteps, totalTimeSS);
			}
		}
	}
}



//
// void testInit()
// {
// 	string filePath = fmt::format("resource/sqls/_stocks_4_4.spaql");
// 	// convert query to class
// 	auto spq = parseSpaqlFromFile(filePath);
// 	if (spq)
// 	{
// 		spq->validate();
// 		unique_ptr<Stat> stat = std::make_unique<Stat>();
// 		stat->analyze(spq->tableName);
// 		size_t N = 10000;
// 		double E = 50;
// 		// N is number of scenarios in order to approx, E - expected package size in sol
// 		// set values that are variables in the query
// 		Bounder bounder(spq, N, E);
// 		deb(spq->executable(), spq);
// 		std::vector<std::string> headers = {"Hardness", "SummarySearchObj", "NaiveObj", "SS-Feas", "Naive-Feas", "TimeSummarySearch", "TimeNaive"};
// 		DataWriter writer("/home/fm2288/StochasticPackageQuery/test/dataStocks_4_4_eps050.csv", headers);

// 		Profiler stopwatchSS;
// 		Profiler stopwatchNaive;

// 		for (int h = -4; h <= -3; h = h + 1)
// 		{
// 			bounder.set(h);
// 			deb(spq);
// 			SummarySearch SS(10000, 1e6, spq, 0.46);

// 			string labelSS = "SummarySearch with Hardness" + std::to_string(h);
// 			stopwatchSS.clock(labelSS);
// 			SolutionMetadata sol = SS.summarySearch(SS.spq, SS.M_hat, SS.M, 5, 1, SS.epsilon);
// 			stopwatchSS.stop(labelSS);
// 			double totalTimeSS = stopwatchSS.getTime(labelSS);

// 			cout << "Summary Search Produces" << endl;
// 			cout << (sol.isFeasible ? "Feasible Solution" : "Infeasible Solution") << endl;
// 			deb(sol.x);

// 			SPQChecker Check(SS.spq);
// 			SolType res;
// 			for (int i = 0; i < SS.NTuples; i++)
// 			{
// 				res[i + 1] = sol.x[i];
// 			}

// 			cout << (Check.feasible(res) ? "Validation Feasible Solution" : "Validation Infeasible Solution") << endl;
// 			cout << "Training Objective Summary Search = " << sol.w << endl;
// 			cout << "Validation Objective Summary Search = " << Check.getObjective(res) << endl;

// 			Naive Nai(10000, 10000, spq);

// 			string labelNaive = "Naive with Hardness" + std::to_string(h);
// 			stopwatchNaive.clock(labelNaive);
// 			Solution solu = Nai.solveNaive(Nai.spq, 100, 10000, 10000);
// 			stopwatchNaive.stop(labelNaive);
// 			double totalTimeNaive = stopwatchNaive.getTime(labelNaive);

// 			SolType resNaive;

// 			for (int i = 0; i < Nai.NTuples; i++)
// 			{
// 				resNaive[i + 1] = solu.x[i];
// 			}
// 			if (solu.isFeasible)
// 			{
// 				deb(solu.x);
// 				cout << "Solution is Feasible" << endl;
// 				cout << "Training Objective Naive  = " << solu.W_q << endl;
// 			}
// 			else
// 			{
// 				cout << "Solution is Infeasible" << endl;
// 			}

// 			cout << (Check.feasible(resNaive) ? "Validation Feasible Solution" : "Validation Infeasible Solution") << endl;
// 			cout << "Training Objective Naive  = " << Check.getObjective(resNaive) << endl;
// 			writer.addRow(h, Check.getObjective(res), Check.getObjective(resNaive), Check.feasible(res), Check.feasible(resNaive), totalTimeSS, totalTimeNaive);
// 		}
// 	}
// }


// void testInit2()
// {
// 	string filePath = fmt::format("resource/sqls/_stocks_3_2.spaql");
// 	// convert query to class
// 	auto spq = parseSpaqlFromFile(filePath);
// 	if (spq)
// 	{
// 		spq->validate();
// 		unique_ptr<Stat> stat = std::make_unique<Stat>();
// 		stat->analyze(spq->tableName);
// 		size_t N = 10000;
// 		double E = 50;
// 		// N is number of scenarios in order to approx, E - expected package size in sol
// 		// set values that are variables in the query
// 		Bounder bounder(spq, N, E);
// 		deb(spq->executable(), spq);
// 		std::vector<std::string> headers = {"Hardness", "SummarySearchObj", "NaiveObj", "SS-Feas", "Naive-Feas", "TimeSummarySearch", "TimeNaive"};
// 		DataWriter writer("/home/fm2288/StochasticPackageQuery/test/dataStocks_4_4_eps046.csv", headers);

// 		for (int h = -4; h <= 3; h = h + 1)
// 		{
// 			double summarySearchObj;
// 			double naiveObj;
// 			bool summarySearchFeas;
// 			bool naiveFeas;
// 			bounder.set(h);
// 			deb(spq);
// 			SummarySearch SS(100, 1e6, spq, 0.46);

// 			SolutionMetadata sol = SS.summarySearch(SS.spq, SS.M_hat, SS.M, 5, 1, 0.46);

// 			cout << "Summary Search Produces" << endl;
// 			cout << (sol.isFeasible ? "Feasible Solution" : "Infeasible Solution") << endl;
// 			deb(sol.x);
// 			SPQChecker Check(SS.spq);
			
// 			if(sol.isFeasible)
// 			{
// 				SolType res;
// 				for (int i = 0; i < SS.NTuples; i++)
// 				{
// 					res[i + 1] = sol.x[i];
// 				}
// 				summarySearchObj = Check.getObjective(res);
// 				summarySearchFeas = Check.feasible(res);
// 				cout << (summarySearchFeas ? "Validation Feasible Solution" : "Validation Infeasible Solution") << endl;
// 				cout << "Training Objective Summary Search = " << sol.w << endl;
// 				cout << "Validation Objective Summary Search = " << summarySearchObj << endl;
// 			}

// 			Naive Nai(100, 100, spq);

// 			Solution solu = Nai.solveNaive(Nai.spq, 100, 100, 100);

// 			SolType resNaive;

// 			for (int i = 0; i < Nai.NTuples; i++)
// 			{
// 				resNaive[i + 1] = solu.x[i];
// 			}
// 			if (solu.isFeasible)
// 			{
// 				deb(solu.x);
// 				cout << "Solution is Feasible" << endl;
// 				cout << "Training Objective Naive  = " << solu.W_q << endl;
// 			}
// 			else
// 			{
// 				cout << "Solution is Infeasible" << endl;
// 			}

// 			naiveObj = Check.getObjective(resNaive);
// 			naiveFeas = Check.feasible(resNaive);
// 			cout << (naiveFeas ? "Validation Feasible Solution" : "Validation Infeasible Solution") << endl;
// 			cout << "Training Objective Naive  = " << naiveObj << endl;
// 			writer.addRow(h, summarySearchObj, naiveObj, summarySearchFeas, naiveFeas, sol.Runtime, solu.Runtime);
// 		}
// 	}
// }




// class summarysearch input spq, options (arctan, different epsilon)

// void test(){
// 	INIT(pro);
// 	CLOCK("a");
// 	for (int i = 6; i <= 6; i ++){
// 		string filePath = fmt::format("resource/sqls/_stocks_{}_2.spaql", i);
// 		auto spq = parseSpaqlFromFile(filePath);
// 		if (spq){
// 			// cout << "Success!\n" << spq;
// 			// deb(spq->validate());
// 			spq->validate();
// 			unique_ptr<Stat> stat = std::make_unique<Stat>();
// 			stat->analyze(spq->tableName);
// 			size_t N = 10000;
// 			double E = 50;
// 			CLK(pro, "initBounder");
// 			Bounder bounder (spq, N, E);
// 			STP(pro, "initBounder");
// 			vector<double> hards;
// 			for (double i = -10; i <= 10; i ++) hards.push_back(i);
// 			CLK(pro, "hard");
// 			bounder.generate(hards);
// 			for (auto h : hards){
// 				bounder.set(h);
// 				STP(pro, "hard");
// 				deb(spq->executable(), spq);
// 				CLK(pro, "taylorinit");
// 				Taylor taylor (spq, {}, {
// 					{"soft_deterministic_constraint", false},
// 					{"max_number_of_iterations", 50},
// 					{"dependency_var", true}});
// 				STP(pro, "taylorinit");
// 				CLK(pro, "taylor");
// 				taylor.solve();
// 				STP(pro, "taylor");
// 				deb(filePath, h, taylor.status);
// 				SPQChecker chk (spq);
// 				chk.display(taylor.getSol());
// 			}
// 			// chk.display({{104,36.7424},{984,38.2576}});
// 		}
// 		STOP("a");
// 		// PRINT(pro);
// 	}
// }

#include "oneapi/tbb/concurrent_unordered_map.h"
using oneapi::tbb::concurrent_unordered_map;
#include <pcg_random.hpp>
#include "util/uconfig.hpp"

void testTBB()
{
	// concurrent_unordered_map<int, int> cmap;
	// cmap[0] = 1;
	// deb(cmap);
	// pcg32 gen (Config::getInstance()->seed());
	// vector<unsigned int> seeds;
	// for (int i = 0; i < 80; i ++) seeds.push_back(gen());
	// deb(seeds);
}

void testOmp()
{
	size_t N = 10000000;
	size_t R = 10;
	std::vector<double> v(N);
	std::iota(v.begin(), v.end(), 0);
	INIT(pro);
	for (size_t r = 0; r < R; ++r)
	{
		CLK(pro, "a");
		double norm = 0;
		std::vector<double> vv(N);
		for (size_t i = 0; i < N; ++i)
		{
			norm += v[i] * v[i];
		}
		norm = sqrt(norm);
		for (size_t i = 0; i < N; ++i)
		{
			vv[i] = v[i] / norm;
		}
		STP(pro, "a");
	}
	int core = 80;
	for (size_t r = 0; r < R; ++r)
	{
		CLK(pro, "b");
		std::vector<double> vv(N);
		double norm;
#pragma omp parallel num_threads(core)
		{
			double norm_ = 0;
#pragma omp for nowait
			for (size_t i = 0; i < N; ++i)
			{
				norm_ += v[i] * v[i];
			}
#pragma omp atomic
			norm += norm_;
#pragma omp barrier
#pragma omp master
			{
				norm = sqrt(norm);
			}
#pragma omp barrier
#pragma omp for nowait
			for (size_t i = 0; i < N; ++i)
			{
				vv[i] = v[i] / norm;
			}
		}
		STP(pro, "b");
	}
	PRINT(pro);
}

void testNumeric()
{
	std::vector<double> v = {0, 0};
	deb(normalize(v));
	deb(v);
}

void testgb()
{
	GRBenv *env = NULL;
	GRBmodel *model = NULL;
	ckg(GRBemptyenv(&env), env);
	ckg(GRBsetintparam(env, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_OFF), env);
	ckg(GRBsetintparam(env, GRB_INT_PAR_SIFTING, 2), env);
	ckg(GRBsetintparam(env, GRB_INT_PAR_METHOD, GRB_METHOD_DUAL), env);
	// ckgb(GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0), env);
	ckg(GRBstartenv(env), env);
	int numvars = 2;
	std::vector<double> obj = {1, 1};
	std::vector<double> lb = {0, 1};
	std::vector<double> ub = {4, GRB_INFINITY};
	std::vector<char> vtype = {GRB_CONTINUOUS, GRB_CONTINUOUS};
	ckgb(GRBnewmodel(env, &model, NULL, numvars, obj.data(), lb.data(), ub.data(), vtype.data(), NULL), env, model);
	ckgb(GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE), env, model);
	int numconstrs = 4;
	int numnz = 7;
	std::vector<int> cbeg = {0, 1, 3, 5, 7};
	std::vector<int> cind = {1, 0, 1, 0, 1, 0, 1};
	std::vector<double> cval = {1, 1, 2, 1, 2, 3, 2};
	std::vector<char> sense = {GRB_LESS_EQUAL, GRB_GREATER_EQUAL, GRB_LESS_EQUAL, GRB_GREATER_EQUAL};
	std::vector<double> rhs = {7, 5, 15, 6};
	ckgb(GRBaddconstrs(model, numconstrs, numnz, cbeg.data(), cind.data(), cval.data(), sense.data(), rhs.data(), NULL), env, model);
	ckgb(GRBupdatemodel(model), env, model);
	ckgb(GRBoptimize(model), env, model);
	int status;
	ckgb(GRBgetintattr(model, GRB_INT_ATTR_STATUS, &status), env, model);
	if (status == GRB_OPTIMAL)
	{
		std::vector<double> sol(numvars);
		ckgb(GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, numvars, sol.data()), env, model);
		deb(sol);
	}
	if (model)
		GRBfreemodel(model);
	if (env)
		GRBfreeenv(env);
}

int main()
{
	// testgb();
	// testNumeric();
	// analyzeAll();
	// test();

	//testSSvsNaive();

	testStochDualRed();
	//experimentEpsilonvsObj(2);
	//testNaive();

	//debugNaive();

	// testOmp();
	// testTBB();
	// testHighs();
	// map<string, Option> ok = {{"ok","1"},{"ok1",1.5},{"ok2",false}};
	// deb(ok);
	// for (auto p : ok) deb(p.second.which());
}