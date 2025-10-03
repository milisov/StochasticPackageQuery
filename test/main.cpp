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
#include "nlohmann/json.hpp"
#include "spq/formulator.hpp"
#include "spq/rsformulator.hpp"
#include "solver/RobustSatisficing.hpp"
#include "util/data.hpp"

using std::map;
using std::vector;
using json = nlohmann::json;

#include "solver/taylor.hpp"

// void testNaive(string path, int M_input, string outPath)
// {
// 	int M = M_input;
// 	int M_hat = 1e6;
// 	double epsilon = 0.46;

// 	string filePath = path;
// 	auto spq = parseSpaqlFromFile(filePath);


// 	if (spq)
// 	{
// 		Data::init(spq);  
//     	Data::getInstance().fetchData();
// 		spq->validate();
// 		unique_ptr<Stat> stat = std::make_unique<Stat>();
// 		stat->analyze(spq->tableName);
// 		size_t N = 10000;
// 		double E = 50;
// 		// N is number of scenarios in order to approx, E - expected package size in sol
// 		// set values that are variables in the query
// 		Bounder bounder(spq, N, E);
// 		deb(spq->executable(), spq);
// 		std::vector<std::string> headers = {"Hardness", "Naive-objective", "Naive-Feas", "RuntimeNaive"};
// 		string output = "/home/fm2288/StochasticPackageQuery/test/Experiments2/Naive/Naive" + outPath + ".csv";
// 		DataWriter writer(output, headers);

// 		Profiler stopwatchNaive;

// 		for (int h = -4; h <= -2; h = h + 1)
// 		{
// 			double NaiveObjective;
// 			bool NaiveFeas;
// 			bounder.set(h);
// 			deb(spq);
// 			Naive Nai(M, M_hat, spq);

// 			string labelNaive = "Naive with Hardness" + std::to_string(h);
// 			stopwatchNaive.clock(labelNaive);

// 			FormulateOptions formOptions;
// 			formOptions.reduced = false;

// 			DecisionVarOptions decVarOptions;

// 			decVarOptions.lb = 0.0;
// 			decVarOptions.ub = 1.0;
// 			decVarOptions.obj = 0.0;
// 			decVarOptions.varType = GrbVarType::Binary;
// 			formOptions.decisionVarOptions = decVarOptions;
// 			Solution solu = Nai.solveNaive(Nai.spq, 10, M, M_hat, formOptions);
// 			stopwatchNaive.stop(labelNaive);
// 			double totalTimeNaive = stopwatchNaive.getTime(labelNaive);

// 			SPQChecker Check(Nai.spq);
// 			SolType resNaive;
// 			double distance = 0;
// 			if (solu.x.size() > 0)
// 			{
// 				for (int i = 0; i < Nai.NTuples; i++)
// 				{
// 					resNaive[i + 1] = solu.x[i];
// 				}
// 				NaiveObjective = Check.getObjective(resNaive);
// 				NaiveFeas = Check.feasible(resNaive, distance);
// 			}
// 			else
// 			{
// 				NaiveObjective = -1;
// 				NaiveFeas = 0;
// 			}
// 			//{"Hardness","Naive-objective","Naive-Feas","TimeNaive"}
// 			writer.addRow(h, NaiveObjective, distance, totalTimeNaive);
// 		}
// 	}
// }

// void testSummarySearch(string path, int M_input, string outPath)
// {
// 	int M = M_input;
// 	int M_hat = 1e6;
// 	double epsilon = 0.46;

// 	map<string, Option> countConstraintOptions = {
// 		{"omit count constraint", false},
// 		{"scale down", false},
// 		{"scale factor", -1.0},
// 	};

// 	map<string, Option> curveFitOptions = {
// 		{"arctan", false},
// 		{"binarySearch", true}};

// 	string filePath = path;
// 	auto spq = parseSpaqlFromFile(filePath);
	
// 	if (spq)
// 	{
// 		Data::init(spq);  
//     	Data::getInstance().fetchData();
// 		spq->validate();
// 		unique_ptr<Stat> stat = std::make_unique<Stat>();
// 		stat->analyze(spq->tableName);
// 		size_t N = 10000;
// 		double E = 50;
// 		// N is number of scenarios in order to approx, E - expected package size in sol
// 		// set values that are variables in the query
// 		Bounder bounder(spq, N, E);
// 		deb(spq->executable(), spq);
// 		std::vector<std::string> headers = {"Hardness", "SS-objective", "SS-Feas", "Z", "RuntimeSummarySearch"};
// 		string output = "/home/fm2288/StochasticPackageQuery/test/Experiments2/SummarySearch/SS" + outPath + ".csv";
// 		DataWriter writer(output, headers);

// 		Profiler stopwatchSS;

// 		vector<int> reducedIds;
// 		bool reduced = false;
// 		for (int h = 4; h <= 4; h = h + 1)
// 		{
// 			double SSObjective;
// 			bool SSFeas;
// 			int Z;
// 			bounder.set(h);
// 			deb(spq);
// 			SummarySearch SS(M, spq, epsilon);

// 			string labelSS = "SummarySearch with Hardness" + std::to_string(h);
// 			stopwatchSS.clock(labelSS);
// 			SSFormulator formulator(SS.spq);
// 			FormulateOptions formOptions;
// 			formOptions.reduced = false;
// 			formOptions.Z = 1;
// 			DecisionVarOptions decVarOptions;
// 			setDecisionVarOptions(decVarOptions, 0.0, 1.0, 0.0, GrbVarType::Integer);
// 			formOptions.decisionVarOptions = decVarOptions;

// 			SolutionMetadata<int> sol = SS.summarySearch<int>(SS.spq, formulator, formOptions, curveFitOptions, 1);
// 			stopwatchSS.stop(labelSS);
// 			double totalTimeSS = stopwatchSS.getTime(labelSS);

// 			SPQChecker Check(SS.spq);
// 			double distance;
// 			if (sol.x.size() > 0)
// 			{
// 				deb(sol.x);
// 				SolType res;
// 				for (int i = 0; i < SS.NTuples; i++)
// 				{
// 					res[i + 1] = sol.x[i];
// 				}
// 				SSFeas = Check.feasible(res, distance);
// 				SSObjective = Check.getObjective(res); // we need to change this to validation objective
// 				Z = sol.Z;
// 			}
// 			else
// 			{
// 				SSObjective = -1;
// 				SSFeas = 0;
// 				Z = sol.Z;
// 			}
// 			writer.addRow(h, SSObjective, distance, Z, totalTimeSS);
// 		}
// 	}
// }

void testRS(string path, int M_input, string outPath)
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
		Data::init(spq);  
		Data::getInstance().fetchData();
		spq->validate();
		unique_ptr<Stat> stat = std::make_unique<Stat>();
		stat->analyze(spq->tableName);
		size_t N = 10000;
		double E = 50;

		Bounder bounder(spq, N, E);
		//deb(spq->executable(), spq);
		std::vector<std::string> headers = {"Hardness", "RS-objective", "RS-Feas", "Z", "q","nudge","bestEps", "finalEps","RuntimeRS","TimeStage1", "TimeStage2"};
		string output = "/home/fm2288/StochasticPackageQuery/test/Experiments2/RobustSatisficing/RS" + outPath + ".csv";
		DataWriter writer(output, headers);

		Profiler stopwatchRS;

		vector<int> reducedIds;
		bool reduced = false;

		double RSObjective;
		bool RSFeas;
		int Z; 
		int q;


		for (int h = 2; h <= 4; h = h + 1)
		{
			double RSObjective;
			bool RSFeas;
			int Z;
			bounder.set(h);
			RobustSatisficing RS(M, spq, epsilon);

			string labelRS = "Robust Satisficing with Hardness" + std::to_string(h);
			stopwatchRS.clock(labelRS);
			for(double nudge = 1; nudge <= 1; nudge *= 10)
			{
				SolutionMetadata<int> sol = RS.stochasticDualReducer(RS.spq, curveFitOptions, nudge);
				stopwatchRS.stop(labelRS);
				double totalTimeRS = stopwatchRS.getTime(labelRS);
				double timeStage1 = sol.timeStage1;
				double timeStage2 = sol.timeStage2;
				SPQChecker Check(RS.spq);
				double distance;
				if (sol.x.size() > 0)
				{
					deb(sol.x);
					SolType res;
					for (int i = 0; i < RS.NTuples; i++)
					{
						res[i + 1] = sol.x[i];
					}
					RSFeas = Check.feasible(res, distance);
					RSObjective = Check.getObjective(res);
					Z = sol.Z;
					q = sol.qSz;
				}
				else
				{
					RSObjective = -1;
					RSFeas = 0;
					Z = sol.Z;
					q = sol.qSz;
				}
				writer.addRow(h, RSObjective, distance, Z, q, nudge,sol.bestEps, sol.objConsValue, totalTimeRS, timeStage1, timeStage2);
			}
		}
	}
}

// void testStochDualRed(string path, int M_input, string outPath)
// {
// 	map<string, Option> countConstraintOptions = {
// 		{"omit count constraint", true},
// 		{"scale down", false},
// 		{"scale factor", -1.0},
// 	};

// 	map<string, Option> curveFitOptions = {
// 		{"arctan", false},
// 		{"binarySearch", true}};

// 	int M = M_input;
// 	int M_hat = 1e6;
// 	double epsilon = 0.46;
// 	int qSz = 500;

// 	string filePath = path;
// 	// convert query to class
// 	auto spq = parseSpaqlFromFile(filePath);
// 	if (spq)
// 	{
// 		Data::init(spq);  
// 		Data::getInstance().fetchData();
// 		spq->validate();
// 		unique_ptr<Stat> stat = std::make_unique<Stat>();
// 		stat->analyze(spq->tableName);
// 		size_t N = 10000;
// 		double E = 50;

// 		Bounder bounder(spq, N, E);
// 		deb(spq->executable(), spq);
// 		std::vector<std::string> headers = {"hardness", "SDR-Objective", "SDR-feas", "SDR-optimal", "Z", "qSz", "Binary Search Steps", "Runtime"};
// 		string output = "/home/fm2288/StochasticPackageQuery/test/Experiments2/SDR/SDR" + outPath + ".csv";
// 		cout << output << endl;
// 		DataWriter writer(output, headers);

// 		Profiler stopwatchSDR;

// 		vector<int> reducedIds;
// 		bool reduced = false;
// 		for (int h = 4; h <= 4; h = h + 1)
// 		{
// 			bounder.set(h);
// 			deb(spq);
// 			SummarySearch SS(M, spq, epsilon);

// 			FormulateOptions formOptions;
// 			formOptions.reduced = false;
// 			formOptions.Z = 1;
// 			DecisionVarOptions decVarOptions;
// 			setDecisionVarOptions(decVarOptions, 0.0, 1.0, 0.0, GrbVarType::Binary);
// 			formOptions.decisionVarOptions = decVarOptions;
// 			formOptions.qSz = qSz;

// 			string labelSS = "SummarySearch with Hardness" + std::to_string(h);
// 			stopwatchSDR.clock(labelSS);
// 			SolutionMetadata<int> sol = SS.stochDualReducer(SS.spq, formOptions, curveFitOptions);
// 			stopwatchSDR.stop(labelSS);
// 			double totalTimeSS = stopwatchSDR.getTime(labelSS);

// 			if (sol.x.size() > 0)
// 			{
// 				deb(sol.x);

// 				SPQChecker Check(SS.spq);
// 				SolType res;
// 				for (int i = 0; i < SS.NTuples; i++)
// 				{
// 					res[i + 1] = sol.x[i];
// 				}

// 				double distance;
// 				bool feas = Check.feasible(res, distance);
// 				double obj = Check.getObjective(res);

// 				cout << (feas ? "Validation Feasible Solution" : "Validation Infeasible Solution") << endl;
// 				cout << "Training Objective Summary Search = " << sol.w << endl;
// 				cout << "Validation Objective Summary Search = " << obj << endl;

// 				writer.addRow(h, obj, distance, sol.isOptimal, sol.Z, sol.qSz, sol.binarySearchSteps, totalTimeSS);
// 			}
// 		}
// 	}
// }

void generateQuerywithHardness(string path, string outPath, int n, int m)
{
	string filePath = path;
	auto spq = parseSpaqlFromFile(filePath);

	if (spq)
	{
		string output = "/home/fm2288/StochasticPackageQuery/test/Queries/" + outPath + "/" + outPath;
		spq->validate();
		unique_ptr<Stat> stat = std::make_unique<Stat>();
		stat->analyze(spq->tableName);
		size_t N = 10000;
		double E = 50;
		Bounder bounder(spq, N, E);

		for (int h = -4; h <= 4; h = h + 1)
		{
			bounder.set(h);
			std::string filename = output + "_" + std::to_string(h) + ".spaql";
			std::ofstream outFile(filename);

			if (!outFile)
			{
				std::cerr << "Error: Could not open " << filename << " for writing.\n";
				continue; // Skip this iteration if file couldn't be opened
			}
			std::string queryString = static_cast<std::string>(*spq);
			outFile << queryString << std::endl;
			outFile.close();
			std::cout << "Saved query for h = " << h << " to " << filename << std::endl;
		}
	}
}

void validation(string filepath, string outpath)
{
	std::vector<std::string> headers = {"hardness", "objective", "feas", "Runtime"};

	// Open JSON file
	std::ifstream file(filepath);
	if (!file)
	{
		std::cerr << "Error: Could not open file!" << std::endl;
		return;
	}

	json jsonData;
	file >> jsonData;

	// Check if jsonData is an array (multiple records) or a single object
	if (!jsonData.is_array())
	{
		jsonData = json::array({jsonData}); // Wrap single object in an array for uniform processing
	}

	string output = "/home/fm2288/StochasticPackageQuery/test/Experiments2/RCL/RCL" + outpath + ".csv";
	cout << output << endl;
	DataWriter writer(output, headers);

	for (const auto &record : jsonData)
	{
		std::cout << "Processing new record..." << std::endl;

		string query = record["Query"];
		string dbPath = fmt::format("resource/sqls/_{}.spaql", query);
		int h = record["Hardness"];

		auto spq = parseSpaqlFromFile(dbPath);
		if (!spq)
		{
			std::cerr << "Error: Could not parse SPARQL file!" << std::endl;
			continue;
		}

		spq->validate();
		unique_ptr<Stat> stat = std::make_unique<Stat>();
		stat->analyze(spq->tableName);
		size_t N = 10000;
		double E = 50;
		Bounder bounder(spq, N, E);
		bounder.set(h);
		deb(spq);

		SPQChecker Check(spq);
		SolType res;
		PgManager pg;
		int nTuples = pg.getTableSize(spq->tableName);
		std::cout << "Number of Tuples: " << nTuples << std::endl;

		for (int i = 0; i < nTuples; i++)
		{
			res[i + 1] = 0;
		}

		for (const auto &[key, value] : record["Package"].items())
		{
			int intKey = std::stoi(key);
			res[intKey] = static_cast<int>(value);
		}
		double distance;
		bool feas = Check.feasible(res, distance);
		double objective = Check.getObjective(res);
		std::cout << (feas ? "Validation Feasible Solution" : "Validation Infeasible Solution") << std::endl;
		std::cout << "Validation Objective = " << objective << std::endl;
		writer.addRow(h, objective, distance, record["Runtime"]);
	}
}

int main(int argc, char *argv[])
{

	if (argc < 4)
	{
		std::cerr << "Usage: " << argv[0] << " <N> <M> <algorithm>" << std::endl;
		return 1;
	}

	// validation(path);
	// return 1;
	int N = std::stoi(argv[1]);
	int M = std::stoi(argv[2]);
	std::string algorithm = argv[3];

	std::string dbPath = fmt::format("resource/sqls/_stocks_{}_{}.spaql", N, M);
	std::string outPath = fmt::format("stocks_{}_{}", N, M);

	// std::cout << "File Path: " << dbPath << std::endl;
	// std::cout << "Algorithm: " << algorithm << std::endl;

	if (algorithm == "Naive")
	{
		//testNaive(dbPath, M, outPath);
	}
	else if (algorithm == "SummarySearch")
	{
		//testSummarySearch(dbPath, M, outPath);
	}
	else if (algorithm == "SDR")
	{
		//testStochDualRed(dbPath, M, outPath);
	}
	else if (algorithm == "RS")
	{
		testRS(dbPath, M, outPath);
	}
	else if (algorithm == "generate")
	{
		generateQuerywithHardness(dbPath, outPath, N, M);
	}
	else if (algorithm == "validate")
	{
		string path = fmt::format("/home/fm2288/StochasticPackageQuery/stochastic-sketchrefine/{}.json", outPath);
		validation(path, outPath);
	}
	else
	{
		cout << "Please enter a correct algorithm name next time" << endl;
	}

	return 0;
}