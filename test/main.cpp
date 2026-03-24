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
#include "solver/StochDualReducer.hpp"
#include "solver/Naive.hpp"
#include "gurobi_c++.h"
#include <gurobi_c.h>
#include "nlohmann/json.hpp"
#include "spq/formulator.hpp"
#include "spq/rsformulator.hpp"
#include "solver/RobustSatisficing.hpp"
#include "util/data.hpp"
#include <sstream>
#include <stdexcept>
#include <dirent.h>
#include <sys/stat.h> // For mkdir
#include <cerrno>
#include <iostream>
#include <string>
#include <fstream>
#include <memory>
// Required for creating directories on POSIX systems (like Linux/macOS)
#include <sys/stat.h>
#include <sys/types.h>
#include <chrono>
#include <sys/file.h> // For flock
#include <fcntl.h>	  // For open, O_CREAT, etc.
#include <unistd.h>	  // For close, write, fork
#include <sys/wait.h>
#include <thread>
#include <fstream>
#include "nlohmann/json.hpp"
using std::map;
using std::vector;
using json = nlohmann::json;

void testRS(string path, int M_input, string outPath)
{
	throw("uncomment");
	// int M = M_input;
	// int M_hat = 1e6;
	// double epsilon = 0.46;

	// map<string, Option> countConstraintOptions = {
	// 	{"omit count constraint", false},
	// 	{"scale down", false},
	// 	{"scale factor", -1.0},
	// };

	// map<string, Option> curveFitOptions = {
	// 	{"arctan", false},
	// 	{"binarySearch", true}};

	// string filePath = path;
	// auto spq = parseSpaqlFromFile(filePath);

	// if (spq)
	// {
	// 	Data::init(spq);
	// 	Data::getInstance().fetchData();
	// 	spq->validate();
	// 	unique_ptr<Stat> stat = std::make_unique<Stat>();
	// 	deb("analyzing:", spq->tableName);
	// 	string tableForValidation = fmt::format("{}_validate", spq->tableName);
	// 	stat->analyze(spq->tableName);
	// 	stat->analyze(tableForValidation);
	// 	size_t N = 10000;
	// 	double E = 50;

	// 	Bounder bounder(spq, N, E);
	// 	// deb(spq->executable(), spq);
	// 	std::vector<std::string> headers = {"Hardness", "Objective", "feas", "surplus", "bestEps", "RuntimeRS", "solutions"};
	// 	string output = "/home/fm2288/StochasticPackageQuery/test/Experiments2/RobustSatisficing/RS" + outPath + ".csv";
	// 	DataWriter writer(output, headers);

	// 	Profiler stopwatchRS;

	// 	vector<int> reducedIds;
	// 	bool reduced = false;

	// 	double RSObjective;
	// 	bool RSFeas;
	// 	int Z;
	// 	int q;
	// 	for (int h = 0; h <= 8; h = h + 1)
	// 	{
	// 		double RSObjective;
	// 		bool RSFeas;
	// 		int Z;
	// 		bounder.set(h);
	// 		deb(spq);
	// 		RobustSatisficing RS(M, spq, epsilon);

	// 		string labelRS = "Robust Satisficing with Hardness" + std::to_string(h);
	// 		stopwatchRS.clock(labelRS);

	// 		auto start_time = std::chrono::steady_clock::now();
	// 		double timeout_seconds = 600;
	// 		SolutionMetadata<int> sol = RS.stochasticDualReducer(RS.spq, curveFitOptions, timeout_seconds);
	// 		stopwatchRS.stop(labelRS);
	// 		double totalTimeRS = stopwatchRS.getTime(labelRS);
	// 		double totalValidationTime = gpro.getTime("optimValidation");
	// 		gpro.reset("optimValidation");
	// 		totalTimeRS -= totalValidationTime;
	// 		double timeStage1 = sol.timeStage1;
	// 		double timeStage2 = sol.timeStage2;
	// 		SPQChecker Check(RS.spq);
	// 		vector<double> feasibility;
	// 		vector<double> surpluses;
	// 		if (sol.x.size() > 0)
	// 		{
	// 			deb(sol.x);
	// 			SolType res;
	// 			for (int i = 0; i < RS.NTuples; i++)
	// 			{
	// 				res[i + 1] = sol.x[i];
	// 			}
	// 			RSFeas = Check.feasible(res, feasibility, surpluses);
	// 			RSObjective = Check.getObjective(res);
	// 			Z = sol.Z;
	// 			q = sol.qSz;
	// 		}
	// 		else
	// 		{
	// 			RSObjective = -1;
	// 			RSFeas = 0;
	// 			Z = sol.Z;
	// 			q = sol.qSz;
	// 		}
	// 		deb("writing in a file");
	// 		writer.addRow(h, RSObjective, feasibility, surpluses, sol.bestEps, totalTimeRS, sol.solutions);
	// 		deb("done writing");
	// 	}
	// }
}

void testRSSeed(int M_input, int N_input)
{
	int M = M_input;
	int N = N_input;
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

	std::string outPath = fmt::format("stocks_{}_{}", N, M);
	std::vector<std::string> headers = {"Hardness", "Seed", "Objective", "feas", "surplus", "tuples", "scenarios", "RuntimeRS","solutions"};
	string output = "/home/fm2288/StochasticPackageQuery/test/ExperimentsSeedTimeBudgetNoObj/RobustSatisficing/RS" + outPath + ".csv";
	DataWriter writer(output, headers);
	writer.writeHeaders();

	int max_concurrent = 15;
	int active_processes = 0;

	for (int h = 8; h <= 8; h++)
	{
		for (int seed = 1; seed <= 10; seed = seed + 1)
		{
			// check if any existing children finished
			int status;
			while (waitpid(-1, &status, WNOHANG) > 0)
			{
				active_processes--;
			}

			// If we are full, we MUST pause here until at least one child finishes.
			if (active_processes >= max_concurrent)
			{
				wait(NULL);
				active_processes--;
			}

			pid_t pid = fork();
			if (pid < 0)
			{
				// fork failed
				exit(1);
			}
			else if (pid == 0)
			{
				// child process;
				Profiler stopwatchRS;
				double RSObjective;
				bool RSFeas;
				int Z;
				int q;
				int qScenarios;
				std::string queryPath = fmt::format("resource/sqls/_stocks_{}_{}_seeded_{}.spaql", N, M, seed);
				string tableForHardness = fmt::format("stocks_{}_{}_validate", N, M);

				auto spq = parseSpaqlFromFile(queryPath);
				deb(spq);
				if (spq)
				{
					double fetchingTime = 0.0;
					// we need to consider the fetching into effective runtime
					{
						gpro.clock("effectiveRuntime");
						Data::reset();
						Data::init(spq);
						Data::getInstance().fetchData();
						gpro.stop("effectiveRuntime");
						fetchingTime = gpro.getTime("effectiveRuntime");

						deb(fetchingTime);


						spq->validate();
						unique_ptr<Stat> stat = std::make_unique<Stat>();

						size_t N = 10000;
						double E = 50;
						string originalTableName = spq->tableName;
						spq->setTableName(tableForHardness);
						Bounder bounder(spq, N, E);
						bounder.set(h);
						deb(spq);
						spq->setTableName(originalTableName);
						string tableForValidation = fmt::format("{}_validate", spq->tableName);
						deb(spq);
					}
					gpro.clock("effectiveRuntime");
					RobustSatisficing RS(M, spq, epsilon);
					double timeout_seconds = 10 * 60;
					

					//this is for the time budget experiment
					SolveOptions solveOptions;
					solveOptions.timeout_seconds = timeout_seconds;
					solveOptions.timeBudgeted= true;
					solveOptions.timeBudget = timeout_seconds;

					SolutionMetadata<int> sol = RS.stochasticDualReducer(RS.spq, solveOptions);
					gpro.stop("effectiveRuntime");
					double totalTimeRS = gpro.getTime("effectiveRuntime");

					{
						SPQChecker Check(RS.spq);
						vector<double> feasibility;
						vector<double> surpluses;
						if (sol.x.size() > 0)
						{
							deb(sol.x);
							SolType res;
							for (int i = 0; i < RS.NTuples; i++)
							{
								if (sol.x[i] > 0)
								{
									res[i + 1] = sol.x[i];
								}
							}
							RSFeas = Check.feasible(res, feasibility, surpluses);
							RSObjective = Check.getObjective(res);
							deb("printing-gpro");
							gpro.print();
							Z = sol.Z;
							q = sol.qSz;
							qScenarios = sol.qScenarios;
						}
						else
						{
							RSObjective = -1;
							RSFeas = 0;
							Z = sol.Z;
							q = sol.qSz;
						}
						int lock_fd = open(output.c_str(), O_WRONLY | O_CREAT | O_APPEND, 0666);
						if (lock_fd != -1)
						{
							flock(lock_fd, LOCK_EX); // LOCK
							double gurobiTime = gpro.getTime("gurobi");
							//{"Hardness", "Seed", "Objective", "feas", "surplus", "tuples", "scenarios", "RuntimeRS","solutions"};
							writer.addRow(h, seed, RSObjective, feasibility, surpluses, q, qScenarios, totalTimeRS, sol.solutions);

							flock(lock_fd, LOCK_UN); // UNLOCK
							close(lock_fd);
						}
						else
						{
							std::cerr << "Could not open file for locking!" << std::endl;
						}
						gpro.reset("effectiveRuntime");
					}
				}
				deb("Child finished successfully. Exiting.");
				exit(0);
			}
			else
			{
				active_processes++;
			}
		}
	}
	// the loop is done but 40 children are still running. Wait for them
	while (active_processes > 0)
	{
		wait(NULL);
		active_processes--;
	}
}

// void testRSSeed(int M_input, int N_input)
// {
// 	deb("here");
// 	int M = M_input;
// 	int N = N_input;
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

// 	std::string outPath = fmt::format("stocks_{}_{}", N, M);
// 	std::vector<std::string> headers = {"Hardness", "Seed", "Objective", "feas", "surplus", "tuples", "bestEps", "RuntimeRS", "solutions"};
// 	string output = "/home/fm2288/StochasticPackageQuery/test/ExperimentsSeed/RobustSatisficing/RS" + outPath + ".csv";
// 	DataWriter writer(output, headers);
// 	for (int h = 0; h <= 8; h++)
// 	{
// 		for (int seed = 1; seed <= 10; seed = seed + 1)
// 		{
// 			Profiler stopwatchRS;
// 			double RSObjective;
// 			bool RSFeas;
// 			int Z;
// 			int q;
// 			std::string queryPath = fmt::format("resource/sqls/_stocks_{}_{}_seeded_{}.spaql", N, M, seed);
// 			string tableForHardness = fmt::format("stocks_{}_{}_validate", N, M);

// 			auto spq = parseSpaqlFromFile(queryPath);
// 			deb(spq);
// 			if (spq)
// 			{
// 				//we need to consider the fetching into effective runtime
// 				gpro.clock("effectiveRuntime");
// 				Data::reset();
// 				Data::init(spq);
// 				Data::getInstance().fetchData();
// 				gpro.stop("effectiveRuntime");

// 				deb("validating");
// 				spq->validate();
// 				deb("done");
// 				// unique_ptr<Stat> stat = std::make_unique<Stat>();

// 				size_t N = 10000;
// 				double E = 50;
// 				string originalTableName = spq->tableName;
// 				deb(tableForHardness);
// 				// stat->analyze(tableForHardness);
// 				spq->setTableName(tableForHardness);
// 				Bounder bounder(spq, N, E);
// 				deb("bounding");
// 				bounder.set(h);
// 				deb(spq);
// 				spq->setTableName(originalTableName);
// 				string tableForValidation = fmt::format("{}_validate", spq->tableName);
// 				// stat->analyze(spq->tableName);
// 				deb(spq->tableName);
// 				// stat->analyze(tableForValidation);
// 				deb(tableForValidation);
// 				deb(spq);

// 				gpro.clock("effectiveRuntime");
// 				RobustSatisficing RS(M, spq, epsilon);
// 				double timeout_seconds = 30*60;
// 				SolutionMetadata<int> sol = RS.stochasticDualReducer(RS.spq, curveFitOptions, timeout_seconds);
// 				gpro.stop("effectiveRuntime");
// 				double totalTimeRS = gpro.getTime("effectiveRuntime");

// 				SPQChecker Check(RS.spq);
// 				vector<double> feasibility;
// 				vector<double> surpluses;
// 				if (sol.x.size() > 0)
// 				{
// 					deb(sol.x);
// 					SolType res;
// 					for (int i = 0; i < RS.NTuples; i++)
// 					{
// 						if(sol.x[i] > 0)
// 						{
// 							res[i + 1] = sol.x[i];
// 						}
// 					}
// 					RSFeas = Check.feasible(res, feasibility, surpluses);
// 					RSObjective = Check.getObjective(res);
// 					deb("printing-gpro");
// 					gpro.print();
// 					Z = sol.Z;
// 					q = sol.qSz;
// 				}
// 				else
// 				{
// 					RSObjective = -1;
// 					RSFeas = 0;
// 					Z = sol.Z;
// 					q = sol.qSz;
// 				}
// 				deb("writing in a file");
// 				writer.addRow(h, seed, RSObjective, feasibility, surpluses, q, sol.bestEps, totalTimeRS, sol.solutions);
// 				gpro.reset("effectiveRuntime");
// 				deb("done writing");
// 			}
// 		}
// 	}
// }

#include <fstream>
#include "nlohmann/json.hpp" // Make sure to have this header available

// For convenience
using json = nlohmann::json;

void testRSExp2(const std::string &basePath, int M_input, int hardness, const std::string &outPath)
{
	// for(int h = hardness; h <= 4; h++)
	// {
	// 	// --- 1. Construct the target directory path to read queries from ---
	// 	std::string targetDirectory = basePath + "/RS/" + std::to_string(h);
	// 	std::cout << "Reading queries from directory: " << targetDirectory << std::endl;

	// 	// --- 2. Set up the single CSV writer for all results ---
	// 	std::vector<std::string> headers = {"Hardness", "p", "RS-objective", "RS-Feas", "bestEps", "RuntimeRS"};
	// 	std::string outputCsvFile = "/home/fm2288/StochasticPackageQuery/test/ExperimentsAUC/RobustSatisficing/RS" + outPath + ".csv";
	// 	DataWriter writer(outputCsvFile, headers);
	// 	std::cout << "Writing results to: " << outputCsvFile << std::endl;

	// 	// --- 3. Open and iterate through the target directory ---
	// 	DIR *dir = opendir(targetDirectory.c_str());
	// 	if (dir == nullptr)
	// 	{
	// 		std::cerr << "Error: Could not open directory " << targetDirectory << std::endl;
	// 		return;
	// 	}

	// 	struct dirent *entry;
	// 	while ((entry = readdir(dir)) != nullptr)
	// 	{
	// 		std::string filename = entry->d_name;

	// 		// Process only files that end with ".spaql"
	// 		if (filename.length() <= 6 || filename.substr(filename.length() - 6) != ".spaql")
	// 		{
	// 			continue;
	// 		}

	// 		std::string fullQueryPath = targetDirectory + "/" + filename;
	// 		std::cout << "\nProcessing file: " << fullQueryPath << std::endl;
	// 		double p;
	// 		try
	// 		{
	// 			size_t last_underscore = filename.rfind('_');
	// 			size_t dot_pos = filename.rfind(".spaql");
	// 			if (last_underscore == std::string::npos || dot_pos == std::string::npos)
	// 			{
	// 				throw std::invalid_argument("Filename format incorrect");
	// 			}
	// 			std::string p_str = filename.substr(last_underscore + 1, dot_pos - (last_underscore + 1));
	// 			p = std::stod(p_str) / 100.0;
	// 		}
	// 		catch (const std::exception &e)
	// 		{
	// 			std::cerr << "Warning: Could not parse probability from filename '" << filename << "'. Skipping." << std::endl;
	// 			continue;
	// 		}

	// 		// --- 5. Core Solving Logic (applied to each query file) ---
	// 		auto spq = parseSpaqlFromFile(fullQueryPath);
	// 		deb(spq);
	// 		if (!spq)
	// 		{
	// 			std::cerr << "Warning: Failed to parse query from file '" << fullQueryPath << "'. Skipping." << std::endl;
	// 			continue;
	// 		}

	// 		// Initialize data and validate the parsed query
	// 		Data::init(spq);
	// 		Data::getInstance().fetchData();
	// 		spq->validate();

	// 		int M = M_input;
	// 		double epsilon = 0.46;
	// 		map<string, Option> curveFitOptions = {{"arctan", false}, {"binarySearch", true}};

	// 		RobustSatisficing RS(M, spq, epsilon);
	// 		Profiler stopwatchRS;
	// 		string labelRS = "RobustSatisficing_h" + std::to_string(h) + "_p" + std::to_string((int)(p * 100));

	// 		stopwatchRS.clock(labelRS);
	// 		SolutionMetadata<int> sol = RS.stochasticDualReducer(RS.spq, curveFitOptions);
	// 		stopwatchRS.stop(labelRS);
	// 		double totalTimeRS = stopwatchRS.getTime(labelRS);

	// 		SPQChecker Check(RS.spq);
	// 		double RSObjective = -1.0;
	// 		double feasScore = 0.0;
	// 		bool RSFeas = false;

	// 		if (!sol.x.empty())
	// 		{
	// 			SolType res;
	// 			for (size_t i = 0; i < sol.x.size(); ++i)
	// 			{
	// 				res[i + 1] = sol.x[i];
	// 			}
	// 			deb(res);
	// 			RSFeas = Check.feasible(res, feasScore);
	// 			RSObjective = Check.getObjective(res);

	// 			// --- 6. Save results to a JSON file ---
	// 			json j;
	// 			j["basePath"] = basePath;
	// 			j["hardness"] = h;
	// 			j["p"] = p;

	// 			// Convert the SolType map to a JSON object
	// 			// The nlohmann/json library can't directly handle a map with integer keys,
	// 			// so we convert the keys to strings.
	// 			json res_json;
	// 			for (const auto& pair : res) {
	// 				res_json[std::to_string(pair.first)] = pair.second;
	// 			}
	// 			j["res"] = res_json;

	// 			// Construct a unique filename for the JSON output
	// 			std::string jsonOutputFilename = "RS_solution_h" + std::to_string(h) + "_p" + std::to_string(static_cast<int>(p * 100)) + ".json";
	// 			std::string jsonOutputPath = "/home/fm2288/StochasticPackageQuery/test/ExperimentsAUC/RobustSatisficing/Solutions/" + jsonOutputFilename;

	// 			// Write the JSON object to a file
	// 			std::ofstream o(jsonOutputPath);
	// 			o << std::setw(4) << j << std::endl;
	// 		}
	// 		writer.addRow(h, p, RSObjective, feasScore, sol.bestEps, totalTimeRS);
	// 	}

	// 	closedir(dir);
	// 	std::cout << "\nBatch processing complete." << std::endl;
	// }
}

std::string trim_leading_whitespace(const std::string &str)
{
	// Find the first character that is not a whitespace
	auto first_char = std::find_if(str.begin(), str.end(), [](unsigned char ch)
								   { return !std::isspace(ch); });
	// Return the substring from that character to the end
	return std::string(first_char, str.end());
}

// Helper function to remove the COUNT(*) BETWEEN constraint from the query string
std::string remove_count_constraint(const std::string &original_query)
{
	std::stringstream input_ss(original_query);
	std::string line;
	std::stringstream result_ss;
	bool first_such_that = true;

	while (std::getline(input_ss, line))
	{
		// Trim leading whitespace to standardize the line for checks
		std::string trimmed_line = trim_leading_whitespace(line);

		// Skip the line containing the COUNT(*) BETWEEN constraint
		if (trimmed_line.find("COUNT(*) BETWEEN") == 0)
		{
			continue;
		}

		// Handle the "SUCH THAT" line to avoid it being orphaned if it's the only constraint
		if (trimmed_line.rfind("SUCH THAT", 0) == 0 && first_such_that)
		{
			// If the next line is the count constraint, we might need to adjust.
			// A simple approach is to just add it and let the next check skip the count line.
			result_ss << line << "\n";
			first_such_that = false;
		}
		else
		{
			result_ss << line << "\n";
		}
	}
	return result_ss.str();
}

std::string transform_query(
	const std::string &original_query,
	bool includeCountConstraint,
	bool addRepeatZero)
{
	std::stringstream input_ss(original_query);
	std::string line;
	std::vector<std::string> original_lines;

	// 1. Read input query line by line
	while (std::getline(input_ss, line))
	{
		if (!line.empty())
		{
			original_lines.push_back(line);
		}
	}

	if (original_lines.empty())
	{
		throw std::runtime_error("Input query is empty.");
	}

	std::stringstream result_ss;
	bool such_that_header_added = false;

	// 2. Process each line
	for (const auto &current_line : original_lines)
	{
		std::string trimmed_line = trim_leading_whitespace(current_line);

		// Rule for the "SELECT" line
		if (trimmed_line.rfind("SELECT", 0) == 0)
		{
			size_t package_start = trimmed_line.find("PACKAGE(");
			size_t package_end = trimmed_line.find(")", package_start);
			if (package_start != std::string::npos && package_end != std::string::npos)
			{
				trimmed_line.replace(package_start + 8, package_end - (package_start + 8), "*");
			}

			size_t such_that_pos = trimmed_line.find(" SUCH THAT");
			if (such_that_pos != std::string::npos)
			{
				trimmed_line.replace(such_that_pos, 11, " as P");
			}

			result_ss << trimmed_line << "\n";
			if (addRepeatZero)
			{
				result_ss << "REPEAT 0\n";
			}
			result_ss << "SUCH THAT\n";
			such_that_header_added = true;
		}
		// Rule for the standalone "SUCH THAT" line, which signals the start of constraints
		else if (trimmed_line == "SUCH THAT")
		{
			if (!such_that_header_added)
			{
				result_ss << "SUCH THAT\n";
				such_that_header_added = true;
			}
		}
		// Rule specifically for the "COUNT(*) BETWEEN" line
		else if (trimmed_line.rfind("COUNT(*) BETWEEN", 0) == 0)
		{
			// Only process this line if the flag is true
			if (includeCountConstraint)
			{
				if (!such_that_header_added)
				{
					result_ss << "SUCH THAT\n";
					such_that_header_added = true;
				}
				std::stringstream between_ss(trimmed_line);
				std::string keyword, between_word, lower_bound, and_word, upper_bound;
				between_ss >> keyword >> between_word >> lower_bound >> and_word >> upper_bound;

				result_ss << "COUNT(*) >= " << lower_bound << " AND\n";
				result_ss << "COUNT(*) <= " << upper_bound << " AND\n";
			}
			// If includeCountConstraint is false, we do nothing, effectively skipping/removing the line.
		}
		// Rule for all other lines
		else
		{
			result_ss << trimmed_line << "\n";
		}
	}

	std::string final_query = result_ss.str();
	if (!final_query.empty() && final_query.back() == '\n')
	{
		final_query.pop_back();
	}

	return final_query;
}

bool create_nested_directories(const std::string &path)
{
	mode_t mode = 0755;
	std::string current_path = "";
	std::string dir_name;
	std::stringstream ss(path);

	// Handle absolute paths starting with '/'
	if (!path.empty() && path[0] == '/')
	{
		current_path = "/";
	}

	while (std::getline(ss, dir_name, '/'))
	{
		if (dir_name.empty())
			continue;

		// Append the next directory name, avoiding double slashes
		if (current_path.length() > 1 && current_path.back() != '/')
		{
			current_path += "/";
		}
		current_path += dir_name;

		if (mkdir(current_path.c_str(), mode) != 0)
		{
			if (errno != EEXIST)
			{
				std::cerr << "Error creating directory " << current_path << ": " << strerror(errno) << std::endl;
				return false;
			}
		}
	}
	return true;
}

void generateQueryWithBounds(
	const std::string &path,
	const std::string &outPath,
	int n,
	int m,
	const std::string &experimentName, // New parameter for the experiment directory
	bool includeCountConstraint,	   // Flag to control the count constraint in RS
	bool addRepeatZero				   // Flag to control "REPEAT 0" in RCL
)
{

	std::string baseOutputDir = "/home/fm2288/StochasticPackageQuery/test/" + experimentName + "/" + outPath;
	std::string rsPath = baseOutputDir + "/RS";
	std::string rclPath = baseOutputDir + "/RCL";

	if (!create_nested_directories(rsPath) || !create_nested_directories(rclPath))
	{
		std::cerr << "Failed to create base directories. Aborting." << std::endl;
		return;
	}
	for (int h = 2; h <= 4; ++h)
	{

		auto spq = parseSpaqlFromFile(path);
		spq->validate();
		deb(spq);
		auto stat = std::make_unique<Stat>();
		stat->analyze(spq->tableName);
		size_t N = 10000;
		double E = 50;

		if (!spq)
		{
			std::cerr << "Error: Could not parse file " << path << std::endl;
			return;
		}
		Bounder bounder(spq, N, E);
		bounder.set(h);
		deb(spq);
		// --- 3. Create Hardness Level Subdirectories inside RS and RCL ---
		std::string rsHardnessPath = rsPath + "/" + std::to_string(h);
		std::string rclHardnessPath = rclPath + "/" + std::to_string(h);
		create_nested_directories(rsHardnessPath);
		create_nested_directories(rclHardnessPath);

		for (double p = 0.80; p < 1; p += 0.02)
		{
			for (auto &con : spq->cons)
			{
				std::shared_ptr<ProbConstraint> probCon;
				std::shared_ptr<AttrConstraint> attrCon;
				if (isStochastic(con, probCon, attrCon))
				{
					probCon->p = p;
				}
			}
			deb(spq);
			// --- 4. Construct the correct filename and full file paths ---
			std::string baseFilename = outPath + "_" + std::to_string(static_cast<int>(p * 100)) + ".spaql";
			std::string rsFilePath = rsHardnessPath + "/" + baseFilename;
			std::string rclFilePath = rclHardnessPath + "/" + baseFilename;

			// --- 5. Generate and Save Both Queries ---
			std::string originalQueryString = static_cast<std::string>(*spq);

			// Conditionally remove the count constraint for the RS query
			if (!includeCountConstraint)
			{
				originalQueryString = remove_count_constraint(originalQueryString);
			}

			std::ofstream rsOutFile(rsFilePath);
			if (!rsOutFile)
			{
				std::cerr << "Error: Could not open " << rsFilePath << " for writing.\n";
			}
			else
			{
				rsOutFile << originalQueryString;
				rsOutFile.close();
			}

			// Pass the flag to the transform function
			std::string transformedQueryString = transform_query(static_cast<std::string>(*spq), includeCountConstraint, addRepeatZero);
			std::ofstream rclOutFile(rclFilePath);
			if (!rclOutFile)
			{
				std::cerr << "Error: Could not open " << rclFilePath << " for writing.\n";
			}
			else
			{
				rclOutFile << transformedQueryString;
				rclOutFile.close();
			}

			int percentage = std::round(p * 100);
			std::cout << "Saved queries for h=" << h << ", p=" << percentage << "% to " << rsHardnessPath << " and " << rclHardnessPath << "\n";
		}
	}
}

// Assume the necessary classes and functions like parseSpaqlFromFile, Bounder,
// and transform_query are defined and included elsewhere.

void generateQuerywithHardness(const std::string &path, const std::string &outPath, int n, int m)
{
	std::string filePath = path;
	auto spq = parseSpaqlFromFile(filePath);

	if (spq)
	{
		// 1. Define and create the separate output directories for RS and RCL
		std::string baseOutputDir = "/home/fm2288/StochasticPackageQuery/test/Queries/" + outPath;
		std::string rsDir = baseOutputDir + "/RS/";
		std::string rclDir = baseOutputDir + "/RCL/";

		// Create the directories if they don't exist.
		// mkdir returns 0 on success. The mode 0775 gives read/write/execute permissions to owner/group.
		mkdir(baseOutputDir.c_str(), 0775); // Create the base directory first
		mkdir(rsDir.c_str(), 0775);			// Create the RS subdirectory
		mkdir(rclDir.c_str(), 0775);		// Create the RCL subdirectory

		std::cout << "Output directory for RS (untransformed) queries: " << rsDir << std::endl;
		std::cout << "Output directory for RCL (transformed) queries: " << rclDir << std::endl;

		// 2. Set up the Bounder to modify query constraints based on hardness
		spq->validate();
		std::unique_ptr<Stat> stat = std::make_unique<Stat>();
		stat->analyze(spq->tableName);
		size_t N = 10000;
		double E = 50;
		Bounder bounder(spq, N, E);

		// 3. Loop through each hardness level
		for (int h = -4; h <= -1; ++h)
		{
			// Apply the hardness setting, which modifies the spq object
			bounder.set(h);

			// --- 4. Handle the UNTRANSFORMED query for RS ---
			std::string untransformedQueryString = static_cast<std::string>(*spq);
			std::string rsFilename = rsDir + outPath + "_" + std::to_string(h) + ".spaql";
			std::ofstream rsOutFile(rsFilename);

			if (!rsOutFile)
			{
				std::cerr << "Error: Could not open " << rsFilename << " for writing.\n";
			}
			else
			{
				rsOutFile << untransformedQueryString << std::endl;
				rsOutFile.close();
				std::cout << "Saved RS query for h=" << h << " to " << rsFilename << std::endl;
			}

			// --- 5. Handle the TRANSFORMED query for RCL ---
			std::string transformedQueryString = transform_query(untransformedQueryString, true, true);
			std::string rclFilename = rclDir + outPath + "_" + std::to_string(h) + ".spaql";
			std::ofstream rclOutFile(rclFilename);

			if (!rclOutFile)
			{
				std::cerr << "Error: Could not open " << rclFilename << " for writing.\n";
			}
			else
			{
				rclOutFile << transformedQueryString << std::endl;
				rclOutFile.close();
				std::cout << "Saved RCL query for h=" << h << " to " << rclFilename << std::endl;
			}
		}
	}
	else
	{
		std::cerr << "Error: Failed to parse query from file: " << filePath << std::endl;
	}
}

void createSummaryTables(int M_input, int N_input)
{
	for (int h = 0; h <= 8; h++)
	{
		for (int seed = 1; seed <= 10; seed++)
		{
			deb("processing", seed, h);
			std::string queryPath = fmt::format("resource/sqls/_stocks_{}_{}_seeded_{}.spaql", N_input, M_input, seed);
			string tableForHardness = fmt::format("stocks_{}_{}_validate", N_input, M_input);

			auto spq = parseSpaqlFromFile(queryPath);
			deb(spq);
			if (spq)
			{
				spq->validate();
				unique_ptr<Stat> stat = std::make_unique<Stat>();
				size_t N = 10000;
				double E = 50;
				string originalTableName = spq->tableName;
				stat->analyze(tableForHardness);
				deb(tableForHardness);
				spq->setTableName(originalTableName);
				string tableForValidation = fmt::format("{}_validate", spq->tableName);
				stat->analyze(spq->tableName);
				deb(spq->tableName);
				stat->analyze(tableForValidation);
				deb(tableForValidation);
			}
		}
	}
}

void validation(string filepath, string outpath)
{
	// std::vector<std::string> headers = {"hardness", "p", "objective", "feas", "Runtime"};
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
		int h = record["Hardness"];
		// string dbPath = fmt::format("/home/fm2288/StochasticPackageQuery/test/QueriesExp2/{}/RS/{}/{}_{}.spaql", query, h, query, record["p"].dump());
		string dbPath = fmt::format("/home/fm2288/StochasticPackageQuery/test/Queries/{}/RS/{}_{}.spaql", query, query, h);
		// double p = (double)record["p"]/100.0;

		auto spq = parseSpaqlFromFile(dbPath);
		size_t N = 10000;
		double E = 50;
		// Bounder bounder(spq, N, E);
		// bounder.set(h);
		deb(spq);
		if (!spq)
		{
			std::cerr << "Error: Could not parse SPAQL file!" << std::endl;
			continue;
		}
		spq->validate();
		unique_ptr<Stat> stat = std::make_unique<Stat>();
		stat->analyze(spq->tableName);

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
		deb(res);
		vector<double> distances;
		vector<double> surpluses;
		bool feas = Check.feasible(res, distances, surpluses);
		double objective = Check.getObjective(res);
		std::cout << (feas ? "Validation Feasible Solution" : "Validation Infeasible Solution") << std::endl;
		std::cout << "Validation Objective = " << objective << std::endl;
		writer.addRow(h, objective, distances, record["Runtime"]);
	}
}

void generateSeededQueries(int N, int M)
{
	// Construct the experiment name and base output directory
	// Example: stocks_3_2_seeded
	std::string experimentName = fmt::format("stocks_{}_{}_seeded", N, M);
	std::string baseOutputDir = "/home/fm2288/StochasticPackageQuery/test/QueriesSeeds/" + experimentName;

	std::string rsBaseDir = baseOutputDir + "/RS";
	std::string rclBaseDir = baseOutputDir + "/RCL";

	// Create the base RS and RCL directories
	if (!create_nested_directories(rsBaseDir) || !create_nested_directories(rclBaseDir))
	{
		std::cerr << "Error: Failed to create base directories at " << baseOutputDir << std::endl;
		return;
	}

	std::cout << "Generating seeded queries in: " << baseOutputDir << std::endl;

	// Iterate through Hardness Levels (0 to 8)
	for (int h = 3; h <= 8; h++)
	{
		// Create hardness-specific subdirectories inside RS and RCL
		std::string rsHardnessDir = rsBaseDir + "/" + std::to_string(h);
		std::string rclHardnessDir = rclBaseDir + "/" + std::to_string(h);

		create_nested_directories(rsHardnessDir);
		create_nested_directories(rclHardnessDir);

		// Iterate through Seeds (1 to 10)
		for (int seed = 1; seed <= 10; seed++)
		{
			// Path to the specific seeded source file
			std::string queryPath = fmt::format("resource/sqls/_stocks_{}_{}_seeded_{}.spaql", N, M, seed);

			auto spq = parseSpaqlFromFile(queryPath);

			if (spq)
			{
				// --- 1. Apply Hardness Logic (Mirroring testRSSeed) ---
				Data::reset();
				Data::init(spq);
				Data::getInstance().fetchData();
				spq->validate();

				unique_ptr<Stat> stat = std::make_unique<Stat>();

				// Save the original seeded table name (e.g., stocks_3_2_seeded_1)
				string originalTableName = spq->tableName;
				// Switch to the validation table to calculate bounds/hardness
				string validationTableName = fmt::format("stocks_{}_{}_validate", N, M);

				stat->analyze(validationTableName);
				spq->setTableName(validationTableName);

				size_t N_bound = 10000;
				double E_bound = 50;
				Bounder bounder(spq, N_bound, E_bound);

				// Set the hardness on the validation table
				bounder.set(h);

				// Switch back to the seeded table name for the final query output
				spq->setTableName(originalTableName);

				// --- 2. Generate Query Strings ---

				// RS Query: The raw query with the new bounds applied
				std::string rsQueryString = static_cast<std::string>(*spq);

				// RCL Query: Transformed query (include count constraints, add REPEAT 0)
				// params: includeCountConstraint=true, addRepeatZero=true
				std::string rclQueryString = transform_query(rsQueryString, true, true);

				// --- 3. Save to Files ---
				// Filename: stocks_{N}_{M}_seeded_{seed}.spaql
				// The hardness is implicit in the directory structure: RS/h/filename
				std::string filename = fmt::format("stocks_{}_{}_seeded_{}.spaql", N, M, seed);

				// Write RS file
				std::string rsFullPath = rsHardnessDir + "/" + filename;
				std::ofstream rsOut(rsFullPath);
				if (rsOut)
				{
					rsOut << rsQueryString;
					rsOut.close();
				}
				else
				{
					std::cerr << "Error: Could not write to " << rsFullPath << std::endl;
				}

				// Write RCL file
				std::string rclFullPath = rclHardnessDir + "/" + filename;
				std::ofstream rclOut(rclFullPath);
				if (rclOut)
				{
					rclOut << rclQueryString;
					rclOut.close();
				}
				else
				{
					std::cerr << "Error: Could not write to " << rclFullPath << std::endl;
				}
			}
			else
			{
				std::cerr << "Warning: Could not parse source query " << queryPath << std::endl;
			}
		}
		std::cout << "Completed generation for Hardness level " << h << std::endl;
	}
}

int main(int argc, char *argv[])
{

	if (argc < 4)
	{
		std::cerr << "Usage: " << argv[0] << " <N> <M> <algorithm>" << std::endl;
		return 1;
	}

	int N = std::stoi(argv[1]);
	int M = std::stoi(argv[2]);
	std::string algorithm = argv[3];

	std::string dbPath = fmt::format("resource/sqls/_stocks_{}_{}.spaql", N, M);
	std::string outPath = fmt::format("stocks_{}_{}", N, M);

	if (algorithm == "Naive")
	{
		// testNaive(dbPath, M, outPath);
	}
	else if (algorithm == "SummarySearch")
	{
		throw std::runtime_error("SummarySearch testing is currently disabled. Please enable it when needed.");
	}
	else if (algorithm == "RS")
	{
		testRS(dbPath, M, outPath);
	}
	else if (algorithm == "RSExp2")
	{
		int h = std::stoi(argv[4]);
		string path = fmt::format("/home/fm2288/StochasticPackageQuery/test/QueriesExp2/{}", outPath);
		testRSExp2(path, M, h, outPath);
	}
	else if (algorithm == "RSSeed")
	{
		testRSSeed(M, N);
	}
	else if (algorithm == "generate")
	{
		generateQuerywithHardness(dbPath, outPath, N, M);
	}
	else if (algorithm == "genSeeds")
	{
		// New option to generate seeded queries with hardness structure
		generateSeededQueries(N, M);
	}
	else if (algorithm == "validate")
	{
		for (int h = 4; h <= 4; h++)
		{
			string path = fmt::format("/home/fm2288/StochasticPackageQuery/stochastic-sketchrefine/{}.json", outPath);
			validation(path, outPath);
		}
	}
	else if (algorithm == "genBounds")
	{
		generateQueryWithBounds(dbPath, outPath, N, M, "QueriesExp2", true, true);
		generateQueryWithBounds(dbPath, outPath, N, M, "QueriesExp3", false, false);
	}
	else if (algorithm == "summaryTables")
	{
		createSummaryTables(M, N);
	}
	else
	{
		cout << "Please enter a correct algorithm name next time" << endl;
	}

	return 0;
}