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

std::map<std::pair<int, int>, double> deter;

std::map<std::pair<int, int>, double> loadCSV(const std::string &filename)
{
	std::map<std::pair<int, int>, double> data;
	std::ifstream file(filename);
	std::string line;

	std::getline(file, line); // skip header

	while (std::getline(file, line))
	{
		std::stringstream ss(line);
		std::string token;
		int n, hardness;
		double objective;

		std::getline(ss, token, ',');
		n = std::stoi(token);
		std::getline(ss, token, ',');
		hardness = std::stoi(token);
		std::getline(ss, token, ',');
		objective = std::stod(token);

		data[{n, hardness}] = objective;
	}
	return data;
}

double getObjective(const std::map<std::pair<int, int>, double> &deter, int n, int hardness)
{
	auto it = deter.find({n, hardness});
	if (it != deter.end())
		return it->second;
	throw std::runtime_error("No entry found for N=" + std::to_string(n) +
							 ", Hardness=" + std::to_string(hardness));
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

void deterministicSolutionForHardness(int N)
{
	DataWriter writer("/home/fm2288/StochasticPackageQuery/test/Deterministic/deterministic_solutions.csv", {"N", "Hardness", "Objective", "surplus"});
	writer.writeHeaders();
	for (int h = 0; h <= 5; h++)
	{
		std::string queryPath = fmt::format("/home/fm2288/StochasticPackageQuery/test/QueriesSeeds/stocks_{}_{}_seeded/RS/{}/stocks_{}_{}_seeded_1.spaql", N, 10, h, N, 10);
		auto spq = parseSpaqlFromFile(queryPath);

		spq->setTableName(fmt::format("stocks_{}_validate", N));
		Data::reset();
		Data::init(spq);
		Data::getInstance().fetchData();
		string validateTable = fmt::format("stocks_{}", N);
		RobustSatisficing RS(10, spq, 0.46, validateTable);

		SolveOptions solveOptions;
		solveOptions.timeBudgeted = false;
		solveOptions.hyperParams = {{"reducedTuples", 500}, {"reducedScenarios", 100}, {"cap", 10}, {"benchmark", -1}};
		SolutionMetadata<int> sol = RS.solveDeterministic(spq, solveOptions);
		writer.addRow(N, h, sol.w, sol.bestRk);
	}
}

void testRSSeed(int M_input, int N_input, map<string, int> hyperParam = {{"reducedTuples", 512}, {"reducedScenarios", 50}, {"cap", 5}, {"benchmark", -1}}, 
										  map<string, bool> ablate = {{"stage1lcvar", false},{"stage1random", false}, {"stage1", false}, {"stage2", false}, {"ablate", false}})
{
	int M = M_input;
	int N = N_input;
	int M_hat = 1e6;
	deb(M,N);
	double epsilon = 0.46;

	std::string validateTable;
	std::string outPath;
	string output;
	std::vector<std::string> headers = {"Hardness", "Seed", "Objective", "ObjRatio", "deter_feas", "prob_feas", "feas", "surplus", "NTuples", "NScenarios", "Runtime", "solutions"};
	if (hyperParam["benchmark"] == -1 && !ablate["ablate"])
	{
		validateTable = fmt::format("stocks_{}", N);
		outPath = fmt::format("stocks_{}_{}", N, M);
		string dir = fmt::format("/home/fm2288/StochasticPackageQuery/test/RS_MAD/");
		create_nested_directories(dir);
		output = "/home/fm2288/StochasticPackageQuery/test/RS_MAD/RS" + outPath + ".csv";
	}
	else if (ablate["ablate"])
	{
		validateTable = fmt::format("stocks_{}", N);
		outPath = fmt::format("stocks_{}_{}", N, M);
		string study = ablate["stage2"] ? "stage2" : ablate["stage1lcvar"] ? "stage1lcvar" : ablate["stage1random"] ? "stage1random" : "stage1fixedN";
		string dir = fmt::format("/home/fm2288/StochasticPackageQuery/test/AblationStudy2/{}/", study);
		create_nested_directories(dir);
		output = dir + outPath + ".csv";
	}
	else
	{
		validateTable = fmt::format("stocks_{}", N);
		outPath = fmt::format("stocks_{}_{}", N, M);
		string dir = fmt::format("/home/fm2288/StochasticPackageQuery/test/HyperParameterBenchmark{}/N_{}_M_{}_cap_{}/", hyperParam["benchmark"], hyperParam["reducedTuples"], hyperParam["reducedScenarios"], hyperParam["cap"]);
		create_nested_directories(dir);
		output = fmt::format("/home/fm2288/StochasticPackageQuery/test/HyperParameterBenchmark{}/N_{}_M_{}_cap_{}/" + outPath + ".csv", hyperParam["benchmark"], hyperParam["reducedTuples"], hyperParam["reducedScenarios"], hyperParam["cap"], N, M);
	}
	DataWriter writer(output, headers);
	writer.writeHeaders();

	int max_concurrent = 15;
	if (N_input == 3)
	{
		max_concurrent = 30;
	}
	else if (N_input == 5 && M_input == 10000)
	{
		max_concurrent = 6;
	}
	int active_processes = 0;

	for (int h = 0; h <= 5; h++)
	{
		double deterObjective = getObjective(deter, N, h);
		for (int seed = 1; seed <= 10; seed = seed + 1)
		{
			// check if any existing children finished
			int status;
			while (waitpid(-1, &status, WNOHANG) > 0)
			{
				active_processes--;
			}

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
				double RSObjective;
				bool RSFeas;
				int Z;
				int q;
				int qScenarios;
				std::string queryPath = fmt::format("/home/fm2288/StochasticPackageQuery/test/QueriesSeeds/stocks_{}_{}_seeded/RS/{}/stocks_{}_{}_seeded_{}.spaql", N, M, h, N, M, seed);
				auto spq = parseSpaqlFromFile(queryPath);
				if (spq->validate())
				{
					double fetchingTime = 0.0;
					// we need to consider the fetching into effective runtime
					{
						gpro.clock("effectiveRuntime");
						gpro.clock("fetchingToEndOfFirstStage");
						Data::reset();
						Data::init(spq);
						Data::getInstance().fetchData();
						gpro.stop("effectiveRuntime");
						fetchingTime = gpro.getTime("effectiveRuntime");
						deb(fetchingTime);
					}
					gpro.clock("effectiveRuntime");
					RobustSatisficing RS(M, spq, epsilon, validateTable);

					double timeout_seconds = 10 * 60;

					// this is for the time budget experiment
					SolveOptions solveOptions;
					solveOptions.hyperParams = hyperParam;
					solveOptions.timeout_seconds = timeout_seconds;
					solveOptions.timeBudgeted = true;
					solveOptions.timeBudget = timeout_seconds;
					solveOptions.randomSeed = static_cast<unsigned int>(h * 10 + seed);
					SolutionMetadata<int> sol = RS.stochasticDualReducer(RS.spq, solveOptions, ablate);
					gpro.stop("effectiveRuntime");
					double totalTimeRS; 
					if(hyperParam["benchmark"] == 2)
					{
						totalTimeRS = gpro.getTime("effectiveRuntime") - gpro.getTime("fetchingToEndOfFirstStage");
					}else
					{
						totalTimeRS = gpro.getTime("effectiveRuntime");
					}

					{
						cout << "FINAL VALIDATION" << endl;
						spq->setTableName(validateTable);
						SPQChecker Check(spq);
						vector<double> feasibility;
						vector<double> surpluses;
						bool deterFeasible = false, probFeasible = false;
						if (sol.x.size() > 0)
						{
							// deb(sol.x);
							SolType res;
							for (int i = 0; i < RS.NTuples; i++)
							{
								if (sol.x[i] > 0)
								{
									res[i + 1] = sol.x[i];
								}
							}
							deb(res, res.size());
							RSFeas = Check.feasible(res, feasibility, surpluses, deterFeasible, probFeasible);
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
							double ratio = RSObjective / deterObjective;
							//{"Hardness", "Seed", "Objective", "ObjRatio", "deter_feas", "prob_feas", "feas", "surplus", "Runtime", "solutions"};
							for (int i = 0; i < surpluses.size(); i++)
							{
								if (surpluses[i] < 0)
								{
									probFeasible = false;
								}
							}
							writer.addRow(h, seed, RSObjective, ratio, deterFeasible, probFeasible, feasibility, surpluses, sol.qSz, sol.qScenarios, totalTimeRS, sol.solutions);

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

#include <fstream>
#include "nlohmann/json.hpp" // Make sure to have this header available

// For convenience
using json = nlohmann::json;

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

void runAllRSSeed()
{
	std::vector<int> Ns = {3, 4, 5};
	std::vector<int> Ms = {10, 100, 10000};
	for (int N : Ns)
	{
		for (int M : Ms)
		{
			testRSSeed(M, N);
		}
	}
}

void runAblationStudy()
{
	std::vector<int> Ns = {3,4,5};
	std::vector<int> Ms = {10,100,10000};
	for (int M : Ms)
	{
		for (int N : Ns)
		{
			testRSSeed(M, N, {{"reducedTuples", 512}, {"reducedScenarios", 50}, {"cap", 5}, {"benchmark", -1}}, {{"stage1lcvar", false},{"stage1random", false}, {"stage1", true}, {"stage2", false}, {"ablate", true}}); //original algorithm early stop in stage 2
			testRSSeed(M, N, {{"reducedTuples", 512}, {"reducedScenarios", 50}, {"cap", 5}, {"benchmark", -1}}, {{"stage1lcvar", true},{"stage1random", false}, {"stage1", true}, {"stage2", false}, {"ablate", true}}); //stage 1 with lcvar
			testRSSeed(M, N, {{"reducedTuples", 512}, {"reducedScenarios", 50}, {"cap", 5}, {"benchmark", -1}}, {{"stage1lcvar", false},{"stage1random", true}, {"stage1", true}, {"stage2", false}, {"ablate", true}}); //random in stage 1
			testRSSeed(M, N, {{"reducedTuples", 512}, {"reducedScenarios", 50}, {"cap", 5}, {"benchmark", -1}}, {{"stage1lcvar", false},{"stage1random", false}, {"stage1", false}, {"stage2", true}, {"ablate", true}}); //stage 2 it will runNaiveFinalStage
			testRSSeed(M, N); //full feedback original algorithm			
		}
	}
}

void runBenchmarkScenariosCap(int N_input)
{
	std::vector<int> reducedTuples = {128};
	std::vector<int> reducedScenariosValues = {10, 25, 50, 100};
	std::vector<int> capValues = {1, 3, 5, 10};
	std::vector<int> M = {10, 100, 10000};

	for (int m : M)
	{
		// Iterate through all combinations of hyperparameters
		for (int reducedTuples : reducedTuples)
		{
			for (int reducedScenarios : reducedScenariosValues)
			{
				for (int cap : capValues)
				{
					std::map<std::string, int> hyperParams = {
						{"reducedTuples", reducedTuples},
						{"reducedScenarios", reducedScenarios},
						{"cap", cap},
						{"benchmark", 1}};
					testRSSeed(m, N_input, hyperParams);
				}
			}
		}
	}

	std::cout << "Completed hyperparameter benchmark." << std::endl;
}

void allBenchmarkTuples()
{
	vector<int> Ns = {5};
	vector<int> Ms = {10, 100, 10000};

	for (int m : Ms)
	{
		for (int n : Ns)
		{
			for (int i = 7; i <= 16; i++)
			{
				int reducedTuples = pow(2, i);
				std::map<std::string, int> hyperParams = {
					{"reducedTuples", reducedTuples},
					{"reducedScenarios", 50},
					{"cap", 5},
					{"benchmark", 2}};
				testRSSeed(m, n, hyperParams);
			}
		}
	}
}

void allBenchmarkScenariosCap()
{
	std::vector<int> Ns = {3, 4};
	for (int N : Ns)
	{
		runBenchmarkScenariosCap(N);
	}
}

void createSummaryTables(int N_input)
{
	string tableForHardness = fmt::format("stocks_{}_validate", N_input);
	unique_ptr<Stat> stat = std::make_unique<Stat>();
	size_t N = 10000;
	double E = 50;
	deb(tableForHardness);
	stat->analyze(tableForHardness);
}

void generateSeededQueries(int N, int M)
{
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
	for (int h = 0; h <= 5; h++)
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
				spq->validate();
				deb("validated");
				unique_ptr<Stat> stat = std::make_unique<Stat>();

				// Save the original seeded table name (e.g., stocks_3_2_seeded_1)
				string originalTableName = spq->tableName;
				// Switch to the validation table to calculate bounds/hardness
				string validationTableName = fmt::format("stocks_{}_validate", N);

				stat->analyze(validationTableName);
				spq->setTableName(validationTableName);

				size_t N_bound = 10000;
				double E_bound = 50;
				Bounder bounder(spq, N_bound, E_bound);

				// Set the hardness on the validation table
				deb("bounding");
				bounder.set(h);
				deb("bounded");

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
	deter = loadCSV("/home/fm2288/StochasticPackageQuery/test/Deterministic/deterministic_solutions.csv");
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

	if (algorithm == "allRSSeed")
	{
		runAllRSSeed();
	}
	else if (algorithm == "RSSeed")
	{
		testRSSeed(M, N);
	}
	else if (algorithm == "deter")
	{
		deterministicSolutionForHardness(N);
	}
	else if (algorithm == "genSeeds")
	{
		// New option to generate seeded queries with hardness structure
		generateSeededQueries(N, M);
	}
	else if (algorithm == "allBenchmarkScenariosCap")
	{
		allBenchmarkScenariosCap();
	}
	else if (algorithm == "allBenchmarkTuples")
	{
		allBenchmarkTuples();
	}
	else if (algorithm == "summaryTables")
	{
		createSummaryTables(N);
	}
	else if (algorithm == "ablation")
	{
		runAblationStudy();
	}
	else
	{
		cout << "Please enter a correct algorithm name next time" << endl;
	}

	return 0;
}