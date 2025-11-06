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

using std::map;
using std::vector;
using json = nlohmann::json;



void testSDR(string path, int M_input, string outPath)
{
	int M = M_input;
	int M_hat = 1e6;
	double epsilon = 0.46;

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
		// N is number of scenarios in order to approx, E - expected package size in sol
		// set values that are variables in the query
		Bounder bounder(spq, N, E);
		deb(spq->executable(), spq);
		std::vector<std::string> headers = {"Hardness", "SDR-objective", "SDR-Feas", "Z", "Q", "RuntimeSDR"};
		string output = "/home/fm2288/StochasticPackageQuery/test/Experiments2/SDR/SDR" + outPath + ".csv";
		DataWriter writer(output, headers);

		Profiler stopwatchSS;
      
		for (int h = 0; h <= 0; h = h + 1)
		{
			double SDRObjective;
			bool SDRFeas;
			int Z;
			bounder.set(h);
			deb(spq);
			StochDualReducer SDR(M, spq, epsilon);

			string labelSDR = "SDR with Hardness" + std::to_string(h);

			auto start_time = std::chrono::steady_clock::now();
			double timeout_seconds = 600;
			FormulateOptions formOptions;
			formOptions.reduced = false;
			formOptions.qSz = 500;
			formOptions.Z = 1;
			stopwatchSS.clock(labelSDR);
			SolutionMetadata<int> sol = SDR.stochDualReducer(SDR.spq, formOptions, curveFitOptions, start_time, timeout_seconds);
			stopwatchSS.stop(labelSDR);
			double totalTimeSDR = stopwatchSS.getTime(labelSDR);

			SPQChecker Check(SDR.spq);
			double distance;
			if (sol.x.size() > 0)
			{
				deb(sol.x, sol.x.size());
				SolType res;
				for (int i = 0; i < SDR.NTuples; i++)
				{
					res[i + 1] = sol.x[i];
				}
				deb(SDR.NTuples);
				SDRFeas = Check.feasible(res, distance);
				SDRObjective = Check.getObjective(res); // we need to change this to validation objective
				Z = sol.Z;
			}
			else
			{
				SDRObjective = -1;
				SDRFeas = 0;
				Z = sol.Z;
			}
			writer.addRow(h, SDRObjective, distance, Z, sol.qSz, totalTimeSDR);
		}
	}
}

void testSummarySearch(string path, int M_input, string outPath)
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
		// N is number of scenarios in order to approx, E - expected package size in sol
		// set values that are variables in the query
		Bounder bounder(spq, N, E);
		deb(spq->executable(), spq);
		std::vector<std::string> headers = {"Hardness", "SS-objective", "SS-Feas", "Z", "RuntimeSummarySearch"};
		string output = "/home/fm2288/StochasticPackageQuery/test/Experiments2/Deter/Deter" + outPath + ".csv";
		DataWriter writer(output, headers);

		Profiler stopwatchSS;

		vector<int> reducedIds;
		bool reduced = false;
		for (int h = -4; h <= 8; h = h + 1)
		{
			double SSObjective;
			bool SSFeas;
			int Z;
			bounder.set(h);
			deb(spq);
			SummarySearch SS(M, spq, epsilon);

			string labelSS = "SummarySearch with Hardness" + std::to_string(h);
			stopwatchSS.clock(labelSS);
			SSFormulator formulator(SS.spq);
			FormulateOptions formOptions;
			formOptions.reduced = false;
			formOptions.Z = 1;
			DecisionVarOptions decVarOptions;
			setDecisionVarOptions(decVarOptions, 0.0, 1.0, 0.0, GrbVarType::Integer);
			formOptions.decisionVarOptions = decVarOptions;

			auto start_time = std::chrono::steady_clock::now();
			double timeout_seconds = 600;
			SolutionMetadata<int> sol = SS.summarySearch<int>(SS.spq, formulator, formOptions, curveFitOptions, 1, start_time, timeout_seconds);
			stopwatchSS.stop(labelSS);
			double totalTimeSS = stopwatchSS.getTime(labelSS);

			SPQChecker Check(SS.spq);
			double distance;
			if (sol.x.size() > 0)
			{
				deb(sol.x);
				SolType res;
				for (int i = 0; i < SS.NTuples; i++)
				{
					res[i + 1] = sol.x[i];
				}
				SSFeas = Check.feasible(res, distance);
				SSObjective = Check.getObjective(res); // we need to change this to validation objective
				Z = sol.Z;
			}
			else
			{
				SSObjective = -1;
				SSFeas = 0;
				Z = sol.Z;
			}
			writer.addRow(h, SSObjective, distance, Z, totalTimeSS);
		}
	}
}


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
		// deb(spq->executable(), spq);
		std::vector<std::string> headers = {"Hardness", "RS-objective", "RS-Feas", "bestEps", "RuntimeRS"};
		string output = "/home/fm2288/StochasticPackageQuery/test/Experiments2/RobustSatisficing/RS" + outPath + ".csv";
		DataWriter writer(output, headers);

		Profiler stopwatchRS;

		vector<int> reducedIds;
		bool reduced = false;

		double RSObjective;
		bool RSFeas;
		int Z;
		int q;
		for (int h = 0; h <= 8; h = h + 1)
		{
			double RSObjective;
			bool RSFeas;
			int Z;
			bounder.set(h);
			RobustSatisficing RS(M, spq, epsilon);

			string labelRS = "Robust Satisficing with Hardness" + std::to_string(h);
			stopwatchRS.clock(labelRS);

			SolutionMetadata<int> sol = RS.stochasticDualReducer(RS.spq, curveFitOptions);
			stopwatchRS.stop(labelRS);
			double totalTimeRS = stopwatchRS.getTime(labelRS);
			double timeStage1 = sol.timeStage1;
			double timeStage2 = sol.timeStage2;
			SPQChecker Check(RS.spq);
			double feasScore;
			if (sol.x.size() > 0)
			{
				deb(sol.x);
				SolType res;
				for (int i = 0; i < RS.NTuples; i++)
				{
					res[i + 1] = sol.x[i];
				}
				RSFeas = Check.feasible(res, feasScore);
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
			writer.addRow(h, RSObjective, feasScore, sol.bestEps, totalTimeRS);
		}
	}
}

#include <fstream>
#include "nlohmann/json.hpp" // Make sure to have this header available

// For convenience
using json = nlohmann::json;

void testRSExp2(const std::string &basePath, int M_input, int hardness, const std::string &outPath)
{
	for(int h = hardness; h <= 4; h++)
	{
		// --- 1. Construct the target directory path to read queries from ---
		std::string targetDirectory = basePath + "/RS/" + std::to_string(h);
		std::cout << "Reading queries from directory: " << targetDirectory << std::endl;

		// --- 2. Set up the single CSV writer for all results ---
		std::vector<std::string> headers = {"Hardness", "p", "RS-objective", "RS-Feas", "bestEps", "RuntimeRS"};
		std::string outputCsvFile = "/home/fm2288/StochasticPackageQuery/test/ExperimentsAUC/RobustSatisficing/RS" + outPath + ".csv";
		DataWriter writer(outputCsvFile, headers);
		std::cout << "Writing results to: " << outputCsvFile << std::endl;

		// --- 3. Open and iterate through the target directory ---
		DIR *dir = opendir(targetDirectory.c_str());
		if (dir == nullptr)
		{
			std::cerr << "Error: Could not open directory " << targetDirectory << std::endl;
			return;
		}

		struct dirent *entry;
		while ((entry = readdir(dir)) != nullptr)
		{
			std::string filename = entry->d_name;

			// Process only files that end with ".spaql"
			if (filename.length() <= 6 || filename.substr(filename.length() - 6) != ".spaql")
			{
				continue;
			}

			std::string fullQueryPath = targetDirectory + "/" + filename;
			std::cout << "\nProcessing file: " << fullQueryPath << std::endl;
			double p;
			try
			{
				size_t last_underscore = filename.rfind('_');
				size_t dot_pos = filename.rfind(".spaql");
				if (last_underscore == std::string::npos || dot_pos == std::string::npos)
				{
					throw std::invalid_argument("Filename format incorrect");
				}
				std::string p_str = filename.substr(last_underscore + 1, dot_pos - (last_underscore + 1));
				p = std::stod(p_str) / 100.0;
			}
			catch (const std::exception &e)
			{
				std::cerr << "Warning: Could not parse probability from filename '" << filename << "'. Skipping." << std::endl;
				continue;
			}

			// --- 5. Core Solving Logic (applied to each query file) ---
			auto spq = parseSpaqlFromFile(fullQueryPath);
			deb(spq);
			if (!spq)
			{
				std::cerr << "Warning: Failed to parse query from file '" << fullQueryPath << "'. Skipping." << std::endl;
				continue;
			}

			// Initialize data and validate the parsed query
			Data::init(spq);
			Data::getInstance().fetchData();
			spq->validate();

			int M = M_input;
			double epsilon = 0.46;
			map<string, Option> curveFitOptions = {{"arctan", false}, {"binarySearch", true}};

			RobustSatisficing RS(M, spq, epsilon);
			Profiler stopwatchRS;
			string labelRS = "RobustSatisficing_h" + std::to_string(h) + "_p" + std::to_string((int)(p * 100));

			stopwatchRS.clock(labelRS);
			SolutionMetadata<int> sol = RS.stochasticDualReducer(RS.spq, curveFitOptions);
			stopwatchRS.stop(labelRS);
			double totalTimeRS = stopwatchRS.getTime(labelRS);

			SPQChecker Check(RS.spq);
			double RSObjective = -1.0;
			double feasScore = 0.0;
			bool RSFeas = false;

			if (!sol.x.empty())
			{
				SolType res;
				for (size_t i = 0; i < sol.x.size(); ++i)
				{
					res[i + 1] = sol.x[i];
				}
				deb(res);
				RSFeas = Check.feasible(res, feasScore);
				RSObjective = Check.getObjective(res);

				// --- 6. Save results to a JSON file ---
				json j;
				j["basePath"] = basePath;
				j["hardness"] = h;
				j["p"] = p;

				// Convert the SolType map to a JSON object
				// The nlohmann/json library can't directly handle a map with integer keys,
				// so we convert the keys to strings.
				json res_json;
				for (const auto& pair : res) {
					res_json[std::to_string(pair.first)] = pair.second;
				}
				j["res"] = res_json;

				// Construct a unique filename for the JSON output
				std::string jsonOutputFilename = "RS_solution_h" + std::to_string(h) + "_p" + std::to_string(static_cast<int>(p * 100)) + ".json";
				std::string jsonOutputPath = "/home/fm2288/StochasticPackageQuery/test/ExperimentsAUC/RobustSatisficing/Solutions/" + jsonOutputFilename;

				// Write the JSON object to a file
				std::ofstream o(jsonOutputPath);
				o << std::setw(4) << j << std::endl;
			}
			writer.addRow(h, p, RSObjective, feasScore, sol.bestEps, totalTimeRS);
		}

		closedir(dir);
		std::cout << "\nBatch processing complete." << std::endl;
	}
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
    const std::string &original_query,\
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
    const std::string& experimentName, // New parameter for the experiment directory
    bool includeCountConstraint,       // Flag to control the count constraint in RS
    bool addRepeatZero                 // Flag to control "REPEAT 0" in RCL
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
            if (!includeCountConstraint) {
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

void generateQuerywithHardness(const std::string& path, const std::string& outPath, int n, int m)
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
		mkdir(rsDir.c_str(), 0775);         // Create the RS subdirectory
		mkdir(rclDir.c_str(), 0775);        // Create the RCL subdirectory

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

void validation(string filepath, string outpath)
{
	//std::vector<std::string> headers = {"hardness", "p", "objective", "feas", "Runtime"};
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
		//string dbPath = fmt::format("/home/fm2288/StochasticPackageQuery/test/QueriesExp2/{}/RS/{}/{}_{}.spaql", query, h, query, record["p"].dump());
		string dbPath = fmt::format("/home/fm2288/StochasticPackageQuery/test/Queries/{}/RS/{}_{}.spaql", query, query, h);
		//double p = (double)record["p"]/100.0;
		
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
		testSummarySearch(dbPath, M, outPath);
	}
	else if (algorithm == "SDR")
	{
		testSDR(dbPath, M, outPath);
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
	else if (algorithm == "generate")
	{
		generateQuerywithHardness(dbPath, outPath, N, M);
	}
	else if (algorithm == "validate")
	{
		for(int h = 4; h <= 4; h++)
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
	else
	{
		cout << "Please enter a correct algorithm name next time" << endl;
	}

	return 0;
}