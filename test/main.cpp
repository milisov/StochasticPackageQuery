#include <iostream>

#include "util/udebug.hpp"
#include "spq/spq.hpp"
#include "core/stat.hpp"
#include "spq/bounder.hpp"
#include "core/checker.hpp"

void analyzeAll(){
	vector<int> nStocks = {3, 30, 300, 3000};
	vector<int> nPaths = {100, 10000};
	for (auto nStock : nStocks){
		for (auto nPath : nPaths){
			INIT(pro);
			CLOCK(pro);
			string filePath = fmt::format("resource/sqls/_stocks_{}_{}.spaql", nStock, nPath);
			auto spq = parseSpaqlFromFile(filePath);
			if (spq){
				cout << "Success!\n" << spq;
				deb(spq->validate());
				unique_ptr<Stat> stat = std::make_unique<Stat>();
				stat->analyze(spq->tableName);
			}
			STOP(pro);
			PRINT(pro);
		}
	}
}

#include "solver/taylor.hpp"

void test(){
	INIT(pro);
	CLOCK(pro);
	string filePath = "resource/sqls/_stocks_30_100.spaql";
	auto spq = parseSpaqlFromFile(filePath);
	if (spq){
		cout << "Success!\n" << spq;
		deb(spq->validate());
		unique_ptr<Stat> stat = std::make_unique<Stat>();
		stat->analyze(spq->tableName);
		size_t N = 1000;
		double E = 50;
		CLK(pro, "initBounder");
		Bounder bounder (spq, N, E);
		STP(pro, "initBounder");
		vector<double> hards;
		// for (double i = -10; i <= 10; i ++) hards.push_back(i);
		CLK(pro, "hard");
		// bounder.generate(hards);
		bounder.set(3);
		STP(pro, "hard");
		deb(spq->executable(), spq);
		CLK(pro, "taylorinit");
		Taylor taylor (spq, {}, {
			{"soft_deterministic_constraint", false},
			{"max_number_of_iterations", 100}});
		STP(pro, "taylorinit");
		CLK(pro, "taylor");
		taylor.solve();
		STP(pro, "taylor");
	}
	STOP(pro);
	// PRINT(pro);
}

#include "oneapi/tbb/concurrent_unordered_map.h"
using oneapi::tbb::concurrent_unordered_map;
#include <pcg_random.hpp>
#include "util/uconfig.hpp"

void testTBB(){
	// concurrent_unordered_map<int, int> cmap;
	// cmap[0] = 1;
	// deb(cmap);
	// pcg32 gen (Config::getInstance()->seed());
	// vector<unsigned int> seeds;
	// for (int i = 0; i < 80; i ++) seeds.push_back(gen());
	// deb(seeds);
}

#include <Highs.h>
#include <vector>

void testHighs(){
	HighsModel model;
	model.lp_.num_col_ = 1;
	model.lp_.num_col_++;
	model.lp_.num_row_ = 3;
	model.lp_.sense_ = ObjSense::kMinimize;
	model.lp_.offset_ = 3;
	model.lp_.col_cost_ = {1.0};
	model.lp_.col_cost_.push_back(1.0);
	model.lp_.col_lower_ = {0.0, 1.0};
	model.lp_.col_upper_ = {4.0, 1.0e30};
	model.lp_.row_lower_ = {-1.0e30, 5.0, 6.0};
	model.lp_.row_upper_ = {7.0, 15.0, 1.0e30};
	//
	// Here the orientation of the matrix is column-wise
	model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
	// a_start_ has num_col_1 entries, and the last entry is the number
	// of nonzeros in A, allowing the number of nonzeros in the last
	// column to be defined
	model.lp_.a_matrix_.start_ = {0, 2, 5};
	model.lp_.a_matrix_.index_ = {1, 2, 0, 1, 2};
	model.lp_.a_matrix_.value_ = {1.0, 3.0, 1.0, 2.0, 2.0};
	//
	// Create a Highs instance
	Highs highs;
	HighsStatus return_status;
	// highs.setOptionValue("output_flag", false); // This should turn off logging output
	// highs.setOptionValue("log_to_console", false); // This should turn off logging output
	highs.setOptionValue("random_seed", abs(static_cast<int>(Config::getInstance()->seed()))); // This should turn off logging output
	highs.setOptionValue("simplex_strategy", 3); // This should turn off logging output
	highs.setOptionValue("infinite_bound", POS_INF);
	//
	// Pass the model to HiGHS
	return_status = highs.passModel(model);
	assert(return_status==HighsStatus::kOk);
	// If a user passes a model with entries in
	// model.lp_.a_matrix_.value_ less than (the option)
	// small_matrix_value in magnitude, they will be ignored. A logging
	// message will indicate this, and passModel will return
	// HighsStatus::kWarning
	//
	// Get a const reference to the LP data in HiGHS
	const HighsLp& lp = highs.getLp();
	// Solve the model
	return_status = highs.run();
	assert(return_status==HighsStatus::kOk);
	//
	// Get the model status
	const HighsModelStatus& model_status = highs.getModelStatus();
	assert(model_status==HighsModelStatus::kOptimal);
	cout << "Model status: " << highs.modelStatusToString(model_status) << endl;
	//
	// Get the solution information
	const HighsInfo& info = highs.getInfo();
	cout << "Simplex iteration count: " << info.simplex_iteration_count << endl;
	cout << "Objective function value: " << info.objective_function_value << endl;
	cout << "Primal  solution status: " << highs.solutionStatusToString(info.primal_solution_status) << endl;
	cout << "Dual    solution status: " << highs.solutionStatusToString(info.dual_solution_status) << endl;
	cout << "Basis: " << highs.basisValidityToString(info.basis_validity) << endl;
	const bool has_values = info.primal_solution_status;
	const bool has_duals = info.dual_solution_status;
	const bool has_basis = info.basis_validity;
	//
	// Get the solution values and basis
	const HighsSolution& solution = highs.getSolution();
	const HighsBasis& basis = highs.getBasis();
	//
	// Report the primal and solution values and basis
	for (int col=0; col < lp.num_col_; col++) {
		cout << "Column " << col;
		if (has_values) cout << "; value = " << solution.col_value[col];
		if (has_duals) cout << "; dual = " << solution.col_dual[col];
		if (has_basis) cout << "; status: " << highs.basisStatusToString(basis.col_status[col]);
		cout << endl;
	}
	for (int row=0; row < lp.num_row_; row++) {
		cout << "Row    " << row;
		if (has_values) cout << "; value = " << solution.row_value[row];
		if (has_duals) cout << "; dual = " << solution.row_dual[row];
		if (has_basis) cout << "; status: " << highs.basisStatusToString(basis.row_status[row]);
	
		cout << endl;
	}
}

void testOmp(){
	size_t N = 10000000;
	size_t R = 10;
	vector<double> v (N);
	std::iota(v.begin(), v.end(), 0);
	INIT(pro);
	for (size_t r = 0; r < R; ++r){
		CLK(pro, "a");
		double norm = 0;
		vector<double> vv (N);
		for (size_t i = 0; i < N; ++i){
			norm += v[i]*v[i];
		}
		norm = sqrt(norm);
		for (size_t i = 0; i < N; ++i){
			vv[i] = v[i]/norm;
		}
		STP(pro, "a");
	}
	int core = 80;
	for (size_t r = 0; r < R; ++r){
		CLK(pro, "b");
		vector<double> vv (N);
		double norm;
		#pragma omp parallel num_threads(core)
		{
			double norm_ = 0;
			#pragma omp for nowait
			for (size_t i = 0; i < N; ++i){
				norm_ += v[i]*v[i];
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
			for (size_t i = 0; i < N; ++i){
				vv[i] = v[i]/norm;
			}
		}
		STP(pro, "b");
	}
	PRINT(pro);
}

void testNumeric(){
	vector<double> v = {0, 0};
	deb(normalize(v));
	deb(v);
}

int main() {
	// testNumeric();
	// analyzeAll();
	test();
	// testOmp();
	// testTBB();
	// testHighs();
	// map<string, Option> ok = {{"ok","1"},{"ok1",1.5},{"ok2",false}};
	// deb(ok);
	// for (auto p : ok) deb(p.second.which());
}