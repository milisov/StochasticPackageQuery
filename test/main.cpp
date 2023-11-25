#include <iostream>

#include "util/udebug.hpp"
#include "spq/spq.hpp"
#include "spq/stat.hpp"
#include "spq/bounder.hpp"

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

void test(){
	INIT(pro);
	CLOCK(pro);
	string filePath = "resource/sqls/_stocks_3000_100.spaql";
	auto spq = parseSpaqlFromFile(filePath);
	if (spq){
		cout << "Success!\n" << spq;
		deb(spq->validate());
		unique_ptr<Stat> stat = std::make_unique<Stat>();
		stat->analyze(spq->tableName);
		size_t N = 1000;
		double E = 150;
		CLK(pro, "initBounder");
		Bounder bounder (spq, N, E);
		STP(pro, "initBounder");
		vector<double> hards;
		for (double i = -10; i <= 10; i ++) hards.push_back(i);
		CLK(pro, "hard");
		bounder.generate(hards);
		bounder.set(1);
		STP(pro, "hard");
		deb(spq->executable(), spq);
	}
	STOP(pro);
	PRINT(pro);
}

int main() {
	// analyzeAll();
	test();
}