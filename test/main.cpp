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
		size_t N = 100;
		double E = 10.5;
		Bounder bounder (spq, N, E);
		// spq->setVariable("l1", 10);
		// cout << spq;
	}
	STOP(pro);
	PRINT(pro);
}

int main() {
	// analyzeAll();
	test();
}