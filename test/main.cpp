#include <iostream>

#include "util/udebug.hpp"
#include "spq/spq.hpp"
#include "spq/stat.hpp"

// #include <pcg_random.hpp>
// #include <gfx/timsort.hpp>

int main() {
	INIT(pro);
	CLOCK(pro);
	string filePath = "resource/sqls/_stocks_30_100.spaql";
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