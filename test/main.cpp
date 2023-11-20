#include <iostream>

#include "util/udebug.hpp"
#include "util/uio.hpp"

int main() {
	string filePath = "resource/sqls/ebnf.spaql";
	auto spq = parseSpaqlFromFile(filePath);
	if (spq){
		// cout << "Success!\n" << spq;
		spq->validate();
	}
}