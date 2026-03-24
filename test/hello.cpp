#include <iostream>

#include "util/udebug.hpp"
#include "spq/spq.hpp"
#include "core/stat.hpp"
#include "spq/bounder.hpp"
#include "core/checker.hpp"
#include "util/uio.hpp"
#include "solver/starter.hpp"
#include "dlib/optimization.h"
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <gurobi_c.h>
#include "util/uio.hpp"
#include <fstream>
#include <string>

void dropSeeded()
{
    deb("Deleting Tables...");
    try
    {
        PgManager pg; // connect to DB
        for (int i = 1; i <= 10; i++)
        {
            std::string tableName = fmt::format("stocks_3_2_seeded_{}", i);

            std::string sql = fmt::format(R"(DROP TABLE IF EXISTS {})", tableName);

            auto res = PQexec(pg.conn.get(), sql.c_str());
            ck(pg.conn, res);
            PQclear(res);

            std::cout << "Table '" << tableName << "' dropped.\n";
        }

        for (int i = 1; i <= 10; i++)
        {
            std::string tableName = fmt::format("stocks_3_2_seeded_{}_summary", i);

            std::string sql = fmt::format(R"(DROP TABLE IF EXISTS {})", tableName);

            auto res = PQexec(pg.conn.get(), sql.c_str());
            ck(pg.conn, res);
            PQclear(res);

            std::cout << "Table '" << tableName << "' dropped.\n";
        }

        for (int i = 1; i <= 10; ++i)
        {
            std::string tableName = fmt::format("stocks_3_2_validate_seeded_{}", i);

            std::string sql = fmt::format(R"(DROP TABLE IF EXISTS {})", tableName);

            auto res = PQexec(pg.conn.get(), sql.c_str());
            ck(pg.conn, res);
            PQclear(res);

            std::cout << "Table '" << tableName << "' dropped.\n";
        }
        for (int i = 1; i <= 10; ++i)
        {
            std::string tableName = fmt::format("stocks_3_2_seeded_{}_validate_summary", i);

            std::string sql = fmt::format(R"(DROP TABLE IF EXISTS {})", tableName);

            auto res = PQexec(pg.conn.get(), sql.c_str());
            ck(pg.conn, res);
            PQclear(res);

            std::cout << "Table '" << tableName << "' dropped.\n";
        }

    }
    catch (const std::exception &e)
    {
        std::cerr << "Exception: " << e.what() << "\n";
    }    
}

void createSeededTables()
{
    deb("Creating seeded tables...");
    try
    {
        PgManager pg; // connect to DB

        for (int i = 1; i <= 10; ++i)
        {
            std::string tableName = fmt::format("stocks_3_2_seeded_{}", i);
            std::string seedStr = fmt::format("-seed-{}", i);
            std::string tableToSeedFrom = "stocks_3_4";
            std::string deterAttrTable = "stocks_3_2";

            std::string sql = fmt::format(R"(
                CREATE TABLE IF NOT EXISTS {} AS
                WITH all_profits AS (
                    SELECT unnest(profit) AS p
                    FROM {}
                ),
                sampled_profit AS (
                    SELECT array_agg(p ORDER BY md5(p::text || '{}')) AS full_array
                    FROM all_profits
                )
                SELECT s.id, s.stock, s.day, s.price,
                       full_array[1:100] AS profit
                FROM {} s
                CROSS JOIN sampled_profit;
            )",tableName, tableToSeedFrom, seedStr, deterAttrTable);

            auto res = PQexec(pg.conn.get(), sql.c_str());
            ck(pg.conn, res);
            PQclear(res);

            std::cout << "Table '" << tableName << "' created successfully.\n";
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Exception: " << e.what() << "\n";
    }
}

void createValidateTablesForSeeded()
{
    deb("Creating seeded tables...");
    try
    {
        PgManager pg; // connect to DB

        for (int i = 1; i <= 10; ++i)
        {
            std::string tableName = fmt::format("stocks_3_2_validate", i);
            std::string validTableName = fmt::format("stocks_3_2_seeded_{}_validate", i);

            std::string sql = fmt::format(R"(CREATE TABLE IF NOT EXISTS {} AS SELECT * FROM {})",validTableName, tableName);

            auto res = PQexec(pg.conn.get(), sql.c_str());
            ck(pg.conn, res);
            PQclear(res);

            std::cout << "Table '" << validTableName << "' created successfully.\n";
        }
        for (int i = 1; i <= 10; ++i)
        {
            std::string tableName = fmt::format("stocks_3_2_validate_summary", i);
            std::string validTableName = fmt::format("stocks_3_2_seeded_{}_validate_summary", i);

            std::string sql = fmt::format(R"(CREATE TABLE IF NOT EXISTS {} AS SELECT * FROM {})",validTableName, tableName);

            auto res = PQexec(pg.conn.get(), sql.c_str());
            ck(pg.conn, res);
            PQclear(res);

            std::cout << "Table '" << validTableName << "' created successfully.\n";
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Exception: " << e.what() << "\n";
    }    
}

void createQueriesFiles()
{
    std::string baseTableName = "stocks_3_2";
    std::string directory = "/home/fm2288/StochasticPackageQuery/resource/sqls/seeded/"; 

    for (int i = 1; i <= 10; ++i)
    {
        std::string tableName = baseTableName + "_seeded_" + std::to_string(i);
        std::string fileName = directory + tableName + ".spaql";

        std::ofstream outFile(fileName);

        if (outFile.is_open())
        {
            outFile << "SELECT PACKAGE(*) FROM " << tableName << " SUCH THAT\n";
            outFile << "    COUNT(*) BETWEEN l1 AND u1 AND\n";
            outFile << "    SUM(price) <= u2 AND\n";
            outFile << "    SUM(profit) >= v1 WITH PROBABILITY >= 0.95\n";
            outFile << "MAXIMIZE EXPECTED SUM(profit)";

            outFile.close();
        }
    }
}

bool endsWith(const std::string& str, const std::string& suffix) {
    if (str.length() < suffix.length()) {
        return false;
    }
    return str.compare(str.length() - suffix.length(), suffix.length(), suffix) == 0;
}

void analyzeAllTables()
{
    deb("Starting bulk analysis of tables...");
    try
    {
        PgManager pg;
        unique_ptr<Stat> stat = std::make_unique<Stat>();

        std::string sql = R"(
            SELECT table_name 
            FROM information_schema.tables 
            WHERE table_schema = 'public' 
            AND table_type = 'BASE TABLE'
        )";

        auto res = PQexec(pg.conn.get(), sql.c_str());
        ck(pg.conn, res);

        int rows = PQntuples(res);
        std::cout << fmt::format("Found {} tables in the database.\n", rows);

        for (int i = 0; i < rows; ++i)
        {
            std::string tableName = PQgetvalue(res, i, 0);
            std::cout << fmt::format("Processing table: '{}'\n", tableName);
            // 1. Ignore if it is a summary table
            if (endsWith(tableName, "_summary")) {
                std::cout << fmt::format("Skipping summary table: {}\n", tableName);
                continue;
            }

            // 2. Ignore the metadata stat table itself (usually named 'stat' or defined in config)
            // You might want to adjust this name based on your Config::getInstance()->pt.get<string>("pgmanager.stat_table")
            if (tableName == "stat" || tableName == "pg_stat_statements") {
                continue;
            }

            std::cout << fmt::format("Processing table: '{}' ... ", tableName);

            // 3. Call Analyze
            // Note: This relies on your 'uconfig' file having "partition.<tableName>" entries.
            // If the entry is missing, stat.analyze will usually print an error and return false.
            bool result = stat->analyze(tableName);

            if (result) {
                std::cout << "Analysis Complete.\n";
            } else {
                // Determine if it failed because of missing config or empty table
                std::cout << "Skipped (See error above, likely missing config).\n";
            }
        }

        PQclear(res);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Exception: " << e.what() << "\n";
    }
}




int main()
{
    //dropSeeded();
    analyzeAllTables();
    //createQueriesFiles();
    //createValidateTablesForSeeded();
}
