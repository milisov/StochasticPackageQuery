#ifndef STAT_HPP
#define STAT_HPP

#include <memory>

#include "util/uio.hpp"

using std::unique_ptr;

class Stat{
private:
    static bool init();
    static bool rebuildStat;
    static string statTable;
    unique_ptr<PgManager> pg;
public:
    Stat();
};

#endif