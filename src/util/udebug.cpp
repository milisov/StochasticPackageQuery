#include "udebug.hpp"

#include <omp.h>
#include <fmt/core.h>

const char* RED = "\033[31m";
const char* GREEN = "\033[32m";
const char* RESET = "\033[0m";

Profiler::Profiler(){
}

void Profiler::clock(const string& label){
    timePoints[label] = Clock::now();
}

void Profiler::stop(const string& label){
    if (timePoints.count(label)){
        auto now = Clock::now();
        auto duration = std::chrono::duration<double, std::milli>(now - timePoints[label]);
        if (!clocks.count(label)) clocks[label] = {0.0, 0};
        clocks[label] = clocks[label] + make_pair(duration.count(), 1);
    }
}

void Profiler::add(const Profiler& pro){
    for (const auto& cl : pro.clocks){
        auto label = cl.first;
        if (!clocks.count(label)) clocks[label] = {0.0, 0};
        clocks[label] = clocks[label] + cl.second;
    }
}

void Profiler::print() const{
    for (const auto& cl : clocks){
        auto label = cl.first;
        if (!label.size()) label = "Ø";
        fmt::println("{}[count={} avg={}ms]", label, cl.second.second, cl.second.first/cl.second.second);
    }
}