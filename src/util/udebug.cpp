#include "udebug.hpp"

#include <omp.h>
#include <fmt/core.h>

Profiler::Profiler(){
}

void Profiler::clock(string label){
    timePoints[label] = Clock::now();
}

void Profiler::stop(string label){
    if (timePoints.count(label)){
        auto now = Clock::now();
        auto duration = std::chrono::duration<double, std::milli>(now - timePoints[label]);
        if (!clocks.count(label)) clocks[label] = {0.0, 0};
        clocks[label] = clocks[label] + make_pair(duration.count(), 1);
    }
}

void Profiler::add(Profiler pro){
    for (auto cl : pro.clocks){
        auto label = cl.first;
        if (!clocks.count(label)) clocks[label] = {0.0, 0};
        clocks[label] = clocks[label] + cl.second;
    }
}

void Profiler::print(){
    for (auto cl : clocks){
        auto label = cl.first;
        if (!label.size()) label = "Ø";
        fmt::println("{}[count={} avg={}ms]", label, cl.second.second, cl.second.first/cl.second.second);
    }
}