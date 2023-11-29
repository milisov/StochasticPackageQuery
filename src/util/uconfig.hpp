#ifndef UCONFIG_HPP
#define UCONFIG_HPP

#include <boost/property_tree/ptree.hpp>
#include <memory>
#include <random>
#include <filesystem>

namespace fs = std::filesystem;

using std::shared_ptr;

fs::path getProjectDir();

class Config{
private:
    long long seedMode;
    std::random_device rd;
protected:
    static shared_ptr<Config> config;
public:
    boost::property_tree::ptree pt;
    unsigned int nPhysicalCores, nLogicalCores;
public:
    Config();
    Config(Config &other) = delete;
    void operator=(const Config&) = delete;
    unsigned int seed();
    static shared_ptr<Config> getInstance();
};

#endif