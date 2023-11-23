#ifndef UCONFIG_HPP
#define UCONFIG_HPP

#include <boost/property_tree/ptree.hpp>
#include <memory>
#include <filesystem>

namespace fs = std::filesystem;

using std::shared_ptr;

fs::path getProjectDir();

class Config{
protected:
    static shared_ptr<Config> config;
public:
    boost::property_tree::ptree pt;
    unsigned int nPhysicalCores, nLogicalCores;
public:
    Config();
    Config(Config &other) = delete;
    void operator=(const Config&) = delete;
    static shared_ptr<Config> getInstance();
};

#endif