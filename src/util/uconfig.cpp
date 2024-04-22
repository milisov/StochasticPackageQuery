#include <boost/property_tree/ini_parser.hpp>
#include <thread>

#include "uconfig.hpp"

using std::make_shared;

fs::path getProjectDir(){
	return fs::current_path().parent_path().parent_path().parent_path();
}

Config::Config(){
	fs::path configPath = getProjectDir() / "_config.cfg";
	boost::property_tree::ini_parser::read_ini(configPath.string(), pt);
    auto logicalCores = pt.get_optional<unsigned int>("hardware.logical_cores");
    if (logicalCores) nLogicalCores = *logicalCores;
    else nLogicalCores = std::thread::hardware_concurrency();

    auto physicalCores = pt.get_optional<unsigned int>("hardware.physical_cores");
    if (physicalCores) nPhysicalCores = *physicalCores;
    else nPhysicalCores = std::thread::hardware_concurrency() / 2;

    seedMode = pt.get<long long>("parameters.seed");
    if (seedMode == -1) seedMode = rd();
}

shared_ptr<Config> Config::config = nullptr;

shared_ptr<Config> Config::getInstance(){
    if (config == nullptr){
        config = std::make_shared<Config>();
    }
    return config;
}

unsigned int Config::seed(){
    if (seedMode >= 0) return static_cast<unsigned int>(seedMode);
    return rd();
}