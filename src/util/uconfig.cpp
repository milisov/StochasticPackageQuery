#include <boost/property_tree/ini_parser.hpp>

#include "uconfig.hpp"

using std::make_shared;

fs::path getProjectDir(){
	return fs::current_path().parent_path();
}

Config::Config(){
	fs::path configPath = getProjectDir() / "config.cfg";
	boost::property_tree::ini_parser::read_ini(configPath.string(), pt);
}

shared_ptr<Config> Config::config = nullptr;

shared_ptr<Config> Config::getInstance(){
    if (config == nullptr){
        config = std::make_shared<Config>();
    }
    return config;
}