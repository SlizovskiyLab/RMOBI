// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

#include <fstream>
#include <iostream>
#include "../include/config_loader.h"

Config loadConfig(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open config file: " + filename);
    }

    nlohmann::json j;
    file >> j;

    Config cfg;
    cfg.input_data_path = j.at("input").at("input_data").get<std::string>();

    auto output = j.at("output");
    cfg.output_base       = output.at("base").get<std::string>();
    cfg.output_disease    = output.at("disease_type").get<std::string>();
    cfg.output_mge_group  = output.at("mge_group").get<std::string>();

    auto td = output.at("temporal_dynamics");
    cfg.output_emerge     = td.at("emerge").get<std::string>();
    cfg.output_disappear  = td.at("disappear").get<std::string>();
    cfg.output_transfer   = td.at("transfer").get<std::string>();
    cfg.output_persist    = td.at("persist").get<std::string>();

    auto te = output.at("top_entities");
    cfg.output_top_arg              = te.at("top_arg").get<std::string>();
    cfg.output_top_mge              = te.at("top_mge").get<std::string>();
    cfg.output_top_colocalizations  = te.at("top_colocalizations").get<std::string>();

    cfg.viz_interaction        = j.at("viz").at("interaction_json").get<std::string>();
    cfg.viz_parent             = j.at("viz").at("parent_json").get<std::string>();
    cfg.viz_temporal_dynamics  = j.at("viz").at("temporal_dynamics_disease").get<std::string>();

    return cfg;
}

void createOutputDirectories(const Config& cfg) {
    using namespace std::filesystem;

    create_directories(cfg.output_base);
    create_directories(cfg.output_disease);
    create_directories(cfg.output_mge_group);
        // temporal dynamics parents
    create_directories(path(cfg.output_emerge).parent_path());
    create_directories(path(cfg.output_disappear).parent_path());
    create_directories(path(cfg.output_transfer).parent_path());
    create_directories(path(cfg.output_persist).parent_path());

    // top entities parents
    create_directories(path(cfg.output_top_arg).parent_path());
    create_directories(path(cfg.output_top_mge).parent_path());
    create_directories(path(cfg.output_top_colocalizations).parent_path());

}

