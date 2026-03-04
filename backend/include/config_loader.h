// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

// config_loader.h

#ifndef CONFIG_LOADER_H
#define CONFIG_LOADER_H

#include <string>
#include "../external/json.hpp"

struct Config {
    std::string input_data_path;

    std::string output_base;
    std::string output_disease;
    std::string output_mge_group;
    std::string output_top_entities;
    std::string output_temporal_dynamics;

    std::string output_emerge;
    std::string output_disappear;
    std::string output_transfer;
    std::string output_persist;

    std::string output_top_arg;
    std::string output_top_mge;
    std::string output_top_colocalizations;

    std::string viz_interaction;
    std::string viz_parent;
    std::string viz_temporal_dynamics;
};


Config loadConfig(const std::string& filename);
void createOutputDirectories(const Config& cfg);

#endif // CONFIG_LOADER_H