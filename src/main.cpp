// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

#include <iostream>
#include <filesystem>
#include "../include/Timepoint.h"
#include "../include/graph.h"
#include "../include/parser.h"
#include "../include/id_maps.h"
#include "../include/traversal.h"
#include "../include/analysis.h"
#include "../include/export.h"
#include "../include/graph_utils.h"
#include "../include/export_graph_json.h" 
#include "../include/config_loader.h"  

/* Main entry point: parse arguments, load data, call functions */

namespace fs = std::filesystem;

fs::path data_file;
fs::path interaction_json_path;
fs::path parent_json_path;
fs::path temporal_dynamics_json_path;
fs::path top_entities_output_dir;
fs::path top_colocalizations_output;
fs::path disease_type_output;
fs::path mge_group_output;



int main() {
    try {
        Config cfg = loadConfig("config/paths.json");
        data_file = fs::path(cfg.input_data_path);
        interaction_json_path = fs::path(cfg.viz_interaction);
        parent_json_path = fs::path(cfg.viz_parent);
        temporal_dynamics_json_path = fs::path(cfg.viz_temporal_dynamics);
        top_entities_output_dir = fs::path(cfg.output_base);
        top_colocalizations_output = fs::path(cfg.output_top_colocalizations);


    } catch (const std::exception& e) {
        std::cerr << "Config error: " << e.what() << "\n";
        return 1;
    }
    Graph g;
    std::map<int, std::string> patientToDiseaseMap;
    std::unordered_map<Node, std::unordered_set<Node>> adjacency;

    // parse the data file and construct the graph (true to exclude ARGs requiring SNP confirmation, true to exclude metals)
    parseData(data_file, g, patientToDiseaseMap, true, false);

    addTemporalEdges(g);  
    buildAdjacency(g, adjacency);

    /******************************** Graph Statistics  ************************************/
    writeGraphStatisticsCSV(g, adjacency, "viz/output/graph_statistics.csv");

    /******************************** Traversal of Graph  ************************************/
    std::map<std::pair<int, int>, std::multiset<Timepoint>> colocalizationTimeline;
    traverseAdjacency(g, adjacency, colocalizationTimeline);

    /******************************** Traversal of Graph  ************************************/
    std::map<std::tuple<int, int, int>, std::set<Timepoint>> colocalizationByIndividual;
    traverseGraph(g, colocalizationByIndividual);
    std::map<std::pair<int, int>, std::set<int>> globalPairToPatients;
    

    /********************************* Colocalizations by Timepoints ************************************/
    writeAllDiseasesTemporalDynamicsCounts(colocalizationByIndividual, patientToDiseaseMap);
    writeTemporalDynamicsCountsForMGEGroup(colocalizationByIndividual);
    mostProminentEntities(g);
    getTopARGMGEPairsByFrequencyWODonor(colocalizationByIndividual, 10, patientToDiseaseMap, top_colocalizations_output.string());
    

    exportColocalizations(g, colocalizationByIndividual);

    // /************************************* Graph Visualization ***********************************/

    Graph amrGraphNet = g;
    exportGraphToJsonSimple(amrGraphNet, interaction_json_path.string(), patientToDiseaseMap);
    exportParentGraphToJson(amrGraphNet, parent_json_path.string(), patientToDiseaseMap, true);
    exportColocalizationsToJSONByDisease(colocalizationByIndividual, patientToDiseaseMap, temporal_dynamics_json_path.string());


    return 0;

}






