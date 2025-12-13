// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

#include "../include/parser.h"
#include "../include/id_maps.h"
#include <fstream>
#include <sstream>
#include <filesystem>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <set>
#include <algorithm>
#include <map>
#include <cctype>

/* Read input files (CSV), extract (individual, ARG, MGE, timepoint) data */

void addEdge(Graph& graph, const Node& src, const Node& tgt, bool isColo, int patientID) {
    if (isColo) {
        Node s = std::min(src, tgt);
        Node t = std::max(src, tgt);
        Edge search_edge = {s, t, true, {}, 0};

        auto it = graph.edges.find(search_edge);
        if (it != graph.edges.end()) {
            Edge modified_edge = *it;
            graph.edges.erase(it);
            if (patientID != -1) modified_edge.individuals.insert(patientID);
            graph.edges.insert(modified_edge);
        } else {
            if (patientID != -1) search_edge.individuals.insert(patientID);
            graph.edges.insert(search_edge);
        }
    } else {
        Edge search_edge = {src, tgt, false, {}, 0};
        auto it = graph.edges.find(search_edge);
        if (it != graph.edges.end()) {
            Edge modified_edge = *it;
            graph.edges.erase(it);
            modified_edge.weight++;
            graph.edges.insert(modified_edge);
        } else {
            search_edge.weight = 1;
            graph.edges.insert(search_edge);
        }
    }
}

// This function reads a CSV file containing patient data and constructs a graph.
// It extracts ARG and MGE labels, maps them to IDs, and creates nodes and edges
void parseData(const std::filesystem::path& filename, Graph& graph, std::map<int, std::string>& patientToDiseaseMap, bool includeSNPConfirmationARGs, bool excludeMetals) {
    std::ifstream infile(filename);
    std::string line;
    std::vector<std::string> headers;
    bool isHeader = true;

    std::unordered_map<std::string, Timepoint> columnToTimepoint;

    while (std::getline(infile, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> tokens;

        while (std::getline(ss, token, ',')) {
            tokens.push_back(token);
        }

        if (isHeader) {
            headers = tokens;
            for (const std::string& col_const : headers) {
                std::string col = col_const;
                 while (!col.empty() && isspace(col.back())) {
                    col.pop_back();
                }

                if (col == "Donor") columnToTimepoint[col] = Timepoint::Donor;
                else if (col == "PreFMT") columnToTimepoint[col] = Timepoint::PreFMT;
                else if (col.rfind("PostFMT_", 0) == 0) {
                    std::string day_str = col.substr(8);
                    try {
                        columnToTimepoint[col] = static_cast<Timepoint>(std::stoi(day_str));
                    } catch (const std::invalid_argument& e) {
                        std::cerr << "Warning: Could not parse day from header column: " << col << std::endl;
                    }
                }
            }
            isHeader = false;
            continue;
        }

        if (tokens.size() < 4) continue;

        int patientID = std::stoi(tokens[0]);
        std::string diseaseType = tokens[1];
        patientToDiseaseMap[patientID] = diseaseType; 

        std::string argLabel = tokens[2];
        std::string mgeLabel = tokens[3];

        int argID = getARGId(argLabel);
        int mgeID = getMGEId(mgeLabel);

        if (argID == -1 || mgeID == -1) continue;

        for (size_t i = 4; i < tokens.size(); ++i) {
            std::string colName = headers[i];
            if (columnToTimepoint.count(colName) && (tokens[i] == "1" || tokens[i] == "2")) {
                Timepoint tp = columnToTimepoint.at(colName);
                
                bool requiresSNPConfirmation = argIDSNPConfirmation.count(argID) ? argIDSNPConfirmation.at(argID) : false;
                if (includeSNPConfirmationARGs && requiresSNPConfirmation) continue;
                
                std::string argResistance = argResistanceMap.count(argID) ? argResistanceMap.at(argID) : "Unknown";
                if (excludeMetals && argResistance != "Drugs") continue;

                Node argNode = {argID, true, tp, requiresSNPConfirmation};
                Node mgeNode = {mgeID, false, tp, false};

                graph.nodes.insert(argNode);
                graph.nodes.insert(mgeNode);

                addEdge(graph, argNode, mgeNode, true, patientID);
            }
        }
    }
}



/**
 * This function adds patient-specific temporal edges between nodes.
 * It creates directed edges ONLY between chronologically adjacent timepoints for the same gene within the same patient.
 * @param graph The graph to which temporal edges will be added.
 */
void addTemporalEdges(Graph& graph) {
    std::map<int, std::set<Node>> nodesByPatient;
    for (const auto& edge : graph.edges) {
        if (!edge.isColo) continue;
        for (int patientID : edge.individuals) {
            nodesByPatient[patientID].insert(edge.source);
            nodesByPatient[patientID].insert(edge.target);
        }
    }

    for (const auto& [patientID, nodeSet] : nodesByPatient) {
        std::unordered_map<std::pair<int, bool>, std::vector<Node>> groupedNodes;
        for (const Node& node : nodeSet) {
            groupedNodes[{node.id, node.isARG}].push_back(node);
        }

        for (auto const& [key, nodeGroup] : groupedNodes) {
            auto nodes = nodeGroup;
            std::sort(nodes.begin(), nodes.end());
            for (size_t i = 0; i < nodes.size() - 1; ++i) {
                const Node& sourceNode = nodes[i];
                const Node& targetNode = nodes[i + 1];
                addEdge(graph, sourceNode, targetNode, false, -1);
            }
        }
    }
}