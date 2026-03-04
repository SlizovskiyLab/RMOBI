// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

#ifndef PARSER_H
#define PARSER_H

#include <string>
#include <filesystem>
#include <map> // Required for std::map
#include "graph.h"

namespace std {
    template <>
    struct hash<std::pair<int, bool>> {
        std::size_t operator()(const std::pair<int, bool>& p) const {
            return std::hash<int>()(p.first) ^ (std::hash<bool>()(p.second) << 1);
        }
    };
}


// The signature of parseData is updated to include a map for patient-disease associations.
void parseData(const std::filesystem::path& filename, Graph& graph, std::map<int, std::string>& patientToDiseaseMap, bool includeSNPConfirmationARGs, bool excludeMetals);
void addEdge(Graph& graph, const Node& src, const Node& tgt, bool isColo, int patientID = -1);
void addTemporalEdges(Graph& graph);

#endif // PARSER_H