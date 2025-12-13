// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

#ifndef GRAPH_UTILS_H
#define GRAPH_UTILS_H

#include "graph.h"
#include <string>
#include <unordered_set>
#include <map>

Graph filterGraphByARGName(const Graph& g, const std::string& argName);
Graph filterGraphByMGEName(const Graph& g, const std::string& mgeName);
Graph filterGraphByMGEGroup(const Graph& g, const std::string& groupName);
Graph filterGraphByTimepoint(const Graph& g, const std::string& timepointCategory);
Graph filterGraphByARGAndMGENames(const Graph& g, const std::string& argName, const std::string& mgeName);
Graph filterGraphByDisease(const Graph& g, const std::string& disease, const std::map<int, std::string>& patientToDiseaseMap);

#endif // GRAPH_UTILS_H