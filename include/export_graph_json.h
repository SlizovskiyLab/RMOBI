// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

#pragma once
#include <string>
#include <map>
#include "graph.h"

bool exportGraphToJsonSimple(const Graph& g, const std::string& outPathStr, const std::map<int, std::string>& patientToDiseaseMap);

bool exportParentGraphToJson(const Graph& g, const std::string& outPathStr, const std::map<int, std::string>& patientToDiseaseMap, bool showLabels = true);

void exportColocalizationsToJSONByDisease(
    const std::map<std::tuple<int,int,int>, std::set<Timepoint>>& colocalizationByIndividual,
    const std::map<int, std::string>& patientToDiseaseMap,
    const std::string& jsonOutputPath  // path to the final JSON file
);
