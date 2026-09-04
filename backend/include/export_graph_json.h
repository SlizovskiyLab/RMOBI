// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

#pragma once
#include <string>
#include <map>
#include "graph.h"

/**
 * @file ExportGraphJson.h
 * @brief Contains functions for exporting graph data to JSON format for the web app.
 *
 * This file contains functions for exporting a graph to JSON format, which is used by the web app for visualization and analysis.
 */

bool exportGraphToJsonSimple(const Graph& g, const std::string& outPathStr, const std::map<int, std::string>& patientToDiseaseMap, const std::map<int, std::string>& patientToStudyMap);

bool exportParentGraphToJson(const Graph& g, const std::string& outPathStr, const std::map<int, std::string>& patientToDiseaseMap, const std::map<int, std::string>& patientToStudyMap, bool showLabels = true);

void exportColocalizationsToJSONByDisease(
    const std::map<std::tuple<int,int,int>, std::set<Timepoint>>& colocalizationByIndividual,
    const std::map<int, std::string>& patientToDiseaseMap,
    const std::map<int, std::string>& patientToStudyMap,
    const std::string& jsonOutputPath  // path to the final JSON file
);
