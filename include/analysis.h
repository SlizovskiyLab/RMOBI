// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

#ifndef QUERY_ENGINE_H
#define QUERY_ENGINE_H
#include <map>
#include <set>
#include <tuple>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include "graph.h"
#include "Timepoint.h"

bool isPostFMT(const Timepoint& tp);
bool isPreFMT(const Timepoint& tp);
bool isDonor(const Timepoint& tp);

void getPatientwiseColocalizationsByCriteria(
    const Graph& graph,
    const std::map<std::tuple<int, int, int>, std::set<Timepoint>>& colocalizationByIndividual,
    bool donorStatus,
    bool preFMTStatus,
    bool postFMTStatus,
    const std::string& label,
    const std::string& csvFile = "",
    bool append = false
);

void getTopARGMGEPairsByFrequency(const std::map<std::tuple<int, int, int>, std::set<Timepoint>>& colocalizations, int topN = -1);
void getTopARGMGEPairsByFrequencyWODonor(const std::map<std::tuple<int, int, int>, std::set<Timepoint>>& colocalizations, int topN = -1, 
    const std::map<int, std::string>& patientToDiseaseMap = {}, const std::string& top_colocalizations_output = {});
void getConnectedMGE(const Graph& graph, const std::set<Edge>& edges, int mgeId, const std::string& mgeName);

void getColocalizationsByCriteria(
    const std::map<std::tuple<int, int, int>, std::set<Timepoint>>& colocalizationByIndividual,
    bool donorStatus,
    bool preFMTStatus,
    bool postFMTStatus,
    std::map<std::pair<int, int>, std::set<int>>& globalPairToPatients
);
void getTopARGMGEPairsByUniquePatients(
    const std::map<std::pair<int, int>, std::set<int>>& globalPairToPatients,
    int topN = -1,
    const std::string& label = "All Patients"
);

void writeTemporalDynamicsCountsForDisease(
    const std::string& disease,
    std::map<std::tuple<int, int, int>, std::set<Timepoint>>& colocalizationByIndividual,
    const std::map<int, std::string>& patientToDiseaseMap
);

void writeAllDiseasesTemporalDynamicsCounts(
    std::map<std::tuple<int,int,int>, std::set<Timepoint>>& colocalizationByIndividual,
    const std::map<int, std::string>& patientToDiseaseMap
);

void writeTemporalDynamicsCountsForMGEGroup(
    const std::map<std::tuple<int,int,int>, std::set<Timepoint>>& colocalizationByIndividual
);

void writeColocalizationsToCSV(
    const std::map<std::tuple<int, int, int>, std::set<Timepoint>>& colocs,
    const std::string& filename,
    const std::string& label,
    bool append = false
);

void exportColocalizations(const Graph& g,
    const std::map<std::tuple<int,int,int>, std::set<Timepoint>>& colocalizationByIndividual
);


void writeGraphStatisticsCSV(
    const Graph& g,
    const std::unordered_map<Node, std::unordered_set<Node>>& adjacency,
    const std::string& filename
);

void analyzeColocalizations(const Graph& g, 
                             const std::unordered_map<Node, std::unordered_set<Node>>& adjacency);

void analyzeColocalizationsCollectively(const Graph& g, 
                                          const std::unordered_map<Node, std::unordered_set<Node>>& adjacency);

void mostProminentEntities(const Graph& g);

void writeTopEntitiesToCSV(const std::string& filename, const std::vector<std::pair<int, int>>& entities, bool isARG);

void writeCSV(
    const std::string& filename,
    const std::vector<std::string>& header,
    const std::vector<std::vector<std::string>>& rows, bool append = false
);

#endif