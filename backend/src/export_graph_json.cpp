// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

// export_graph_json.cpp
#include <fstream>
#include <string>
#include <filesystem>
#include <algorithm>
#include <map>
#include <set> 
#include "../include/graph.h"
#include "../include/id_maps.h"
#include "../include/export.h"
#include "../external/json.hpp"
#include "../include/graph.h"
#include "../include/Timepoint.h"
#include "../include/analysis.h"
#include "../include/parser.h" 

using nlohmann::json;
namespace fs = std::filesystem;


// static std::string getTimepointColor(const Timepoint& tp) {
//     int timeValue = static_cast<int>(tp);
//     if (timeValue == 1000) return "yellow";
//     if (timeValue == 0)    return "red";
//     if (timeValue > 0 && timeValue < 31) return "#99D2FF";
//     if (timeValue > 30 && timeValue < 61) return "#4D9DFF";
//     if (timeValue > 60)    return "#3A6EFF";
//     return "green"; // fallback
// }

// static std::string getTimepointCategory(const Timepoint& tp) {
//     int timeValue = static_cast<int>(tp);
//     if (timeValue == 1000) return "donor";
//     if (timeValue == 0)    return "pre";
//     if (timeValue > 0 && timeValue < 31) return "post1";
//     if (timeValue > 30 && timeValue < 61) return "post2";
//     if (timeValue > 60)    return "post3";
//     return "unknown"; // Fallback for any unexpected values
// }


// auto timepointOrder = [](Timepoint tp) -> int {
//     if (tp == Timepoint::Donor) return -1;
//     if (tp == Timepoint::PreFMT) return 0;
//     return static_cast<int>(tp);
// };

//ordering/adjacency decisions (Donor < PreFMT < PostFMT)
// static inline int tpSort(Timepoint tp) {
//     return timepointOrder(tp);
// }

static inline int tpJson(Timepoint tp) {
    const int v = static_cast<int>(tp);
    return v;
}


std::string getLabel(const Node& node) {
    std::string label = node.isARG ? getARGName(node.id) : getMGENameForLabel(node.id);
    return label;
}


bool exportGraphToJsonSimple(const Graph& g,
                             const std::string& outPathStr,
                             const std::map<int, std::string>& patientToDiseaseMap)
{
    json j;
    j["nodes"] = json::array();
    j["links"] = json::array();

    std::unordered_set<Node> active_nodes;
    std::set<std::pair<Node, Node>> processedColoEdges;

    for (const Edge& edge : g.edges) {
        if (edge.source == edge.target) continue;

        active_nodes.insert(edge.source);
        active_nodes.insert(edge.target);

        // For colocalization, deduplicate undirected edges
        if (edge.isColo) {
            auto canon = std::minmax(edge.source, edge.target);
            if (processedColoEdges.count(canon)) continue;
            processedColoEdges.insert(canon);
        }

        json diseases = json::array();
        json diseaseCounts = json::object();

        std::map<std::string, int> diseaseCounter;

        for (int patientID : edge.individuals) {
            auto it = patientToDiseaseMap.find(patientID);
            if (it != patientToDiseaseMap.end()) {
                diseaseCounter[it->second]++;
            }
        }

        for (const auto& [diseaseName, count] : diseaseCounter) {
            diseases.push_back(diseaseName);
            diseaseCounts[diseaseName] = count;
        }

        j["links"].push_back({
            {"source", getNodeName(edge.source)},
            {"target", getNodeName(edge.target)},

            // counts/filters
            {"individualCount", static_cast<int>(edge.individuals.size())},
            {"isColo", edge.isColo},
            {"diseases", diseases},
            {"diseaseCounts", diseaseCounts},

            // timepoints for fast temporal styling in JS
            {"sourceTimepoint", static_cast<int>(edge.source.timepoint)},
            {"targetTimepoint", static_cast<int>(edge.target.timepoint)}
        });
    }

    for (const Node& n : active_nodes) {
        std::string mgeGroup = "";
        if (!n.isARG) {
            mgeGroup = getMGEGroupName(n.id);
        }

        j["nodes"].push_back({
            {"id",        getNodeName(n)},
            {"label",     getLabel(n)},
            {"isARG",     n.isARG},
            {"timepoint", static_cast<int>(n.timepoint)},
            {"mgeGroup",  mgeGroup},
            {"argResistance", n.isARG ? getARGResistance(n.id) : ""},
            {"argResistanceGroup", n.isARG ? getARGGroupName(n.id) : ""}
        });
    }

    // Write to disk
    std::ofstream out(outPathStr);
    if (!out) {
        std::cerr << "[exportGraphToJsonSimple] Cannot open " << outPathStr << " for write\n";
        return false;
    }
    out << j.dump(2) << '\n';

    std::cerr << "[exportGraphToJsonSimple] Wrote nodes=" << j["nodes"].size()
              << " links=" << j["links"].size()
              << " to " << outPathStr << "\n";
    return true;
}
 


bool exportParentGraphToJson(const Graph& g,
                             const std::string& outPathStr,
                             const std::map<int, std::string>& patientToDiseaseMap,
                             bool showLabels)
{
    json j;
    j["nodes"] = json::array();
    j["links"] = json::array();

    auto timepointOrder = [](Timepoint tp) -> int {
        if (tp == Timepoint::Donor) return -1;
        if (tp == Timepoint::PreFMT) return 0;
        return static_cast<int>(tp);
    };

    // (argId, mgeId, tp) -> Parent_X
    int counter = 0;
    std::map<std::tuple<int,int,Timepoint>, std::string> uniqueParents;

    // Aggregation for node payload
    std::map<std::tuple<int,int,Timepoint>, std::string> labelByKey;
    std::map<std::tuple<int,int,Timepoint>, std::string> mgeGroupByKey;
    std::map<std::tuple<int,int,Timepoint>, std::string> argResistanceByKey;
    std::map<std::tuple<int,int,Timepoint>, std::string> argGroupByKey;
    std::map<std::tuple<int,int,Timepoint>, std::map<std::string, std::set<int>>> diseaseToPatientsByKey;
    std::map< std::tuple<int,int,Timepoint,Timepoint>, std::map<std::string, std::set<int>>> temporalDiseaseToPatients;

    // patient -> (arg,mge) -> set<tp>
    std::map<int, std::map<std::pair<int,int>, std::set<Timepoint>>> nodesByPatient;

    // 0) Build nodesByPatient from ALL colocalization edges (same idea as addTemporalEdges)
    for (const Edge& edge : g.edges) {
        if (!edge.isColo) continue;

        // Determine which endpoint is ARG/MGE
        const Node& argNode = edge.source.isARG ? edge.source : edge.target;
        const Node& mgeNode = edge.source.isARG ? edge.target : edge.source;

        for (int patientID : edge.individuals) {
            // Only record the ARG–MGE pair at THIS timepoint for this patient.
            nodesByPatient[patientID][{argNode.id, mgeNode.id}].insert(argNode.timepoint);
        }
    }

    // 1) Create parent nodes for every (arg,mge,tp) implied by patient presence
    //    This guarantees donor nodes exist too, if donor is present for any patient.
    for (const auto& [patientID, pairMap] : nodesByPatient) {
        for (const auto& [pairKey, tpSet] : pairMap) {
            const int argId = pairKey.first;
            const int mgeId = pairKey.second;

            for (Timepoint tp : tpSet) {
                auto key = std::make_tuple(argId, mgeId, tp);
                if (!uniqueParents.count(key)) {
                    uniqueParents[key] = "Parent_" + std::to_string(++counter);
                    mgeGroupByKey[key] = getMGEGroupName(mgeId);
                    argResistanceByKey[key] = getARGResistance(argId);
                    argGroupByKey[key] = getARGGroupName(argId);

                    labelByKey[key] = showLabels
                        ? (getARGName(argId) + "+" + getMGENameForLabel(mgeId))
                        : "";
                }

                // diseaseCounts: count unique patients per disease at this (arg,mge,tp)
                auto itDis = patientToDiseaseMap.find(patientID);
                if (itDis != patientToDiseaseMap.end()) {
                    diseaseToPatientsByKey[key][itDis->second].insert(patientID);
                }
            }
        }
    }

    // 2) Emit nodes (same JSON structure)
    for (const auto& [key, parentName] : uniqueParents) {
        const int argId = std::get<0>(key);
        const int mgeId = std::get<1>(key);
        const Timepoint tp = std::get<2>(key);

        json diseaseCounts = json::object();
        auto itCounts = diseaseToPatientsByKey.find(key);
        if (itCounts != diseaseToPatientsByKey.end()) {
            for (const auto& [diseaseName, patients] : itCounts->second) {
                diseaseCounts[diseaseName] = static_cast<int>(patients.size());
            }
        }

        j["nodes"].push_back({
            {"id", parentName},
            {"label", labelByKey[key]},
            {"argId", argId},
            {"mgeId", mgeId},
            {"timepoint", tpJson(tp)},   // numeric, safe
            {"mgeGroup", mgeGroupByKey[key]},
            {"argResistance", argResistanceByKey[key]},
            {"argGroup", argGroupByKey[key]},
            {"diseaseCounts", diseaseCounts}
        });
    }

        // ------------------------------------------------------------------
        // 3) Build temporal-edge disease support patient-by-patient
        // For each patient and each ARG-MGE pair, create only adjacent temporal transitions:
        //   - Donor -> Pre, - Donor -> firstPost (if Pre missing),   - Pre -> firstPost, //   - post[i] -> post[i+1]
        // ------------------------------------------------------------------
        auto hasTp = [](const std::set<Timepoint>& s, Timepoint tp) {
            return s.find(tp) != s.end();
        };

        for (const auto& [patientID, pairMap] : nodesByPatient) {
            auto itDis = patientToDiseaseMap.find(patientID);
            const bool hasDisease = (itDis != patientToDiseaseMap.end());
            const std::string diseaseName = hasDisease ? itDis->second : "";

            for (const auto& [pairKey, tpSet] : pairMap) {
                const int argId = pairKey.first;
                const int mgeId = pairKey.second;

                const bool hasDonor = hasTp(tpSet, Timepoint::Donor);
                const bool hasPre   = hasTp(tpSet, Timepoint::PreFMT);

                std::vector<Timepoint> postTps;
                for (Timepoint tp : tpSet) {
                    if (tp != Timepoint::Donor && tp != Timepoint::PreFMT) {
                        postTps.push_back(tp);
                    }
                }

                std::sort(postTps.begin(), postTps.end(),
                        [&](Timepoint a, Timepoint b) {
                            return timepointOrder(a) < timepointOrder(b);
                        });

                auto addTemporalSupport = [&](Timepoint srcTp, Timepoint tgtTp) {
                    if (!hasDisease) return;
                    temporalDiseaseToPatients[{argId, mgeId, srcTp, tgtTp}][diseaseName].insert(patientID);
                };

                // Donor -> Pre OR Donor -> firstPost (if Pre missing)
                if (hasDonor) {
                    if (hasPre) {
                        addTemporalSupport(Timepoint::Donor, Timepoint::PreFMT);
                    } else if (!postTps.empty()) {
                        addTemporalSupport(Timepoint::Donor, postTps.front());
                    }
                }

                // Pre -> firstPost
                if (hasPre && !postTps.empty()) {
                    addTemporalSupport(Timepoint::PreFMT, postTps.front());
                }

                // Chain post -> post
                for (size_t i = 0; i + 1 < postTps.size(); ++i) {
                    addTemporalSupport(postTps[i], postTps[i + 1]);
                }
            }
        }

        // ------------------------------------------------------------------
        // 4) Emit deduplicated temporal links
        // ------------------------------------------------------------------
        std::set<std::tuple<std::string,std::string,int,int>> emittedTemporal;

        for (const auto& [tempKey, diseaseMap] : temporalDiseaseToPatients) {
            const int argId      = std::get<0>(tempKey);
            const int mgeId      = std::get<1>(tempKey);
            const Timepoint srcTp = std::get<2>(tempKey);
            const Timepoint tgtTp = std::get<3>(tempKey);

            auto srcKey = std::make_tuple(argId, mgeId, srcTp);
            auto tgtKey = std::make_tuple(argId, mgeId, tgtTp);

            auto itS = uniqueParents.find(srcKey);
            auto itT = uniqueParents.find(tgtKey);
            if (itS == uniqueParents.end() || itT == uniqueParents.end()) continue;

            const std::string& sName = itS->second;
            const std::string& tName = itT->second;

            auto dk = std::make_tuple(sName, tName, timepointOrder(srcTp), timepointOrder(tgtTp));
            if (emittedTemporal.count(dk)) continue;
            emittedTemporal.insert(dk);

            json diseases = json::array();
            json diseaseCounts = json::object();
            int individualCount = 0;

            for (const auto& [diseaseName, patients] : diseaseMap) {
                diseases.push_back(diseaseName);
                diseaseCounts[diseaseName] = static_cast<int>(patients.size());
                individualCount += static_cast<int>(patients.size());
            }

            j["links"].push_back({
                {"source", sName},
                {"target", tName},
                {"isColo", false},
                {"individualCount", individualCount},
                {"diseases", diseases},
                {"diseaseCounts", diseaseCounts},
                {"sourceTimepoint", tpJson(srcTp)},
                {"targetTimepoint", tpJson(tgtTp)}
            });
        }

        // ------------------------------------------------------------------
        // 5) Write JSON to disk
        // ------------------------------------------------------------------
        std::ofstream out(outPathStr);
        if (!out) {
            std::cerr << "[exportParentGraphToJson] Cannot open " << outPathStr << " for write\n";
            return false;
        }
        out << j.dump(2) << '\n';

        std::cerr << "[exportParentGraphToJson] Wrote parent-nodes=" << j["nodes"].size()
                << " links=" << j["links"].size()
                << " to " << outPathStr << "\n";
        return true;
    }



// void exportColocalizationsToJSONByDisease(
//     const Graph& g,
//     const std::map<std::tuple<int,int,int>, std::set<Timepoint>>& colocalizationByIndividual,
//     const std::map<int, std::string>& patientToDiseaseMap,
//     const std::string& jsonOutputDir
// ) {
//     // Create map: disease -> (colocalization pair -> status counts)
//     std::map<std::string, std::map<std::string, std::map<std::string, int>>> diseaseColocCounts;

//     for (const auto& [tuple, tps] : colocalizationByIndividual) {
//         auto [arg, mge, patient] = tuple;

//         auto diseaseIt = patientToDiseaseMap.find(patient);
//         if (diseaseIt == patientToDiseaseMap.end()) continue;
//         const std::string& disease = diseaseIt->second;

//         std::string pairName = getARGName(arg) + "–" + getMGEName(mge);

//         bool hasDonor = std::any_of(tps.begin(), tps.end(), [](const Timepoint& tp) { return isDonor(tp); });
//         bool hasPre = std::any_of(tps.begin(), tps.end(), [](const Timepoint& tp) { return isPreFMT(tp); });
//         bool hasPost = std::any_of(tps.begin(), tps.end(), [](const Timepoint& tp) { return isPostFMT(tp); });


//         std::string status;
//         if (hasPost && !hasPre && !hasDonor) status = "emerged";
//         else if (hasPre && !hasPost && !hasDonor) status = "disappeared";
//         else if (hasDonor && hasPost && !hasPre) status = "transferred";
//         else if (hasPre && hasPost) status = "persisted";
//         else continue; // skip other patterns

//         diseaseColocCounts[disease][pairName][status]++;
//     }

//     // Build JSON structure
//     json rootJson = json::array();

//     for (const auto& [disease, colocMap] : diseaseColocCounts) {
//         json diseaseEntry;
//         diseaseEntry["disease"] = disease;
//         diseaseEntry["data"] = json::array();

//         for (const auto& [pair, statusMap] : colocMap) {
//             for (const auto& [status, count] : statusMap) {
//                 json entry;
//                 entry["colocalization"] = pair;
//                 entry["status"] = status;
//                 entry["patients"] = count;
//                 diseaseEntry["data"].push_back(entry);
//             }
//         }

//         rootJson.push_back(diseaseEntry);

//         // // Write per-disease JSON
//         // std::filesystem::path outFile = std::filesystem::path(jsonOutputDir) / (disease + "_colocalizations.json");
//         // std::ofstream ofs(outFile);
//         // ofs << std::setw(2) << diseaseEntry << std::endl;
//     }

//     std::filesystem::path outFile(jsonOutputDir);
//     std::ofstream all(outFile.string()); 
//     if (!all.is_open()) {
//         throw std::runtime_error("Failed to open output JSON file: " + outFile.string());
//     }
//     all << std::setw(2) << rootJson << std::endl;
//     all.close();
// }


void exportColocalizationsToJSONByDisease(
    const std::map<std::tuple<int,int,int>, std::set<Timepoint>>& colocalizationByIndividual,
    const std::map<int, std::string>& patientToDiseaseMap,
    const std::string& jsonOutputPath  // path to the final JSON file
) {
    std::map<std::string, std::map<std::string, std::map<std::string, int>>> diseaseColocCounts;

    // Store metadata for each pair so we can export it later
    std::map<std::string, std::string> pairToArgResistanceGroup;
    std::map<std::string, std::string> pairToArgResistance;
    std::map<std::string, std::string> pairToMgeGroup;

    // Build counts by disease → colocalization → status
    for (const auto& [tuple, tps] : colocalizationByIndividual) {
        const int patientID = std::get<0>(tuple);
        const int argID     = std::get<1>(tuple);
        const int mgeID     = std::get<2>(tuple);
        std::string disease = patientToDiseaseMap.at(patientID);

        bool hasDonor = std::any_of(tps.begin(), tps.end(), [](const Timepoint& tp){ return isDonor(tp); });
        bool hasPre   = std::any_of(tps.begin(), tps.end(), [](const Timepoint& tp){ return isPreFMT(tp); });
        bool hasPost  = std::any_of(tps.begin(), tps.end(), [](const Timepoint& tp){ return isPostFMT(tp); });

        std::string status;
        if (hasPost && !hasPre && !hasDonor) status = "emerged";
        else if (hasPre && !hasPost && !hasDonor) status = "disappeared";
        else if (hasPre && !hasPost && hasDonor) status = "disappeared";
        else if (hasDonor && hasPost && !hasPre) status = "transferred";
        else if (hasPre && hasPost) status = "persisted";
        else continue; // skip other patterns

        std::string pairName = getARGName(argID) + "–" + getMGEName(mgeID);
        diseaseColocCounts[disease][pairName][status]++;
        // Record metadata once per pair
        pairToArgResistanceGroup[pairName] = getARGGroupName(argID);
        pairToArgResistance[pairName] = getARGResistance(argID);
        pairToMgeGroup[pairName] = getMGEGroupName(mgeID);
    }

    // Build JSON structure
    json rootJson = json::object();  // use object instead of array

    for (const auto& [disease, colocMap] : diseaseColocCounts) {
        json diseaseArray = json::array();

        for (const auto& [pair, statusMap] : colocMap) {
            for (const auto& [status, count] : statusMap) {
                diseaseArray.push_back({
                    {"colocalization", pair},
                    {"argResistance", pairToArgResistanceGroup[pair]},
                    {"argResistanceMap", pairToArgResistance[pair]},
                    {"mgeGroup", pairToMgeGroup[pair]},
                    {"status", status},
                    {"patients", count}
                });
            }
        }

        rootJson[disease] = diseaseArray;
    }

    std::filesystem::path outFile(jsonOutputPath);
    std::ofstream all(outFile.string());  
    if (!all.is_open()) {
        throw std::runtime_error("Failed to open output JSON file: " + outFile.string());
    }
    all << std::setw(2) << rootJson << std::endl;
    all.close();
}
