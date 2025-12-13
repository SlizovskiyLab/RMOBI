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


static std::string getTimepointColor(const Timepoint& tp) {
    int timeValue = static_cast<int>(tp);
    if (timeValue == 1000) return "yellow";
    if (timeValue == 0)    return "red";
    if (timeValue > 0 && timeValue < 31) return "#99D2FF";
    if (timeValue > 30 && timeValue < 61) return "#4D9DFF";
    if (timeValue > 60)    return "#3A6EFF";
    return "green"; // fallback
}

static std::string getTimepointCategory(const Timepoint& tp) {
    int timeValue = static_cast<int>(tp);
    if (timeValue == 1000) return "donor";
    if (timeValue == 0)    return "pre";
    if (timeValue > 0 && timeValue < 31) return "post1";
    if (timeValue > 30 && timeValue < 61) return "post2";
    if (timeValue > 60)    return "post3";
    return "unknown"; // Fallback for any unexpected values
}


auto timepointOrder = [](Timepoint tp) -> int {
    if (tp == Timepoint::Donor) return -1;
    if (tp == Timepoint::PreFMT) return 0;
    return static_cast<int>(tp);
};

std::string getLabel(const Node& node) {
    std::string label = node.isARG ? getARGName(node.id) : getMGENameForLabel(node.id);
    return label;
}


bool exportGraphToJsonSimple(const Graph& g, const std::string& outPathStr, const std::map<int, std::string>& patientToDiseaseMap) {
    json j;
    j["nodes"] = json::array();
    j["links"] = json::array();

    std::unordered_set<Node> active_nodes;
    std::set<std::pair<Node, Node>> processedColoEdges;

    for (const Edge& edge : g.edges) {
        if (edge.source == edge.target) continue;

        active_nodes.insert(edge.source);
        active_nodes.insert(edge.target);
        
        std::string style;
        std::string color;
        double penwidth = 4.0;
        std::string type = "other";
        json diseases = json::array();
        if (edge.isColo) {
            auto canon = std::minmax(edge.source, edge.target);
            if (processedColoEdges.count(canon)) continue;
            processedColoEdges.insert(canon);
            style = "solid";
            color = "#696969";
            type  = "colocalization";
            std::set<std::string> diseaseSet; 
            for (int patientID : edge.individuals) {
                auto it = patientToDiseaseMap.find(patientID);
                if (it != patientToDiseaseMap.end()) {
                    diseaseSet.insert(it->second); 
                }
            }
            for (const auto& diseaseName : diseaseSet) {
                diseases.push_back(diseaseName);
            }
            int count = static_cast<int>(edge.individuals.size());
            if (count > 1) penwidth = 4.0 + (count - 1) * 2.0;
            penwidth = std::min(10.0, penwidth);
        } 
        else if (!edge.isColo) {
            style = "dashed";
            type  = "temporal";
            int w = edge.weight;
            if (w > 1) penwidth = 4.0 + (w - 1) * 2.0;
            penwidth = std::min(10.0, penwidth);
            Timepoint src_tp = edge.source.timepoint;
            Timepoint tgt_tp = edge.target.timepoint;
            bool tgt_is_post = (tgt_tp != Timepoint::Donor && tgt_tp != Timepoint::PreFMT);
            if (src_tp == Timepoint::Donor && tgt_tp == Timepoint::PreFMT)      color = "#006400";
            else if (src_tp == Timepoint::Donor && tgt_is_post)                 color = "#4B0082";
            else if (src_tp == Timepoint::PreFMT && tgt_is_post)                color = "orange";
            else                                                                color = "black";
        } 
        else {
            style = "solid";
            color = "#808080";
        }
        j["links"].push_back({
            {"source", getNodeName(edge.source)},
            {"target", getNodeName(edge.target)},
            {"individualCount", static_cast<int>(edge.individuals.size())},
            {"style", style},
            {"color", color},
            {"penwidth", penwidth},
            {"isColo", edge.isColo},
            {"type", type},
            {"diseases", diseases}
        });
    }

    for (const Node& n : active_nodes) {
        std::string shape;
        std::string mgeGroup = ""; 

        if (n.isARG) {
            shape = "circle";
        } else {
            mgeGroup = getMGEGroupName(n.id);
            shape = getMGEGroupShape(mgeGroup);
        }

        j["nodes"].push_back({
            {"id",                getNodeName(n)},
            {"label",             getLabel(n)},
            {"isARG",             n.isARG},
            {"timepoint",         static_cast<int>(n.timepoint)},
            {"color",             getTimepointColor(n.timepoint)},
            {"shape",             shape},
            {"mgeGroup",          mgeGroup},
            {"timepointCategory", getTimepointCategory(n.timepoint)}
        });
    }

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


bool exportParentGraphToJson(const Graph& g, const std::string& outPathStr, const std::map<int, std::string>& patientToDiseaseMap, bool showLabels) {
    json j;
    j["nodes"] = json::array();
    j["links"] = json::array();

    struct ParentNodeInfo {
        std::string name;
        Timepoint tp;
        int argId;
        int mgeId;
    };

    int colocCounter = 0;
    std::map<std::tuple<int,int,Timepoint>, std::string> uniqueParents;
    std::map<std::pair<int,int>, std::vector<ParentNodeInfo>> colocMap;

    for (const Edge& edge : g.edges) {
        if (!edge.isColo) continue;

        const Node& argNode = edge.source.isARG ? edge.source : edge.target;
        const Node& mgeNode = edge.source.isARG ? edge.target : edge.source;

        int argId = argNode.id;
        int mgeId = mgeNode.id;
        Timepoint tp = argNode.timepoint;

        auto key = std::make_tuple(argId, mgeId, tp);
        if (!uniqueParents.count(key)) {
            std::string parentName = "Parent_" + std::to_string(++colocCounter);
            uniqueParents[key] = parentName;

            std::string color = getTimepointColor(tp);
            std::string groupName = getMGEGroupName(mgeId);
            std::string shape = getMGEGroupShape(groupName);
            std::string label = showLabels ? (getARGName(argId) + "+" + getMGENameForLabel(mgeId)) : "";
            
            std::set<std::string> diseaseSet;
            std::map<std::string, std::set<int>> diseaseToIndividualsLocal; // local per colocalization

            // collect diseases + individuals for this edge only
            for (int patientID : edge.individuals) {
                auto it = patientToDiseaseMap.find(patientID);
                if (it != patientToDiseaseMap.end()) {
                    diseaseSet.insert(it->second);
                    diseaseToIndividualsLocal[it->second].insert(patientID);
                }
            }

            // convert diseases list
            json diseases = json::array();
            for (const auto& diseaseName : diseaseSet) {
                diseases.push_back(diseaseName);
            }

            // convert to JSON: disease → patient count (for this colocalization/timepoint)
            json diseaseCounts = json::object();
            for (const auto& [diseaseName, individuals] : diseaseToIndividualsLocal) {
                diseaseCounts[diseaseName] = static_cast<int>(individuals.size());
            }


            j["nodes"].push_back({
                {"id",                parentName},
                {"label",             label},
                {"argId",             argId},
                {"mgeId",             mgeId},
                {"timepoint",         static_cast<int>(tp)},
                {"color",             color},
                {"shape",             shape},
                {"diseases",          diseases},
                {"diseaseCounts",     diseaseCounts},
                {"mgeGroup",          groupName}, 
                {"timepointCategory", getTimepointCategory(tp)}
            });
        }

        auto pairKey = std::make_pair(argId, mgeId);
        colocMap[pairKey].push_back({uniqueParents[key], tp, argId, mgeId});
    }

    // ... (The temporal link-processing loop remains unchanged) ...
    for (auto& entry : colocMap) {
        auto& parentNodes = entry.second;
        std::sort(parentNodes.begin(), parentNodes.end(),
            [&](const ParentNodeInfo& a, const ParentNodeInfo& b) {
                return timepointOrder(a.tp) < timepointOrder(b.tp);
            });

        for (size_t i = 0; i + 1 < parentNodes.size(); ++i) {
            if (parentNodes[i].tp == parentNodes[i+1].tp || parentNodes[i].name == parentNodes[i+1].name) {
                continue;
            }
            Timepoint src_tp = parentNodes[i].tp;
            Timepoint tgt_tp = parentNodes[i+1].tp;
            bool tgt_is_post = (tgt_tp != Timepoint::Donor && tgt_tp != Timepoint::PreFMT);
            std::string color;
            if (src_tp == Timepoint::Donor && tgt_tp == Timepoint::PreFMT)      color = "#006400";
            else if (src_tp == Timepoint::Donor && tgt_is_post)                 color = "#4B0082";
            else if (src_tp == Timepoint::PreFMT && tgt_is_post)                color = "orange";
            else                                                                color = "black";

            j["links"].push_back({
                {"source", parentNodes[i].name},
                {"target", parentNodes[i+1].name},
                {"style", "dashed"},
                {"color", color},
                {"penwidth", 5.0},
                {"isColo", false},
                {"type", "temporal"}
            });
        }
    }

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
        else if (hasDonor && hasPost && !hasPre) status = "transferred";
        else if (hasPre && hasPost) status = "persisted";
        else continue; // skip other patterns

        std::string pairName = getARGName(argID) + "–" + getMGEName(mgeID);
        diseaseColocCounts[disease][pairName][status]++;
    }

    // Build JSON structure
    json rootJson = json::object();  // use object instead of array

    for (const auto& [disease, colocMap] : diseaseColocCounts) {
        json diseaseArray = json::array();

        for (const auto& [pair, statusMap] : colocMap) {
            for (const auto& [status, count] : statusMap) {
                diseaseArray.push_back({
                    {"colocalization", pair},
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
