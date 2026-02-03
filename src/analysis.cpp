// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <set>
#include "graph.h"
#include "Timepoint.h"
#include "analysis.h"
#include "id_maps.h"
#include "config_loader.h"
#include "traversal.h"
#include <filesystem>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <regex>

std::string getMGEGroupName(int id);
namespace fs = std::filesystem;
fs::path temporal_dynamics_emerge;
fs::path temporal_dynamics_disappear;
fs::path temporal_dynamics_transfer;
fs::path temporal_dynamics_persist;
fs::path top_entities_mge;
fs::path top_entities_arg;
fs::path disease_type_path;
fs::path mge_group_path;

Config cfg = loadConfig("config/paths.json");


bool isPostFMT(const Timepoint& tp) {
    return toString(tp).find("post") != std::string::npos;
}

bool isPreFMT(const Timepoint& tp) {
    return toString(tp).find("pre") != std::string::npos;
}

bool isDonor(const Timepoint& tp) {
    return toString(tp).find("donor") != std::string::npos;
}

static inline bool isPostBin1(Timepoint tp) { int v = static_cast<int>(tp); return v >= 1  && v <= 30; }
static inline bool isPostBin2(Timepoint tp) { int v = static_cast<int>(tp); return v >= 31 && v <= 60; }
static inline bool isPostBin3(Timepoint tp) { int v = static_cast<int>(tp); return v >= 61; }

// In Post stages, highest bin wins if multiple present
static inline int postBinOf(const std::set<Timepoint>& s) {
    bool b1 = false, b2 = false, b3 = false;
    for (auto tp : s) {
        if (!b1 && isPostBin1(tp)) b1 = true;
        else if (!b2 && isPostBin2(tp)) b2 = true;
        else if (!b3 && isPostBin3(tp)) b3 = true;
        if (b1 && b2 && b3) break;
    }
    if (b3) return 3; // post_3 (61+)
    if (b2) return 2; // post_2 (31–60)
    if (b1) return 1; // post_1 (1–30)
    return 0;         // no post
}

/********************************* Patientwise Colocalizations ********************************/
void getPatientwiseColocalizationsByCriteria(
    const Graph& graph,
    const std::map<std::tuple<int, int, int>, std::set<Timepoint>>& colocalizationByIndividual,
    bool donorStatus,
    bool preFMTStatus,
    bool postFMTStatus,
    const std::string& label,
    const std::string& csvFile,
    bool append
) {
    std::map<std::tuple<int, int, int>, std::set<Timepoint>> filteredColocs;

    for (const auto& [tuple, tps] : colocalizationByIndividual) {
        bool hasDonor = std::any_of(tps.begin(), tps.end(), [](const Timepoint& tp) {
            return tp == Timepoint::Donor;
        });
        bool hasPreFMT = std::any_of(tps.begin(), tps.end(), [](const Timepoint& tp) {
            return isPreFMT(tp);
        });
        bool hasPostFMT = std::any_of(tps.begin(), tps.end(), [](const Timepoint& tp) {
            return isPostFMT(tp);
        });

        // Match against provided pattern
        if ((hasDonor == donorStatus) && (hasPreFMT == preFMTStatus) && (hasPostFMT == postFMTStatus)) {
            filteredColocs.insert({tuple, tps});
        }
    }

    std::cout << "Colocalizations (" << label << "): " << filteredColocs.size() << "\n";
    // getTopARGMGEPairsByFrequency(filteredColocs);


    if (!csvFile.empty()) {
        writeColocalizationsToCSV(filteredColocs, csvFile, label, append);
    }
}


/********************************* Prominent Entities ********************************/
void getTopARGMGEPairsByFrequency(
    const std::map<std::tuple<int, int, int>, std::set<Timepoint>>& colocalizations,
    int topN // default: print all
) {
    std::map<std::pair<int, int>, int> countMap;

    // Count (ARG_ID, MGE_ID) occurrences
    for (const auto& [tuple, tps] : colocalizations) {
        int argID = std::get<0>(tuple); // Correct index: ARG is first in tuple
        int mgeID = std::get<1>(tuple); // MGE is second in tuple
        countMap[{argID, mgeID}]++;
    }

    // Convert to vector for sorting
    std::vector<std::pair<std::pair<int, int>, int>> freqList(countMap.begin(), countMap.end());
    std::sort(freqList.begin(), freqList.end(), [](const auto& a, const auto& b) {
        return a.second > b.second;
    });

    std::cout << "Top ARG–MGE pairs by frequency:\n";
    int countPrinted = 0;
    for (const auto& [pair, count] : freqList) {
        if (topN > 0 && countPrinted >= topN) break;

        int argID = pair.first;
        int mgeID = pair.second;

        std::cout << "ARG: ";
        if (argIdMap.count(argID))
            std::cout << getARGName(argID) << " (" << getARGGroupName(argID) << ")";
        else
            std::cout << "Unknown ARG ID " << argID;

        std::cout << ", MGE: ";
        if (mgeIdMap.count(mgeID))
            std::cout << getMGEName(mgeID);
        else
            std::cout << "Unknown MGE ID " << mgeID;

        std::cout << ", Count: " << count << "\n";
        countPrinted++;
    }

    std::cout << "Total unique ARG–MGE pairs: " << freqList.size() << "\n";
}


void getTopARGMGEPairsByFrequencyWODonor(
    const std::map<std::tuple<int, int, int>, std::set<Timepoint>>& colocalizations,
    int topN,
    const std::map<int, std::string>& patientToDiseaseMap, const std::string& top_colocalizations_output)
{
    std::map<std::pair<int, int>, int> countMap;
    std::map<std::pair<int, int>, std::map<std::string, int>> diseaseCountMap;
    std::set<std::string> diseaseSet;

    for (const auto& [tuple, tps] : colocalizations) {
        int patientID = std::get<0>(tuple);
        int argID     = std::get<1>(tuple);
        int mgeID     = std::get<2>(tuple);

        bool hasNonDonor = std::any_of(tps.begin(), tps.end(),
            [](const Timepoint& tp){ return tp != Timepoint::Donor; });
        if (!hasNonDonor) continue;

        auto key = std::make_pair(argID, mgeID);
        countMap[key]++;

        std::string disease = "Unknown";
        if (patientToDiseaseMap.count(patientID))
            disease = patientToDiseaseMap.at(patientID);

        diseaseCountMap[key][disease]++;
        diseaseSet.insert(disease);
    }

    // Sort diseases for consistent column order
    std::vector<std::string> diseases(diseaseSet.begin(), diseaseSet.end());

    // Build header
    std::vector<std::string> header = {
        "ARG_Name","ARG_Group","MGE_Name","TotalCount"
    };
    header.insert(header.end(), diseases.begin(), diseases.end());

    // Sort ARG-MGE pairs by total count
    std::vector<std::pair<std::pair<int,int>,int>> freqList(countMap.begin(), countMap.end());
    std::sort(freqList.begin(), freqList.end(),
              [](auto& a, auto& b){ return a.second > b.second; });

    std::vector<std::vector<std::string>> rows;
    int printed = 0;

    for (const auto& [pair,totalCount] : freqList) {
        if (topN > 0 && printed >= topN) break;

        int argID = pair.first;
        int mgeID = pair.second;

        std::vector<std::string> row = {
            getARGName(argID),
            getARGGroupName(argID),
            getMGEName(mgeID),
            std::to_string(totalCount)
        };

        for (const auto& dz : diseases) {
            int c = 0;
            if (diseaseCountMap[pair].count(dz))
                c = diseaseCountMap[pair].at(dz);

            row.push_back(std::to_string(c));
        }

        rows.push_back(row);
        printed++;
    }

    writeCSV(top_colocalizations_output, header, rows);
}


void mostProminentEntities(const Graph& g) {
    std::map<std::tuple<int, int, int>, std::set<Timepoint>> colocalizationByIndividual;
    traverseGraph(g, colocalizationByIndividual);


    std::vector<std::pair<int, int>> topARGs = getTopKEntities(g, true, static_cast<unsigned int>(10)); // Top 10 ARGs
    {
    std::vector<std::vector<std::string>> rows;
    for (auto& [id,count] : topARGs) {
        rows.push_back({
            std::to_string(id),
            getARGName(id),
            getARGGroupName(id),
            std::to_string(count)
        });
    }

    writeCSV(cfg.output_top_arg,
        {"ARG_ID","ARG_Name","ARG_Group","Count"},
        rows);
    }

    std::vector<std::pair<int, int>> topMGEs = getTopKEntities(g, false, static_cast<unsigned int>(10)); // Top 10 MGEs 
    {
    std::vector<std::vector<std::string>> rows;
    for (auto& [id,count] : topMGEs) {
        rows.push_back({
            std::to_string(id),
            getMGEName(id),
            std::to_string(count)
        });
    }

    writeCSV(cfg.output_top_mge,
        {"MGE_ID","MGE_Name","Count"},
        rows);
    }
}

/***************************************** Connected MGEs *********************************************/
void getConnectedMGEs(const Graph& graph, int argID) {
    std::unordered_set<int> connectedMGEs;

    for (const Edge& edge : graph.edges) {
        if (!edge.isColo) continue;
        if (edge.source.isARG && edge.source.id == argID && !edge.target.isARG) {
            connectedMGEs.insert(edge.target.id);
        } else if (edge.target.isARG && edge.target.id == argID && !edge.source.isARG) {
            connectedMGEs.insert(edge.source.id);
        }
    }

    std::cout << "ARG ID " << argID << " is connected to MGE IDs:\n";
    for (int mge : connectedMGEs) {
        std::cout << "  - MGE " << mge;
        if (mgeIdMap.count(mge)) std::cout << " (" << mgeIdMap.at(mge) << ")";
        std::cout << "\n";
    }
}



/***************************************** Write Functions *********************************************/

/* Write temporal dynamics counts for a specific disease */
void writeTemporalDynamicsCountsForDisease(
    const std::string& disease,
    std::map<std::tuple<int,int,int>,std::set<Timepoint>>& colocalizationByIndividual,
    const std::map<int,std::string>& patientToDiseaseMap)
{
    std::map<std::tuple<int,int,int,int,int>,int> comboCounts;

    for (auto& [key,tps] : colocalizationByIndividual) {
        int patientID = std::get<0>(key);
        if (!patientToDiseaseMap.count(patientID) ||
            patientToDiseaseMap.at(patientID) != disease) continue;

        int argID = std::get<1>(key);
        int mgeID = std::get<2>(key);

        bool donor = std::any_of(tps.begin(),tps.end(),isDonor);
        bool pre   = std::any_of(tps.begin(),tps.end(),isPreFMT);
        int post   = postBinOf(tps);

        comboCounts[{argID,mgeID,donor,pre,post}]++;
    }

    std::vector<std::vector<std::string>> rows;
    for (auto& [k,cnt] : comboCounts) {
        int argID,mgeID,donor,pre,post;
        std::tie(argID,mgeID,donor,pre,post)=k;

        rows.push_back({
            getARGName(argID),
            getMGEName(mgeID),
            std::to_string(donor),
            std::to_string(pre),
            std::to_string(post),
            std::to_string(cnt)
        });
    }

    writeCSV("viz/output/disease_type/" + disease + ".csv",
        {"ARG_ID","MGE_ID","Donor","Pre","Post","PatientCount"},
        rows);
}


/* write a CSV for all diseases */
void writeAllDiseasesTemporalDynamicsCounts(
    std::map<std::tuple<int, int, int>, std::set<Timepoint>>& colocalizationByIndividual,
    const std::map<int, std::string>& patientToDiseaseMap
) {
    std::set<std::string> diseases;
    for (const auto& [pid, dz] : patientToDiseaseMap) diseases.insert(dz);
    for (const auto& dz : diseases)
        writeTemporalDynamicsCountsForDisease(dz, colocalizationByIndividual, patientToDiseaseMap);
}


/* Write temporal dynamics counts for a specific MGE group */
void writeTemporalDynamicsCountsForMGEGroup(const std::map<std::tuple<int,int,int>,std::set<Timepoint>>& colocalizationByIndividual){
    std::unordered_map<std::string,
    std::map<std::tuple<int,int,int,int,int>,int>> groupedCounts;

    for (auto& [key,tps] : colocalizationByIndividual) {
        int argID = std::get<1>(key);
        int mgeID = std::get<2>(key);
        std::string group = getMGEGroupName(mgeID);

        bool donor = std::any_of(tps.begin(),tps.end(),isDonor);
        bool pre   = std::any_of(tps.begin(),tps.end(),isPreFMT);
        bool post  = std::any_of(tps.begin(),tps.end(),isPostFMT);

        groupedCounts[group][{argID,mgeID,donor,pre,post}]++;
    }

    for (auto& [group,comboCounts] : groupedCounts) {
        std::vector<std::vector<std::string>> rows;
        for (auto& [k,cnt] : comboCounts) {
            int argID,mgeID,donor,pre,post;
            std::tie(argID,mgeID,donor,pre,post)=k;

            rows.push_back({
                getARGName(argID),
                getMGEName(mgeID),
                std::to_string(donor),
                std::to_string(pre),
                std::to_string(post),
                std::to_string(cnt)
            });
        }
// remove any filesystem-unfriendly characters not just beginning and end
        std::string filename = group;  // copy, modifiable
        filename = std::regex_replace(filename,std::regex(R"([\/\\:\*\?"<>|])"), "_");

        writeCSV(cfg.output_mge_group + "/"+ filename + ".csv",
            {"ARG_ID","MGE_ID","Donor","Pre","Post","PatientCount"},
            rows);
    }

}

/* Write colocalizations to a CSV file */
void writeColocalizationsToCSV(
    const std::map<std::tuple<int,int,int>, std::set<Timepoint>>& colocs,
    const std::string& filename,
    const std::string& label,
    bool append)
{
    std::map<std::pair<int,int>, std::set<int>> aggregated;

    for (const auto& [tuple,tps] : colocs) {
        int argId = std::get<1>(tuple);
        int mgeId = std::get<2>(tuple);
        int patientId = std::get<0>(tuple);
        aggregated[{argId,mgeId}].insert(patientId);
    }

    std::vector<std::vector<std::string>> rows;
    for (auto& [pair,patients] : aggregated) {
        rows.push_back({
            getARGName(pair.first),
            getMGEName(pair.second),
            std::to_string(patients.size()),
            label
        });
    }

    writeCSV(filename,
        {"ARG_Name","MGE_Name","PatientCount","Label"},
        rows,
        append);
}


/* Export colocalizations to seperate files based on temporal dynamics */
void exportColocalizations(const Graph& g,
    const std::map<std::tuple<int,int,int>, std::set<Timepoint>>& colocalizationByIndividual) 
{
    temporal_dynamics_emerge = fs::path(cfg.output_emerge);
    temporal_dynamics_transfer = fs::path(cfg.output_transfer);
    temporal_dynamics_disappear = fs::path(cfg.output_disappear);
    temporal_dynamics_persist = fs::path(cfg.output_persist);
    // Emerge
    
    getPatientwiseColocalizationsByCriteria(g, colocalizationByIndividual,
        false, false, true, "PostFMT Only", temporal_dynamics_emerge.string());
    // Disappear
    getPatientwiseColocalizationsByCriteria(g, colocalizationByIndividual,
        false, true, false, "PreFMT Only", temporal_dynamics_disappear.string(), false);
    getPatientwiseColocalizationsByCriteria(g, colocalizationByIndividual,
        true, true, false, "Donor & PreFMT Only", temporal_dynamics_disappear.string(), true);
    // Transfer
    getPatientwiseColocalizationsByCriteria(g, colocalizationByIndividual,
        true, false, true, "Donor & PostFMT Only", temporal_dynamics_transfer.string(), true);
    // Persist
    getPatientwiseColocalizationsByCriteria(g, colocalizationByIndividual,
        false, true, true, "PreFMT & PostFMT Only", temporal_dynamics_persist.string(), false);
    getPatientwiseColocalizationsByCriteria(g, colocalizationByIndividual,
        true, true, true, "PreFMT, Donor & PostFMT", temporal_dynamics_persist.string(), true);

    exportDetailedTemporalDynamics(colocalizationByIndividual);
}


/* Write graph node/edge counts to a CSV file */
void writeGraphStatisticsCSV(
    const Graph& g,
    const std::unordered_map<Node,std::unordered_set<Node>>& adjacency,
    const std::string& filename)
{
    size_t total_nodes = g.nodes.size();
    size_t total_edges = g.edges.size();
    size_t arg_count   = std::count_if(g.nodes.begin(),g.nodes.end(),[](auto& n){return n.isARG;});
    size_t mge_count   = total_nodes - arg_count;
    size_t colo_edges  = std::count_if(g.edges.begin(),g.edges.end(),[](auto& e){return e.isColo;});
    size_t temporal_edges = total_edges - colo_edges;

    writeCSV(filename,
        {"TotalNodes","TotalEdges","ARGs","MGEs","ColocalizationEdges","TemporalEdges","AdjacencyNodes"},
        {{
            std::to_string(total_nodes),
            std::to_string(total_edges),
            std::to_string(arg_count),
            std::to_string(mge_count),
            std::to_string(colo_edges),
            std::to_string(temporal_edges),
            std::to_string(adjacency.size())
        }});
}


/* Generic CSV writing function */
void writeCSV(
    const std::string& filename,
    const std::vector<std::string>& header,
    const std::vector<std::vector<std::string>>& rows,
    bool append)
{
    std::ofstream out;
    if (append)
        out.open(filename, std::ios::app);
    else
        out.open(filename);

    if (!out.is_open()) {
        throw std::runtime_error("Failed to open CSV: " + filename);
    }

    if (!append) {
        for (size_t i = 0; i < header.size(); i++)
            out << header[i] << (i+1<header.size() ? "," : "\n");
    }

    for (const auto& row : rows) {
        for (size_t i = 0; i < row.size(); i++)
            out << row[i] << (i+1<row.size() ? "," : "\n");
    }
}



/***************************************** Not used in method********************************************/
void analyzeColocalizations(const Graph& g, 
                            const std::unordered_map<Node, std::unordered_set<Node>>& adjacency) {

    std::map<std::tuple<int, int, int>, std::set<Timepoint>> colocalizationByIndividual;
    traverseGraph(g, colocalizationByIndividual);

    std::cout << "Patientwise Colocalization dynamics over time:\n";
    getPatientwiseColocalizationsByCriteria(g, colocalizationByIndividual, false, false, true, "PostFMT Only");
    getPatientwiseColocalizationsByCriteria(g, colocalizationByIndividual, false, true, false, "PreFMT Only");
    getPatientwiseColocalizationsByCriteria(g, colocalizationByIndividual, true, false, false, "Donor Only");
    getPatientwiseColocalizationsByCriteria(g, colocalizationByIndividual, true, true, false, "Donor & PreFMT Only");
    getPatientwiseColocalizationsByCriteria(g, colocalizationByIndividual, false, true, true, "PreFMT & PostFMT Only");
    getPatientwiseColocalizationsByCriteria(g, colocalizationByIndividual, true, true, true, "PreFMT, Donor & PostFMT");
    getPatientwiseColocalizationsByCriteria(g, colocalizationByIndividual, true, false, true, "Donor & PostFMT Only");
}


// Used in Graph Analysis, Colocalizations collectively global scale. 
/***************************************** Collective Analysis of Colocalizations ****************************************/



void analyzeColocalizationsCollectively(const Graph& g, 
                                         const std::unordered_map<Node, std::unordered_set<Node>>& adjacency) {
    std::map<std::tuple<int, int, int>, std::set<Timepoint>> colocalizationByIndividual;
    traverseGraph(g, colocalizationByIndividual);
    std::map<std::pair<int, int>, std::set<int>> globalPairToPatients;
    std::cout << "Colocalization dynamics over time:\n";
    getColocalizationsByCriteria(colocalizationByIndividual, false, false, true, globalPairToPatients);
    globalPairToPatients.clear(); // Clear the map for the next scenario
    getColocalizationsByCriteria(colocalizationByIndividual, false, true, false, globalPairToPatients);
    globalPairToPatients.clear(); // Clear the map for the next scenario
    getColocalizationsByCriteria(colocalizationByIndividual, true, false, false, globalPairToPatients);
    globalPairToPatients.clear(); // Clear the map for the next scenario
    getColocalizationsByCriteria(colocalizationByIndividual, true, true, false, globalPairToPatients);
    globalPairToPatients.clear(); // Clear the map for the next scenario
    getColocalizationsByCriteria(colocalizationByIndividual, false, true, true, globalPairToPatients);
    globalPairToPatients.clear(); // Clear the map for the next scenario
    getColocalizationsByCriteria(colocalizationByIndividual, true, true, true, globalPairToPatients);
    globalPairToPatients.clear(); // Clear the map for the next scenario
    getColocalizationsByCriteria(colocalizationByIndividual, true, false, true, globalPairToPatients);
    globalPairToPatients.clear(); // Clear the map for the next scenario
}


void getColocalizationsByCriteria(
    const std::map<std::tuple<int, int, int>, std::set<Timepoint>>& colocalizationByIndividual,
    bool donorStatus,
    bool preFMTStatus,
    bool postFMTStatus,
    std::map<std::pair<int, int>, std::set<int>>& globalPairToPatients
) {
    // Step 1: Determine globally valid (ARG, MGE) pairs
    std::map<std::pair<int, int>, std::set<Timepoint>> globalTimepoints;

    for (const auto& [tuple, tps] : colocalizationByIndividual) {
        int argID = std::get<1>(tuple);
        int mgeID = std::get<2>(tuple);
        globalTimepoints[{argID, mgeID}].insert(tps.begin(), tps.end());
    }

    std::set<std::pair<int, int>> validPairs;
    for (const auto& [pair, tps] : globalTimepoints) {
        bool hasDonor   = std::any_of(tps.begin(), tps.end(), [](const Timepoint& tp){ return isDonor(tp); });
        bool hasPreFMT  = std::any_of(tps.begin(), tps.end(), [](const Timepoint& tp){ return isPreFMT(tp); });
        bool hasPostFMT = std::any_of(tps.begin(), tps.end(), [](const Timepoint& tp){ return isPostFMT(tp); });

        if (hasDonor == donorStatus &&
            hasPreFMT == preFMTStatus &&
            hasPostFMT == postFMTStatus) {
            validPairs.insert(pair);
        }
    }

    // Step 2: Count patients per valid (ARG, MGE) pair only if their individual timepoints match the scenario
    for (const auto& [tuple, tps] : colocalizationByIndividual) {
        int patientID = std::get<0>(tuple);
        int argID     = std::get<1>(tuple);
        int mgeID     = std::get<2>(tuple);
        std::pair<int, int> pair = {argID, mgeID};

        if (!validPairs.count(pair)) continue;

        bool hasDonor   = std::any_of(tps.begin(), tps.end(), [](const Timepoint& tp) { return isDonor(tp); });
        bool hasPreFMT  = std::any_of(tps.begin(), tps.end(), [](const Timepoint& tp) { return isPreFMT(tp); });
        bool hasPostFMT = std::any_of(tps.begin(), tps.end(), [](const Timepoint& tp) { return isPostFMT(tp); });

        if (hasDonor == donorStatus &&
            hasPreFMT == preFMTStatus &&
            hasPostFMT == postFMTStatus) {
            globalPairToPatients[pair].insert(patientID);
        }
    }

    getTopARGMGEPairsByUniquePatients(globalPairToPatients, 10, "Filtered by scenario");
}


void getTopARGMGEPairsByUniquePatients(
    const std::map<std::pair<int, int>, std::set<int>>& globalPairToPatients,
    int topN,
    const std::string& label
) {
    std::vector<std::pair<std::pair<int, int>, int>> freqList;
    for (const auto& [pair, patientSet] : globalPairToPatients) {
        freqList.emplace_back(pair, patientSet.size());
    }

    std::sort(freqList.begin(), freqList.end(), [](const auto& a, const auto& b) {
        return a.second > b.second;
    });

    std::cout << "\nTop colocalizations by unique patients (" << label << "):\n";
    int countPrinted = 0;
    for (const auto& [pair, count] : freqList) {
        if (topN > 0 && countPrinted >= topN) break;
        int argID = pair.first;
        int mgeID = pair.second;

        std::cout << "ARG: ";
        if (argIdMap.count(argID)) std::cout << getARGName(argID) << " (" << getARGGroupName(argID) << ")";
        else std::cout << "Unknown ARG ID " << argID;

        std::cout << ", MGE: ";
        if (mgeIdMap.count(mgeID)) std::cout << getMGEName(mgeID);
        else std::cout << "Unknown MGE ID " << mgeID;

        std::cout << ",Patients: " << count << "\n";
    }

    // std::cout << "Total unique colocalizations: " << freqList.size() << "\n";
}

void writeDetailedCSV(
    const std::map<std::tuple<int, int, int>, std::set<Timepoint>>& colocs,
    const std::string& filename,
    const std::string& label,
    bool append
) {
    std::ofstream file;
    if (append) {
        file.open(filename, std::ios::app);
    } else {
        file.open(filename);
        file << "ARG_Name,MGE_Name,Patient_ID,Label\n"; 
    }

    for (const auto& [tuple, tps] : colocs) {
        int patientId = std::get<0>(tuple);
        int argId     = std::get<1>(tuple);
        int mgeId     = std::get<2>(tuple);

        std::string argName = getARGName(argId);
        std::string mgeName = getMGEName(mgeId);

        if (argName.empty()) argName = "Unknown_ARG_" + std::to_string(argId);
        if (mgeName.empty()) mgeName = "Unknown_MGE_" + std::to_string(mgeId);

        file << argName << ","
             << mgeName << ","
             << patientId << "," 
             << label << "\n";
    }
    file.close();
    std::cout << "Detailed list written to " << filename << "\n";
}

void getDetailedColocalizationsByCriteria(
    const std::map<std::tuple<int, int, int>, std::set<Timepoint>>& colocalizationByIndividual,
    bool donorStatus,
    bool preFMTStatus,
    bool postFMTStatus,
    const std::string& label,
    const std::string& csvFile,
    bool append
) {
    std::map<std::tuple<int, int, int>, std::set<Timepoint>> filteredColocs;

    for (const auto& [tuple, tps] : colocalizationByIndividual) {
        bool hasDonor = std::any_of(tps.begin(), tps.end(), [](const Timepoint& tp) {
            return tp == Timepoint::Donor;
        });
        bool hasPreFMT = std::any_of(tps.begin(), tps.end(), [](const Timepoint& tp) {
            return isPreFMT(tp);
        });
        bool hasPostFMT = std::any_of(tps.begin(), tps.end(), [](const Timepoint& tp) {
            return isPostFMT(tp);
        });

        if ((hasDonor == donorStatus) && (hasPreFMT == preFMTStatus) && (hasPostFMT == postFMTStatus)) {
            filteredColocs.insert({tuple, tps});
        }
    }

    if (!csvFile.empty()) {
        writeDetailedCSV(filteredColocs, csvFile, label, append);
    }
}

void exportDetailedTemporalDynamics(
    const std::map<std::tuple<int,int,int>, std::set<Timepoint>>& colocalizationByIndividual) 
{
    auto makeDetailedPath = [](fs::path p) {
        std::string s = p.string();
        size_t dot = s.find_last_of('.');
        if(dot != std::string::npos) s.insert(dot, "_detailed");
        else s += "_detailed";
        return s;
    };

    // Emerge
    getDetailedColocalizationsByCriteria(colocalizationByIndividual,
        false, false, true, "PostFMT Only", makeDetailedPath(temporal_dynamics_emerge), false);
    
    // Disappear
    getDetailedColocalizationsByCriteria(colocalizationByIndividual,
        false, true, false, "PreFMT Only", makeDetailedPath(temporal_dynamics_disappear), false);
    getDetailedColocalizationsByCriteria(colocalizationByIndividual,
        true, true, false, "Donor & PreFMT Only", makeDetailedPath(temporal_dynamics_disappear), true);

    // Transfer
    getDetailedColocalizationsByCriteria(colocalizationByIndividual,
        true, false, true, "Donor & PostFMT Only", makeDetailedPath(temporal_dynamics_transfer), true);

    // Persist
    getDetailedColocalizationsByCriteria(colocalizationByIndividual,
        false, true, true, "PreFMT & PostFMT Only", makeDetailedPath(temporal_dynamics_persist), false);
    getDetailedColocalizationsByCriteria(colocalizationByIndividual,
        true, true, true, "PreFMT, Donor & PostFMT", makeDetailedPath(temporal_dynamics_persist), true);
        
    std::cout << "\n[SUCCESS] Detailed datasets with Patient IDs have been generated.\n";
}



