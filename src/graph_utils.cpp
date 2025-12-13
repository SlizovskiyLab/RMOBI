// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

#include "graph_utils.h"
#include "id_maps.h"
#include <iostream>
#include <unordered_set>

Graph filterGraphByARGAndMGENames(const Graph& g, const std::string& argName, const std::string& mgeName) {
    int argID = getARGId(argName);
    if (argID == -1) {
        std::cerr << "ARG not found: " << argName << "\n";
        return {};
    }

    // This now uses the new function that searches both MGE name maps.
    int mgeID = getMGEIdByName(mgeName); 
    if (mgeID == -1) {
        std::cerr << "MGE not found: " << mgeName << "\n";
        return {};
    }

    Graph subgraph;
    std::unordered_set<Node> relevant_nodes;

    // First pass: Find all colocalization edges between the specific ARG and MGE
    // and collect all the nodes involved.
    for (const Edge& edge : g.edges) {
        if (!edge.isColo) continue;

        bool arg_is_source = edge.source.isARG && edge.source.id == argID;
        bool mge_is_target = !edge.target.isARG && edge.target.id == mgeID;

        bool mge_is_source = !edge.source.isARG && edge.source.id == mgeID;
        bool arg_is_target = edge.target.isARG && edge.target.id == argID;

        if ((arg_is_source && mge_is_target) || (mge_is_source && arg_is_target)) {
            subgraph.edges.insert(edge);
            relevant_nodes.insert(edge.source);
            relevant_nodes.insert(edge.target);
        }
    }

    // Second pass: Add all temporal edges that connect the relevant nodes.
    for (const Edge& edge : g.edges) {
        if (edge.isColo) continue; // Skip colocalization edges this time

        // Check if a temporal edge connects two nodes that we've already identified as relevant
        if (relevant_nodes.count(edge.source) && relevant_nodes.count(edge.target)) {
            subgraph.edges.insert(edge);
        }
    }
    
    // Add all the collected nodes to the subgraph
    for(const auto& node : relevant_nodes) {
        subgraph.nodes.insert(node);
    }

    return subgraph;
}

Graph filterGraphByTimepoint(const Graph& g, const std::string& timepointCategory) {
    Graph subgraph;

    auto matchesCategory = [&](Timepoint tp) {
        if (timepointCategory == "donor") {
            return tp == Timepoint::Donor;
        }
        if (timepointCategory == "pre") {
            return tp == Timepoint::PreFMT;
        }
        if (timepointCategory == "post") {
            return tp != Timepoint::Donor && tp != Timepoint::PreFMT;
        }
        return false;
    };

    for (const auto& edge : g.edges) {
        if (matchesCategory(edge.source.timepoint) && matchesCategory(edge.target.timepoint)) {
            subgraph.nodes.insert(edge.source);
            subgraph.nodes.insert(edge.target);
            subgraph.edges.insert(edge);
        }
    }

    if (subgraph.nodes.empty()) {
        std::cerr << "Warning: No nodes found for timepoint category '" << timepointCategory 
                  << "'. The resulting graph will be empty.\n";
    }

    return subgraph;
}

Graph filterGraphByARGName(const Graph& g, const std::string& argName) {
    int argID = getARGId(argName);
    if (argID == -1) {
        std::cerr << "ARG not found: " << argName << "\n";
        return {};
    }

    Graph subgraph;
    std::unordered_set<Node> relevant_nodes;

    // First pass: Find colocalization edges involving the ARG and collect all connected nodes.
    for (const Edge& edge : g.edges) {
        if (!edge.isColo) continue;

        if ((edge.source.isARG && edge.source.id == argID) ||
            (edge.target.isARG && edge.target.id == argID)) {
            
            subgraph.edges.insert(edge);
            relevant_nodes.insert(edge.source);
            relevant_nodes.insert(edge.target);
        }
    }

    // Second pass: Add temporal edges that connect any of the nodes we've collected.
    for (const Edge& edge : g.edges) {
        if (edge.isColo) continue;

        if (relevant_nodes.count(edge.source) && relevant_nodes.count(edge.target)) {
            subgraph.edges.insert(edge);
        }
    }

    // Add all collected nodes to the subgraph's node set.
    for(const auto& node : relevant_nodes) {
        subgraph.nodes.insert(node);
    }

    return subgraph;
}

Graph filterGraphByMGEName(const Graph& g, const std::string& mgeName) {
    // This now uses the new function that searches both MGE name maps.
    int mgeID = getMGEIdByName(mgeName); 
    if (mgeID == -1) {
        std::cerr << "MGE not found: " << mgeName << "\n";
        return {};
    }

    Graph subgraph;
    std::unordered_set<Node> relevant_nodes;

    // First pass: Find colocalization edges involving the MGE and collect all connected nodes.
    for (const Edge& edge : g.edges) {
        if (!edge.isColo) continue;

        if ((!edge.source.isARG && edge.source.id == mgeID) ||
            (!edge.target.isARG && edge.target.id == mgeID)) {

            subgraph.edges.insert(edge);
            relevant_nodes.insert(edge.source);
            relevant_nodes.insert(edge.target);
        }
    }

    // Second pass: Add temporal edges that connect any of the nodes we've collected.
    for (const Edge& edge : g.edges) {
        if (edge.isColo) continue;

        if (relevant_nodes.count(edge.source) && relevant_nodes.count(edge.target)) {
            subgraph.edges.insert(edge);
        }
    }

    // Add all collected nodes to the subgraph's node set.
    for(const auto& node : relevant_nodes) {
        subgraph.nodes.insert(node);
    }

    return subgraph;
}

Graph filterGraphByMGEGroup(const Graph& g, const std::string& groupName) {
    std::unordered_set<int> mgeIDsInGroup;
    for (const auto& [id, name] : mgeGroupMap) {
        if (name == groupName) {
            mgeIDsInGroup.insert(id);
        }
    }

    if (mgeIDsInGroup.empty()) {
        std::cerr << "No MGEs found for group: " << groupName << "\n";
        return {};
    }

    Graph subgraph;
    std::unordered_set<Node> relevant_nodes;

    // First pass: Find colocalization edges involving the MGE group and collect all connected nodes.
    for (const Edge& edge : g.edges) {
        if (!edge.isColo) continue;

        bool sourceInGroup = !edge.source.isARG && mgeIDsInGroup.count(edge.source.id);
        bool targetInGroup = !edge.target.isARG && mgeIDsInGroup.count(edge.target.id);

        if (sourceInGroup || targetInGroup) {
            subgraph.edges.insert(edge);
            relevant_nodes.insert(edge.source);
            relevant_nodes.insert(edge.target);
        }
    }

    // Second pass: Add temporal edges that connect any of the nodes we've collected.
    for (const Edge& edge : g.edges) {
        if (edge.isColo) continue;

        if (relevant_nodes.count(edge.source) && relevant_nodes.count(edge.target)) {
            subgraph.edges.insert(edge);
        }
    }

    // Add all collected nodes to the subgraph's node set.
    for(const auto& node : relevant_nodes) {
        subgraph.nodes.insert(node);
    }

    return subgraph;
}


Graph filterGraphByDisease(const Graph& g, const std::string& disease, const std::map<int, std::string>& patientToDiseaseMap) {
    // Step 1: Find all patient IDs that match the target disease.
    std::set<int> targetPatientIDs;
    for (const auto& [patientID, diseaseName] : patientToDiseaseMap) {
        if (diseaseName == disease) {
            targetPatientIDs.insert(patientID);
        }
    }

    if (targetPatientIDs.empty()) {
        std::cerr << "Warning: No patients found for disease '" << disease << "'. The resulting graph will be empty.\n";
        return {};
    }

    // Step 2: Identify all nodes that belong to those patients by checking colocalization edges.
    std::unordered_set<Node> relevant_nodes;
    for (const auto& edge : g.edges) {
        if (!edge.isColo) continue;
        for (int patientID : edge.individuals) {
            if (targetPatientIDs.count(patientID)) {
                relevant_nodes.insert(edge.source);
                relevant_nodes.insert(edge.target);
                break; 
            }
        }
    }
    
    // Step 3: Build the subgraph using only the relevant nodes.
    Graph subgraph;
    subgraph.nodes = std::set<Node>(relevant_nodes.begin(), relevant_nodes.end());

    for (const auto& edge : g.edges) {
        // Add any edge (colocalization or temporal) if BOTH its nodes are in our relevant set.
        if (relevant_nodes.count(edge.source) && relevant_nodes.count(edge.target)) {
            subgraph.edges.insert(edge);
        }
    }
    
    return subgraph;
}