// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

/* Traverse the graph based on earliest colocalizations */
#include "../include/traversal.h"
#include "../include/graph.h"
#include "../include/Timepoint.h"
#include "../include/id_maps.h"
#include <map>
#include <tuple>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <queue>

// This function traverses the graph and collects colocalization timepoints for each ARG and MGE pair.
// It uses an adjacency map to represent the graph structure, allowing for efficient traversal of nodes and their neighbors based on the defined edges.
void traverseAdjacency(const Graph& graph, const std::unordered_map<Node, std::unordered_set<Node>>& adjacency, 
    std::map<std::pair<int, int>, std::multiset<Timepoint>>& colocalizationTimeline) {
    for (const auto& node : graph.nodes) {
        // Skip if node is not an ARG
        if (!node.isARG) continue;

        for (const auto& neighbor : adjacency.at(node)) {
            // Skip if neighbor is not an MGE
            if (neighbor.isARG) continue;

            // ARG to MGE, record colocalization
            int argId = node.id;
            int mgeId = neighbor.id;
            Timepoint tp = node.timepoint;

            colocalizationTimeline[{argId, mgeId}].insert(tp);
        }
    }
}


// The function builds a timeline of colocalizations for each individual, ARG, and MGE pair, allowing for further analysis of colocalization patterns over time.
// This allows for efficient tracking of colocalization events for each individual across different ARG and MGE pairs.
void traverseGraph(const Graph& graph, 
    std::map<std::tuple<int, int, int>, std::set<Timepoint>>& colocalizationByIndividual) {
    // std::unordered_map<Node, std::unordered_set<Node>> adjacency;
    for (const auto& edge : graph.edges) {
        if (!edge.isColo) continue;
        int arg_id = edge.source.isARG ? edge.source.id : edge.target.id;
        int mge_id = edge.source.isARG ? edge.target.id : edge.source.id;
        Timepoint tp = edge.source.timepoint; // or target.timepoint

        for (int ind_id : edge.individuals) {
            auto key = std::make_tuple(ind_id, arg_id, mge_id);
            colocalizationByIndividual[key].insert(tp);
        }
    }
}



/********************************************  Traverse by Time 1 (not considering patients) ******************************************/
// This version considers only the first occurrence of colocalization between ARG and MGE
// It builds a timeline of colocalizations for each ARG and MGE pair, without individual
void findFirstOccurrence(const Graph& graph, std::unordered_map<Node, std::unordered_set<Node>>& adjacency,
                     std::map<std::pair<int, int>, Node>& firstOccurrence){
    for (const auto& edge : graph.edges) {
        if (!edge.isColo) continue;

        int arg = edge.source.isARG ? edge.source.id : edge.target.id;
        int mge = edge.source.isARG ? edge.target.id : edge.source.id;
        Node src = edge.source.isARG ? edge.source : edge.target;

        auto key = std::make_pair(arg, mge);

        if (!firstOccurrence.count(key) || src.timepoint < firstOccurrence[key].timepoint)
            firstOccurrence[key] = src;
    }
    // std::cout << "First occurrences: " << firstOccurrence.size() << "\n";
}


// This function performs a BFS-like traversal starting from the first occurrence node
// It explores the graph in a forward-in-time manner, collecting colocalization timepoints for the ARG and MGE pair.
void bfsTemporal(const Node& start, const std::unordered_map<Node, std::unordered_set<Node>>& adjacency, std::map<std::pair<int, int>, std::set<Timepoint>>& colocalizationTimeline){
    std::queue<Node> q;
    std::unordered_set<Node> visited;

    q.push(start);
    visited.insert(start);

    while (!q.empty()) {
        Node curr = q.front();
        q.pop();

        for (const Node& neighbor : adjacency.at(curr)) {
            if (neighbor.timepoint < curr.timepoint) continue; // Enforce forward-in-time
            if (visited.count(neighbor)) continue;

            visited.insert(neighbor);
            q.push(neighbor);

            if (curr.isARG != neighbor.isARG) {
                int arg = curr.isARG ? curr.id : neighbor.id;
                int mge = curr.isARG ? neighbor.id : curr.id;

                colocalizationTimeline[{arg, mge}].insert(neighbor.timepoint);
            }
        }
    }
}

// This function builds a timeline of colocalizations for each ARG and MGE pair
// It first finds the first occurrence of colocalization between ARG and MGE, then performs a BFS traversal starting from each first occurrence node to collect all colocalization timepoints.
void traverseTempGraph(const Graph& graph, std::unordered_map<Node, std::unordered_set<Node>>& adjacency,
                     std::map<std::pair<int, int>, Node>& firstOccurrence, std::map<std::pair<int, int>, std::set<Timepoint>>& colocalizationsByTime){
    findFirstOccurrence(graph, adjacency, firstOccurrence);
    for (const auto& [key, startNode] : firstOccurrence) {
        bfsTemporal(startNode, adjacency, colocalizationsByTime);
    }
}



/************************************************  Traverse by Time 2  **********************************************/
// This version considers individuals and their first occurrence of colocalization
// It builds a timeline of colocalizations for each individual, ARG, and MGE pair
void findFirstOccurrenceByInd(
    const Graph& graph,
    std::unordered_map<Node, std::unordered_set<Node>>& adjacency,
    std::map<std::tuple<int, int, int>, Node>& firstOccurrenceByInd // This will store the first occurrence of the ARG-MGE pair for the individual (e.g{[key], srcNode(ARG)} )
) {
    for (const auto& edge : graph.edges) {
        if (!edge.isColo) continue;
        if (edge.source.isARG == edge.target.isARG) {
            std::cerr << "Warning: Invalid ARG-MGE edge with same type nodes (ID " << edge.source.id << ", " << edge.target.id << ")\n";
            continue;
        }
        
        int arg = edge.source.isARG ? edge.source.id : edge.target.id;
        int mge = edge.source.isARG ? edge.target.id : edge.source.id;
        Node src = edge.source.isARG ? edge.source : edge.target;

        for (int ind : edge.individuals) {
            auto key = std::make_tuple(ind, arg, mge);

            if (!firstOccurrenceByInd.count(key) || static_cast<int>(src.timepoint) < static_cast<int>(firstOccurrenceByInd[key].timepoint)) {
                firstOccurrenceByInd[key] = src;
            }
        }
    }
    std::cout << "First occurrences by individual found: " << firstOccurrenceByInd.size() << "\n";
}



// This function performs a BFS-like traversal starting from the first occurrence node
// It explores the graph in a forward-in-time manner, collecting colocalization timepoints for the specified individual, ARG, and MGE pair.
void temporalTimelineTraversal(
    const Node& start,
    const std::unordered_map<Node, std::unordered_set<Node>>& adjacency,
    const std::map<std::pair<Node, Node>, std::vector<Edge>>& edgeMap,
    int ind, int arg, int mge,
    std::map<std::tuple<int, int, int>, std::set<Timepoint>>& colocalizationTimelineByInd
) {
    std::queue<Node> q;
    std::unordered_set<Node> visited;

    q.push(start);
    visited.insert(start);

    auto key = std::make_tuple(ind, arg, mge);
    colocalizationTimelineByInd[key].insert(start.timepoint);

    while (!q.empty()) {
        Node curr = q.front();
        q.pop();

        for (const Node& neighbor : adjacency.at(curr)) {
            if (neighbor.timepoint < curr.timepoint) continue;
            if (visited.count(neighbor)) continue;

            auto it = edgeMap.find({curr, neighbor});
            if (it == edgeMap.end()) continue;

            for (const Edge& edge : it->second) {
                if (!edge.isColo || !edge.individuals.count(ind)) continue;

                int this_arg = curr.isARG ? curr.id : neighbor.id;
                int this_mge = curr.isARG ? neighbor.id : curr.id;

                if (this_arg == arg && this_mge == mge) {
                    colocalizationTimelineByInd[key].insert(neighbor.timepoint);
                    visited.insert(neighbor);
                    q.push(neighbor);
                    break; // one valid edge is enough
                }
            }
        }
    }
}

// Builds an edge map for efficient lookup of edges between nodes.
// The adjacency map is used to represent the graph structure, allowing for efficient traversal of nodes and their neighbors based on the defined edges.
// The function iterates over the first occurrences of each ARG-MGE pair and performs a BFS traversal starting from each first occurrence node. During the traversal, it collects colocalization timepoints for the specified individual, ARG, and MGE pair, ensuring that only valid colocalizations are recorded.
void traverseGraphByInd(const Graph& graph, std::unordered_map<Node, std::unordered_set<Node>>& adjacency, const std::set<Edge>& edges,
                     std::map<std::tuple<int, int, int>, Node>& firstOccurrenceByInd, std::map<std::tuple<int, int, int>, std::set<Timepoint>>& colocalizationTimelineByInd) {
    findFirstOccurrenceByInd(graph, adjacency, firstOccurrenceByInd);
    std::map<std::pair<Node, Node>, std::vector<Edge>> edgeMap;
    
    for (const auto& edge : graph.edges) {
        edgeMap[{edge.source, edge.target}].push_back(edge);
        edgeMap[{edge.target, edge.source}].push_back(edge); // undirected lookup
    }

    for (const auto& [key, startNode] : firstOccurrenceByInd) {
        auto [ind, arg, mge] = key;
        temporalTimelineTraversal(startNode, adjacency, edgeMap, ind, arg, mge, colocalizationTimelineByInd);
    }
}




/************************************************   **********************************************/
// This function retrieves the top K entities (ARGs or MGEs) based on their frequency of occurrence in the graph
std::vector<std::pair<int, int>> getTopKEntities(const Graph& graph, bool isARG, unsigned int K) {
    std::unordered_map<int, int> countMap;

    for (const Edge& edge : graph.edges) {
        if (!edge.isColo) continue;
        const Node& node = isARG ? edge.source : edge.target;
 
        if (node.isARG == isARG && node.timepoint != Timepoint::Donor) {
            for (size_t i = 0; i < edge.individuals.size(); ++i) {
                // Count occurrences of ARG or MGE by individual
                countMap[node.id]++;
            }
        }
    }

    std::vector<std::pair<int, int>> freqList(countMap.begin(), countMap.end());
    std::sort(freqList.begin(), freqList.end(), [](auto& a, auto& b) {
        return a.second > b.second;
    });

    if (K < freqList.size())
        freqList.resize(K);
    
    return freqList;
}

// Retrieve the timeline of a specific ARG and all MGEs it colocalizes with over time
void getTimelineForARG(const Graph& graph, const std::string& argName) {
    std::map<int, std::set<Timepoint>> timeline;
    int count = 0;
    int argID = getARGId(argName);
    if (argID == -1) {
        std::cerr << "ARG with name '" << argName << "' not found.\n";
        return;
    }
    for (const Edge& edge: graph.edges){
        if (!edge.isColo) continue;

        if (edge.source.id == argID || edge.target.id == argID) {
            int mgeID = edge.source.id == argID ? edge.target.id : edge.source.id;
            for (size_t i = 0; i < edge.individuals.size(); ++i) {
                timeline[mgeID].insert(edge.source.timepoint);
            }
            
        }
    }
    
    std::cout << "Timeline for ARG ID " << argID << ":\n";
    for (const auto& [mgeID, timepoints] : timeline) {
        std::cout << "  MGE ID: " << getMGEName(mgeID) << ", Timepoints: ";
        for (const auto& tp : timepoints) {
            count++;
            std::cout << tp << " ";
        }
        std::cout << "\n";
    }
    std::cout << "Total MGEs colocalized with ARG " << argName << ": " << timeline.size() << "\n";
    std::cout << "Total number of Timepoints for ARG " << argName << ": " << count << "\n";
}

// Retrieve the timeline of specific MGE and all ARGs it colocalizes with over time
void getTimelineForMGE(const Graph& graph, const std::string& mgeName) {
    std::map<int, std::set<Timepoint>> timeline;
    int count = 0;
    int mgeID = getMGEId(mgeName);
    if (mgeID == -1) {
        std::cerr << "MGE with name '" << mgeName << "' not found.\n";
        return;
    }
    for (const Edge& edge: graph.edges){
        if (!edge.isColo) continue;
        if (edge.source.id == mgeID || edge.target.id == mgeID) {
            int argID = edge.source.id == mgeID ? edge.target.id : edge.source.id;
            for (size_t i = 0; i < edge.individuals.size(); ++i) {
                // Insert the timepoint for the ARG-MGE colocalization
                timeline[argID].insert(edge.source.timepoint);
            }
        }
    }
    std::cout << "Timeline for MGE " << mgeName << ":\n";
    for (const auto& [argID, timepoints] : timeline) {
        std::cout << "  ARG ID: " << getARGName(argID) << ", Timepoints: ";
        for (const auto& tp : timepoints) {
            count++;
            std::cout << tp << " ";
        }
        std::cout << "\n";
    }
    std::cout << "Total ARGs colocalized with MGE " << mgeName << ": " << timeline.size() << "\n";
    std::cout << "Total number of Timepoints for MGE " << mgeName << ": " << count << "\n";
}


std::map<Timepoint, int> computeNodeDegreeOverTime(const Graph& graph, bool isARG, const std::string& name) {
    std::map<Timepoint, int> degreeOverTime;

    int nodeID = isARG ? getARGId(name) : getMGEId(name);
    if (nodeID == -1) {
        std::cerr << "Node with name '" << name << "' not found.\n";
        return degreeOverTime;
    }

    for (const Edge& e : graph.edges) {
        if (!e.isColo) continue; // Only count colocalization edges

        if (e.source.id == nodeID && e.source.isARG == isARG) {
            degreeOverTime[e.source.timepoint]++;
        } else if (e.target.id == nodeID && e.target.isARG == isARG) {
            degreeOverTime[e.target.timepoint]++;
        }
    }

    return degreeOverTime;
}
