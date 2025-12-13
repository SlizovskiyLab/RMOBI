// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

#include <fstream>
#include <string>
#include <unordered_set>
#include <set>
#include <map>
#include <utility>
#include <algorithm>
#include "graph.h"
#include "id_maps.h"
#include "export.h"

std::string getNodeName(const Node& node) {
    std::string name = "N_" + std::to_string(node.id);
    name += node.isARG ? "_ARG_" : "_MGE_";
    name += std::to_string(static_cast<int>(node.timepoint));
    return name;
}

std::string getNodeLabel(const Node& node) {
    std::string label = node.isARG ? getARGName(node.id) : getMGENameForLabel(node.id);
    label += "\\n" + toString(node.timepoint);
    return label;
}

std::string getTimepointColor(const Timepoint& tp) {
    int timeValue = static_cast<int>(tp);
    if (timeValue == 1000) return "\"yellow\"";
    if (timeValue == 0)  return "\"red\"";
    if (timeValue > 0 && timeValue < 31) return "\"#99D2FF\""; 
    if (timeValue > 30 && timeValue < 61) return "\"#4D9DFF\""; 
    if (timeValue > 60) return "\"#3A6EFF\"";
    return "\"green\""; // Fallback for any unexpected values
}

bool isTemporalEdge(const Edge& edge) {
    return edge.source.id == edge.target.id &&
           edge.source.isARG == edge.target.isARG &&
           edge.source.timepoint != edge.target.timepoint;
}

// auto timepointOrder = [](Timepoint tp) -> int {
//     if (tp == Timepoint::Donor) return -1;       // Donor comes first
//     if (tp == Timepoint::PreFMT) return 0;       // PreFMT next
//     return static_cast<int>(tp);               
// };


std::string getMGEGroupShape(const std::string& groupName) {
    if (groupName == "plasmid" || groupName == "Colicin_plasmid" || groupName == "Inc_plasmid") {
        return "diamond";
    }
    if (groupName == "prophage") {
        return "hexagon";
    }
    if (groupName == "virus") {
        return "triangle";
    }
    if (groupName == "ICE" || groupName == "ICEberg") {
        return "octagon";
    }
    if (groupName == "replicon") {
        return "parallelogram";
    }
    if (groupName == "likely IS/TE") {
        return "trapezium";
    }
    return "box"; // Default for UNCLASSIFIED or others
}

// void exportToDot(const Graph& g, const std::string& filename, bool showLabels) {
//     std::ofstream file(filename);
//     file << "digraph G {\n";
//     file << "  layout=sfdp;\n";
//     file << "  graph [nodesep=2.0, ranksep=2.0, overlap=false];\n";
//     file << "  node [style=filled];\n";

//     std::unordered_set<Node> active_nodes;
//     for (const Edge& edge : g.edges) {
//         active_nodes.insert(edge.source);
//         active_nodes.insert(edge.target);
//     }

//     for (const Node& node : active_nodes) {
//         std::string nodeName = getNodeName(node);
//         std::string label = showLabels ? getNodeLabel(node) : "";
//         std::string color = getTimepointColor(node.timepoint);
        
//         std::string shape;
//         if (node.isARG) {
//             shape = "circle";
//         } else {
//             std::string groupName = getMGEGroupName(node.id);
//             shape = getMGEGroupShape(groupName);
//         }

//         file << "  " << nodeName
//              << " [label=\"" << label << "\", shape=" << shape
//              << ", fixedsize=true, width=0.5, height=0.5, fillcolor=" << color << "]\n";
//     }

//     std::set<std::pair<Node, Node>> processedColoEdges;

//     for (const Edge& edge : g.edges) {
//         std::string sourceName = getNodeName(edge.source);
//         std::string targetName = getNodeName(edge.target);

//         std::string color;
//         std::string style;
//         std::string extraAttributes;
//         double penwidth = 4.0; 
//         std::string edgeLabel = "";

//         if (edge.isColo) {
//             auto canonical_pair = std::minmax(edge.source, edge.target);
//             if (processedColoEdges.count(canonical_pair)) {
//                 continue;
//             }
//             processedColoEdges.insert(canonical_pair);

//             color = "\"#696969\"";
//             style = "solid";
//             double count = edge.individuals.size();
//             if (count > 1) {
//                 penwidth = 4.0 + (count - 1) * 2.0;
//             }
//             penwidth = std::min(10.0, penwidth);
//             if (showLabels) {
//                 edgeLabel = std::to_string((int)count);
//             }
//             extraAttributes = "dir=both";

//         } else if (isTemporalEdge(edge)) {
//             style = "dashed";
//             double count = edge.weight;
//             if (count > 1) {
//                 penwidth = 4.0 + (count - 1) * 2.0;
//             }
//             penwidth = std::min(10.0, penwidth);
//             if (showLabels) {
//                 edgeLabel = std::to_string((int)count);
//             }
            
//             Timepoint src_tp = edge.source.timepoint;
//             Timepoint tgt_tp = edge.target.timepoint;

//             bool src_is_post = (src_tp != Timepoint::Donor && src_tp != Timepoint::PreFMT);
//             bool tgt_is_post = (tgt_tp != Timepoint::Donor && tgt_tp != Timepoint::PreFMT);

//             if (src_tp == Timepoint::Donor && tgt_tp == Timepoint::PreFMT) {
//                 color =  "\"#006400\"";
//             } else if (src_tp == Timepoint::Donor && tgt_is_post) {
//                 color = "\"#4B0082\"";
//             } else if (src_tp == Timepoint::PreFMT && tgt_is_post) {
//                 color = "\"orange\"";
//             } else if (src_is_post && tgt_is_post) {
//                 color = "\"black\"";
//             } else {
//                 color = "\"black\"";
//             }
            
//         } else {
//             color = "\"#808080\"";
//             style = "solid";
//         }

//         file << "  " << sourceName << " -> " << targetName
//              << " [style=" << style
//              << ", color=" << color
//              << ", arrowsize=0.3"
//              << ", penwidth=" << penwidth
//              << ", label=\"" << edgeLabel << "\"";

//         if (!extraAttributes.empty()) {
//             file << ", " << extraAttributes;
//         }

//         file << "]\n";
//     }

//     file << "}\n";
// }


// void exportParentTemporalGraphDot(const Graph& g, const std::string& filename, bool showLabels) {
//     std::ofstream file(filename);
//     file << "digraph ParentTemporalGraph {\n";
//     file << "  rankdir=LR;\n"; // Left to Right layout
//     file << "  graph [nodesep=2.0, ranksep=2.0, overlap=false];\n";
//     file << "  edge [splines=true, minlen=2, arrowsize=0.6, penwidth=2];\n";
//     if (showLabels) {
//         file << "  label=\"Parent Temporal Graph\";\n";
//     } else {
//         file << "  label=\"\";\n"; // No label if showLabels is false
//     }
//     file << "  layout=sfdp;\n";
//     file << "  node [style=filled];\n";

//     struct ParentNodeInfo {
//         std::string name;
//         Timepoint tp;
//     };

//     int colocCounter = 0;
//     std::map<std::pair<int,int>, std::vector<ParentNodeInfo>> colocMap;
//     std::map<std::tuple<int,int,Timepoint>, std::string> uniqueParents;

//     // --- Create parent nodes ---
//     for (const Edge& edge : g.edges) {
//         if (edge.isColo) {
//             const Node& argNode = edge.source.isARG ? edge.source : edge.target;
//             const Node& mgeNode = edge.source.isARG ? edge.target : edge.source;

//             int argId = argNode.id;
//             int mgeId = mgeNode.id;
//             Timepoint tp = argNode.timepoint;

//             auto uniqueKey = std::make_tuple(argId, mgeId, tp);

//             // Avoid duplicate parent nodes
//             if (!uniqueParents.count(uniqueKey)) {
//                 std::string colocName = "Parent_" + std::to_string(++colocCounter);
//                 uniqueParents[uniqueKey] = colocName;

//                 // Color and shape
//                 std::string color = getTimepointColor(argNode.timepoint);
//                 std::string groupName = getMGEGroupName(mgeId);
//                 std::string shape = getMGEGroupShape(groupName);
//                 std::string label = getARGName(argId) + "\n" + getMGENameForLabel(mgeId) 
//                                     + "\n" + toString(tp);

//                 // Write node
//                 file << "  \"" << colocName << "\" [label=\"" << label
//                      << "\", shape=" << shape
//                      << ", fixedsize=true, width=0.6, height=0.6, fillcolor="
//                      << color << "];\n";
//             }

//             auto pairKey = std::make_pair(argId, mgeId);
//             colocMap[pairKey].push_back({uniqueParents[uniqueKey], tp});
//         }
//     }

//     // --- Draw temporal edges for each ARG–MGE pair ---
//     for (auto& entry : colocMap) {
//     auto& parentNodes = entry.second;

//     // Custom sort: Donor → PreFMT → PostFMT
//     std::sort(parentNodes.begin(), parentNodes.end(), [&](const ParentNodeInfo& a, const ParentNodeInfo& b) {
//         return timepointOrder(a.tp) < timepointOrder(b.tp);
//     });

//     // Connect consecutive temporal edges
//     for (size_t i = 0; i + 1 < parentNodes.size(); ++i) {
//         if (parentNodes[i].tp == parentNodes[i + 1].tp ||
//             parentNodes[i].name == parentNodes[i + 1].name) {
//             continue; // skip self or same timepoint
//         }

//         Timepoint src_tp = parentNodes[i].tp;
//         Timepoint tgt_tp = parentNodes[i + 1].tp;

//         bool src_is_post = (src_tp != Timepoint::Donor && src_tp != Timepoint::PreFMT);
//         bool tgt_is_post = (tgt_tp != Timepoint::Donor && tgt_tp != Timepoint::PreFMT);

//         std::string color;
//         if (src_tp == Timepoint::Donor && tgt_tp == Timepoint::PreFMT) {
//             color = "\"#006400\"";
//         } else if (src_tp == Timepoint::Donor && tgt_is_post) {
//             color = "\"#4B0082\"";
//         } else if (src_tp == Timepoint::PreFMT && tgt_is_post) {
//             color = "\"orange\"";
//         } else if (src_is_post && tgt_is_post) {
//             color = "\"black\"";
//         } else {
//             color = "\"black\"";
//         }

//         file << "  \"" << parentNodes[i].name << "\" -> \"" << parentNodes[i + 1].name
//              << "\" [style=dashed, color=" << color
//              << ", penwidth=5.0];\n";
//     }
// }
// file << "}\n";
// }

