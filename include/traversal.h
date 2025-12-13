// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

#ifndef TRAVERSAL_H
#define TRAVERSAL_H

#include <map>
#include <unordered_map>
#include <tuple>
#include <set>
#include "Timepoint.h"
#include "graph.h"

// std::map<std::pair<int, int>, std::set<Timepoint>> colocalizationTimeline;
void traverseAdjacency(const Graph& graph, const std::unordered_map<Node, std::unordered_set<Node>>& adjacency, 
    std::map<std::pair<int, int>, std::multiset<Timepoint>>& colocalizationTimeline);

void traverseGraph(const Graph& graph, 
    std::map<std::tuple<int, int, int>, std::set<Timepoint>>& colocalizationByIndividual);

void findFirstOccurrence(const Graph& graph, std::unordered_map<Node, std::unordered_set<Node>>& adjacency,
                      std::map<std::pair<int, int>, Node>& firstOccurrence);

void bfsTemporal(const Node& startNode, const std::unordered_map<Node, std::unordered_set<Node>>& adjacency,
                 std::map<std::pair<int, int>, std::set<Timepoint>>& colocalizationsByTime);

void traverseTempGraph(const Graph& graph, std::unordered_map<Node, std::unordered_set<Node>>& adjacency,
                     std::map<std::pair<int, int>, Node>& firstOccurrence, std::map<std::pair<int, int>, std::set<Timepoint>>& colocalizationsByTime);


void findFirstOccurrenceByInd(const Graph& graph, std::unordered_map<Node, std::unordered_set<Node>>& adjacency,
    std::map<std::tuple<int, int, int>, Node>& firstOccurrenceByInd);
    
// void bfsTemporalByInd(const Node& start, const std::unordered_map<Node, std::unordered_set<Node>>& adjacency, const std::set<Edge>& edges, int ind,
//     std::map<std::tuple<int, int, int>, std::set<Timepoint>>& colocalizationTimelineByInd);

void traverseGraphByInd(const Graph& graph, std::unordered_map<Node, std::unordered_set<Node>>& adjacency, const std::set<Edge>& edges,
                     std::map<std::tuple<int, int, int>, Node>& firstOccurrenceByInd, std::map<std::tuple<int, int, int>, std::set<Timepoint>>& colocalizationTimelineByInd);

void temporalTimelineTraversal(
    const Node& start,
    const std::unordered_map<Node, std::unordered_set<Node>>& adjacency,
    const std::map<std::pair<Node, Node>, std::vector<Edge>>& edgeMap,
    int ind, int arg, int mge,
    std::map<std::tuple<int, int, int>, std::set<Timepoint>>& colocalizationTimelineByInd
);

std::vector<std::pair<int, int>> getTopKEntities(const Graph& graph, bool isARG, unsigned int K);

void getTimelineForARG(const Graph& graph, const std::string& argName);
void getTimelineForMGE(const Graph& graph, const std::string& mgeName);

std::map<Timepoint, int> computeNodeDegreeOverTime(const Graph& graph, bool isARG, const std::string& name);


inline std::ostream& operator<<(std::ostream& os, const Timepoint& tp) {
    switch (tp) {
        case Timepoint::Donor: return os << "donor";
        case Timepoint::PreFMT:   return os << "pre";
        case Timepoint::PostFMT_001:  return os << "post_001";
        case Timepoint::PostFMT_002:  return os << "post_002";
        case Timepoint::PostFMT_003:  return os << "post_003";
        case Timepoint::PostFMT_006:  return os << "post_006";
        case Timepoint::PostFMT_007:  return os << "post_007";
        case Timepoint::PostFMT_012:  return os << "post_012";
        case Timepoint::PostFMT_013:  return os << "post_013";              
        case Timepoint::PostFMT_014:  return os << "post_014";
        case Timepoint::PostFMT_015:  return os << "post_015";
        case Timepoint::PostFMT_016:  return os << "post_016";
        case Timepoint::PostFMT_020:  return os << "post_020";
        case Timepoint::PostFMT_021:  return os << "post_021";
        case Timepoint::PostFMT_028:  return os << "post_028";
        case Timepoint::PostFMT_029:  return os << "post_029";
        case Timepoint::PostFMT_030:  return os << "post_030";
        case Timepoint::PostFMT_031:  return os << "post_031";
        case Timepoint::PostFMT_035:  return os << "post_035";
        case Timepoint::PostFMT_036:  return os << "post_036";
        case Timepoint::PostFMT_040:  return os << "post_040";
        case Timepoint::PostFMT_041:  return os << "post_041";
        case Timepoint::PostFMT_042:  return os << "post_042";
        case Timepoint::PostFMT_044:  return os << "post_044";
        case Timepoint::PostFMT_054:  return os << "post_054";
        case Timepoint::PostFMT_056:  return os << "post_056";
        case Timepoint::PostFMT_059:  return os << "post_059";
        case Timepoint::PostFMT_061:  return os << "post_061";
        case Timepoint::PostFMT_063:  return os << "post_063";
        case Timepoint::PostFMT_064:  return os << "post_064";
        case Timepoint::PostFMT_065:  return os << "post_065";
        case Timepoint::PostFMT_068:  return os << "post_068";
        case Timepoint::PostFMT_081:  return os << "post_081";
        case Timepoint::PostFMT_084:  return os << "post_084";
        case Timepoint::PostFMT_090:  return os << "post_090";
        case Timepoint::PostFMT_094:  return os << "post_094";
        case Timepoint::PostFMT_095:  return os << "post_095";
        case Timepoint::PostFMT_097:  return os << "post_097";
        case Timepoint::PostFMT_098:  return os << "post_098";
        case Timepoint::PostFMT_111:  return os << "post_111";
        case Timepoint::PostFMT_112:  return os << "post_112";
        case Timepoint::PostFMT_120:  return os << "post_120";
        case Timepoint::PostFMT_135:  return os << "post_135";
        case Timepoint::PostFMT_140:  return os << "post_140";
        case Timepoint::PostFMT_150:  return os << "post_150";
        case Timepoint::PostFMT_179:  return os << "post_179";
        case Timepoint::PostFMT_180:  return os << "post_180";
        case Timepoint::PostFMT_195:  return os << "post_195";
        case Timepoint::PostFMT_365:  return os << "post_365";
        case Timepoint::PostFMT_384:  return os << "post_384";
        case Timepoint::PostFMT_408:  return os << "post_408";
        case Timepoint::PostFMT_730:  return os << "post_730";
        default:               return os << "unknown";
    }
}

#endif // End of traversal.h
