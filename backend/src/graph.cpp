// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

#include <iostream>
#include "../include/graph.h"

void buildAdjacency(const Graph& g, std::unordered_map<Node, std::unordered_set<Node>>& adj) {
    // std::unordered_map<int, std::vector<Node>> nodesById;

    // for (const Node& node : g.nodes) {
    //     nodesById[node.id].push_back(node);
    // }
    for (const Edge& edge : g.edges) {
        const Node& src = edge.source;
        const Node& tgt = edge.target;

        if (edge.isColo) {
            // Bidirectional colocalization edge
            adj[src].insert(tgt);
            adj[tgt].insert(src);
        } else {
            // Directional temporal edge (already ensured during edge creation)
            adj[src].insert(tgt);
            // // Directional temporal edge
            // if (static_cast<int>(src.timepoint) < static_cast<int>(tgt.timepoint)) {
            //     adj[src].insert(tgt); // temporal forward
            // }else if (static_cast<int>(src.timepoint) > static_cast<int>(tgt.timepoint)) {
            //     adj[tgt].insert(src);
            // }
        }
    }
}
