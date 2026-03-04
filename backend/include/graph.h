// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

#ifndef GRAPH_H
#define GRAPH_H

#pragma once

#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <vector>
#include <tuple>
#include <ostream>
#include "Timepoint.h"

// ------------------ Node ------------------

struct Node {
    int id;
    bool isARG;
    Timepoint timepoint;
    bool requiresSNPConfirmation;

    bool operator<(const Node& other) const {
        // First, compare by the entity's ID and type
        if (id != other.id) return id < other.id;
        if (isARG != other.isARG) return isARG < other.isARG;

        // If they are the same entity, use custom chronological sorting for timepoints
        Timepoint tp1 = timepoint;
        Timepoint tp2 = other.timepoint;

        // If one is Donor and the other isn't, Donor always comes first.
        if (tp1 == Timepoint::Donor && tp2 != Timepoint::Donor) return true;
        if (tp1 != Timepoint::Donor && tp2 == Timepoint::Donor) return false;

        // Otherwise, sort by the standard numerical order (PreFMT < PostFMT).
        return static_cast<int>(tp1) < static_cast<int>(tp2);
    }

    bool operator==(const Node& other) const {
        return id == other.id && isARG == other.isARG && timepoint == other.timepoint;
    }

    friend std::ostream& operator<<(std::ostream& os, const Node& node) {
        // Print relevant Node info, e.g.:
        os << "Node(" << node.id << ")";
        return os;
    }
};

// Hash function for Node to use in unordered_map/set
namespace std {
    template <>
    struct hash<Node> {
        std::size_t operator()(const Node& n) const {
            return hash<int>()(n.id) ^ (hash<bool>()(n.isARG) << 1) ^ (hash<int>()(static_cast<int>(n.timepoint)) << 2);
        }
    };
}

// ------------------ Edge ------------------

struct Edge {
    Node source;
    Node target;
    bool isColo;
    std::set<int> individuals;
    int weight = 0;
    

    bool operator<(const Edge& other) const {
        return std::tie(source, target, isColo, weight) < std::tie(other.source, other.target, other.isColo, other.weight);
    }

    bool operator==(const Edge& other) const {
        return source == other.source &&
               target == other.target &&
               isColo == other.isColo &&
               weight == other.weight;
    }

};

// ------------------ Graph ------------------

struct Graph {
    std::set<Node> nodes;
    std::set<Edge> edges;
};

// Function declaration
void buildAdjacency(const Graph& g, std::unordered_map<Node, std::unordered_set<Node>>& adj);


#endif // GRAPH_H
