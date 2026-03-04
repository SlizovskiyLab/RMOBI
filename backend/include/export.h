// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

#ifndef EXPORT_H
#define EXPORT_H

#include <string>
#include "graph.h"

std::string getNodeName(const Node& node);
std::string getNodeLabel(const Node& node);
// std::string getTimepointColor(const Timepoint& tp);
bool isTemporalEdge(const Edge& edge);  
std::string getMGEGroupShape(const std::string& groupName);
void exportToDot(const Graph& g, const std::string& filename, bool showLabels = true);
void exportParentTemporalGraphDot(const Graph& g, const std::string& filename, bool showLabels=true);
#endif