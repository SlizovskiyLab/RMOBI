// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

#ifndef ID_MAPS_H
#define ID_MAPS_H

#include <unordered_map>
#include <string>

extern const std::unordered_map<int, std::string> argIdMap;
extern const std::unordered_map<int, std::string> mgeIdMap;
extern const std::unordered_map<int, std::string> mgeNameMap;
extern const std::unordered_map<int, std::string> argGroupMap;
extern const std::unordered_map<int, bool> argIDSNPConfirmation;
extern const std::unordered_map<int, std::string> argResistanceMap;
extern const std::unordered_map<int, std::string> mgeGroupMap;

std::string getARGName(int id);
std::string getMGEName(int id);
std::string getMGENameForLabel(int id);
std::string getARGGroupName(int id);
std::string getARGResistance(const std::string& name);
std::string getMGEGroupName(int id);

int getARGId(const std::string& name);
int getMGEId(const std::string& name);
int getMGEIdByName(const std::string& name);
int getARGGroupId(const std::string& name);

#endif
