// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

#pragma once
#include <string>
#include <iostream>
#include <ostream>
// enum class Timepoint {
//     Donor = -1,
//     Pre = 0,
//     Day4 = 4,
//     Day7 = 7,
//     Day14 = 14
// };


enum class Timepoint {
    Donor = 1000,
    PreFMT = 0,
    PostFMT_001 = 1,
    PostFMT_002 = 2,
    PostFMT_003 = 3,
    PostFMT_006 = 6,
    PostFMT_007 = 7,
    PostFMT_012 = 12,
    PostFMT_013 = 13,
    PostFMT_014 = 14,
    PostFMT_015 = 15,
    PostFMT_016 = 16,
    PostFMT_020 = 20,
    PostFMT_021 = 21,
    PostFMT_028 = 28,
    PostFMT_029 = 29,
    PostFMT_030 = 30,
    PostFMT_031 = 31,
    PostFMT_035 = 35,
    PostFMT_036 = 36,
    PostFMT_040 = 40,
    PostFMT_041 = 41,
    PostFMT_042 = 42,
    PostFMT_044 = 44,
    PostFMT_054 = 54,
    PostFMT_056 = 56,
    PostFMT_059 = 59,
    PostFMT_061 = 61,
    PostFMT_063 = 63,
    PostFMT_064 = 64,
    PostFMT_065 = 65,
    PostFMT_068 = 68,
    PostFMT_081 = 81,
    PostFMT_084 = 84,
    PostFMT_090 = 90,
    PostFMT_094 = 94,
    PostFMT_095 = 95,
    PostFMT_097 = 97,
    PostFMT_098 = 98,
    PostFMT_111 = 111,
    PostFMT_112 = 112,
    PostFMT_120 = 120,
    PostFMT_135 = 135,
    PostFMT_140 = 140,
    PostFMT_150 = 150,
    PostFMT_179 = 179,
    PostFMT_180 = 180,
    PostFMT_195 = 195,
    PostFMT_365 = 365,
    PostFMT_384 = 384,
    PostFMT_408 = 408,
    PostFMT_730 = 730,
};

inline std::string toString(Timepoint tp) {
    switch (tp) {
        case Timepoint::Donor: return "donor";
        case Timepoint::PreFMT: return "pre";
        case Timepoint::PostFMT_001: return "post_001";
        case Timepoint::PostFMT_002: return "post_002";
        case Timepoint::PostFMT_003: return "post_003";
        case Timepoint::PostFMT_006: return "post_006";
        case Timepoint::PostFMT_007: return "post_007";
        case Timepoint::PostFMT_012: return "post_012";
        case Timepoint::PostFMT_013: return "post_013";
        case Timepoint::PostFMT_014: return "post_014";
        case Timepoint::PostFMT_015: return "post_015";
        case Timepoint::PostFMT_016: return "post_016";
        case Timepoint::PostFMT_020: return "post_020";
        case Timepoint::PostFMT_021: return "post_021";
        case Timepoint::PostFMT_028: return "post_028";
        case Timepoint::PostFMT_029: return "post_029";
        case Timepoint::PostFMT_030: return "post_030";
        case Timepoint::PostFMT_031: return "post_031";
        case Timepoint::PostFMT_035: return "post_035";
        case Timepoint::PostFMT_036: return "post_036";
        case Timepoint::PostFMT_040: return "post_040";
        case Timepoint::PostFMT_041: return "post_041";
        case Timepoint::PostFMT_042: return "post_042";
        case Timepoint::PostFMT_044: return "post_044";
        case Timepoint::PostFMT_054: return "post_054";
        case Timepoint::PostFMT_056: return "post_056";
        case Timepoint::PostFMT_059: return "post_059";
        case Timepoint::PostFMT_061: return "post_061";
        case Timepoint::PostFMT_063: return "post_063";
        case Timepoint::PostFMT_064: return "post_064";
        case Timepoint::PostFMT_065: return "post_065";
        case Timepoint::PostFMT_068: return "post_068";
        case Timepoint::PostFMT_081: return "post_081";
        case Timepoint::PostFMT_084: return "post_084";
        case Timepoint::PostFMT_090: return "post_090";
        case Timepoint::PostFMT_094: return "post_094";
        case Timepoint::PostFMT_095: return "post_095";
        case Timepoint::PostFMT_097: return "post_097";
        case Timepoint::PostFMT_098: return "post_098";
        case Timepoint::PostFMT_111: return "post_111";
        case Timepoint::PostFMT_112: return "post_112";
        case Timepoint::PostFMT_120: return "post_120";
        case Timepoint::PostFMT_135: return "post_135";
        case Timepoint::PostFMT_140: return "post_140";
        case Timepoint::PostFMT_150: return "post_150";
        case Timepoint::PostFMT_179: return "post_179";
        case Timepoint::PostFMT_180: return "post_180";
        case Timepoint::PostFMT_195: return "post_195";
        case Timepoint::PostFMT_365: return "post_365";
        case Timepoint::PostFMT_384: return "post_384";
        case Timepoint::PostFMT_408: return "post_408";
        case Timepoint::PostFMT_730: return "post_730";
        default: return "unknown";
    }
}