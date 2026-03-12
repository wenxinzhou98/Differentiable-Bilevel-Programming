#pragma once

#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <string>
#include <utility>
#include <functional>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>

// Link: directed edge (from_node, to_node)
using Link = std::pair<int, int>;

// OD: origin-destination pair
using OD = std::pair<int, int>;

// Path: sequence of node IDs
using Path = std::vector<int>;

// Maps keyed by Link
using LinkMapD = std::map<Link, double>;

// Hash for Path (to use as key in unordered_map)
struct PathHash {
    std::size_t operator()(const Path& p) const {
        std::size_t seed = p.size();
        for (auto& v : p) {
            seed ^= std::hash<int>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

// Path flow: OD -> {path -> flow}
using PathFlowMap = std::map<OD, std::unordered_map<Path, double, PathHash>>;

// Path set: OD -> set of path strings
using PathSetMap = std::map<OD, std::set<std::string>>;

// Adjacency list: node -> list of neighbor nodes
using AdjList = std::map<int, std::vector<int>>;

// OD tree: origin -> list of destinations
using OdTree = std::map<int, std::vector<int>>;

// Demand map: OD -> demand value
using DemandMap = std::map<OD, double>;

// Convert path to string representation (for pathset)
inline std::string path_to_string(const Path& p) {
    std::string s;
    for (size_t i = 0; i < p.size(); ++i) {
        if (i > 0) s += "-";
        s += std::to_string(p[i]);
    }
    return s;
}
