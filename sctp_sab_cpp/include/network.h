#pragma once

#include "types.h"
#include <Eigen/Sparse>

class Network {
public:
    std::string name;
    AdjList Innode;
    AdjList Outnode;
    std::vector<Link> Link_list;
    OdTree Odtree;

    // UE/SO solution state
    PathFlowMap pathflow;
    LinkMapD flow;
    int path_number = 0;
    PathSetMap pathset;
    DemandMap demand_map;

    // Path enumeration data
    std::map<OD, std::map<Path, int, std::less<Path>>> path_id;
    std::map<OD, int> od2num;
    int numberoflink = 0;
    int numberofod = 0;
    int numberofpath = 0;

    // Sparse matrix data (Eigen)
    Eigen::SparseMatrix<double> path_edge;    // |links| x |paths|
    Eigen::SparseMatrix<double> path_demand;  // |ODs| x |paths|

    Network(const std::string& name,
            const AdjList& Innode,
            const AdjList& Outnode,
            const std::vector<Link>& Link_list,
            const OdTree& Odtree);

    // Label-Correcting shortest path from origin o
    std::pair<std::map<int, double>, std::map<int, int>>
    LC(int o, const LinkMapD& weight) const;

    // Solve User Equilibrium
    void solve_ue(const DemandMap& demand, const LinkMapD& fftime,
                  const LinkMapD& capacity, int K0, double e0,
                  bool warm_start = false, const LinkMapD* toll = nullptr);

    // Solve System Optimal
    void solve_so(const DemandMap& demand, const LinkMapD& fftime,
                  const LinkMapD& capacity, int K0, double e0,
                  bool warm_start = false);

    // Enumerate paths and build index structures
    void path_enumeration();

    // Generate Eigen sparse matrices for path_edge and path_demand
    void generate_sparse_matrix();
};
