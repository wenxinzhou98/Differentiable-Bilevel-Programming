#pragma once

#include "types.h"
#include "read.h"
#include <Eigen/Dense>

class Instance {
public:
    // Raw network data (dict-based)
    LinkMapD Capacity;
    LinkMapD Length;
    LinkMapD Fftime;
    DemandMap Demand;
    std::vector<Link> Link_list;

    // Eigen vector data
    Eigen::VectorXd demand_vec;
    Eigen::VectorXd cap;
    Eigen::VectorXd tfree;

    // Toll link indices (0-based into Link_list)
    std::vector<int> toll_link;

    // Convergence parameters
    double epsilon;
    double xi;
    double eps_eva;

    // Full network data for building Network object
    AdjList Innode;
    AdjList Outnode;
    OdTree Odtree;

    Instance(const std::string& network, const std::string& base_dir);
};
