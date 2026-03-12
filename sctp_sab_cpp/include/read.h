#pragma once

#include "types.h"
#include <string>

struct NetworkData {
    int node_count;
    int link_count;
    double flow_count;

    LinkMapD Capacity;
    LinkMapD Length;
    LinkMapD Fftime;
    AdjList Innode;
    AdjList Outnode;
    std::vector<Link> Link_list;
    OdTree Odtree;
    DemandMap Demand;
};

// Read network topology from .net file
void read_net(const std::string& network, const std::string& base_dir,
              LinkMapD& Capacity, LinkMapD& Length, LinkMapD& Fftime,
              AdjList& Innode, AdjList& Outnode, std::vector<Link>& Link_list,
              int& node_count, int& link_count);

// Read trip (OD demand) from .trp file
void read_trp(const std::string& network, const std::string& base_dir,
              OdTree& Odtree, DemandMap& Demand, double& flow_count);

// Combined read without toll
NetworkData read_notoll(const std::string& network, const std::string& base_dir);
