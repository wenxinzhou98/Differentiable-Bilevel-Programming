#include "network.h"
#include "gp.h"
#include <deque>
#include <algorithm>

Network::Network(const std::string& name,
                 const AdjList& Innode,
                 const AdjList& Outnode,
                 const std::vector<Link>& Link_list,
                 const OdTree& Odtree)
    : name(name), Innode(Innode), Outnode(Outnode),
      Link_list(Link_list), Odtree(Odtree) {}

std::pair<std::map<int, double>, std::map<int, int>>
Network::LC(int o, const LinkMapD& weight) const {
    std::map<int, double> impedance;
    std::map<int, int> predecessor;

    for (auto& [node, _] : Outnode) {
        impedance[node] = (node == o) ? 0.0 : 1e30;
        predecessor[node] = -1;
    }
    // Also initialize nodes that only appear in Innode
    for (auto& [node, _] : Innode) {
        if (impedance.find(node) == impedance.end()) {
            impedance[node] = 1e30;
            predecessor[node] = -1;
        }
    }

    std::deque<int> Q;
    Q.push_back(o);

    while (!Q.empty()) {
        int i = Q.front();
        Q.pop_front();

        auto it = Outnode.find(i);
        if (it == Outnode.end()) continue;

        for (int j : it->second) {
            Link l = {i, j};
            auto wit = weight.find(l);
            if (wit == weight.end()) continue;

            if (impedance[j] > impedance[i] + wit->second) {
                impedance[j] = impedance[i] + wit->second;
                predecessor[j] = i;
                if (std::find(Q.begin(), Q.end(), j) == Q.end()) {
                    Q.push_back(j);
                }
            }
        }
    }

    return {impedance, predecessor};
}

void Network::solve_ue(const DemandMap& demand, const LinkMapD& fftime,
                       const LinkMapD& capacity, int K0, double e0,
                       bool warm_start, const LinkMapD* toll) {
    auto result = GP(*this, demand, fftime, capacity, K0, e0, warm_start, toll);
    pathflow = std::move(result.pathflow);
    flow = std::move(result.flow);
    path_number = result.numberofpath;
    pathset = std::move(result.pathset);
    demand_map = demand;

    int ii = 0;
    for (auto& [od, d] : demand) {
        od2num[od] = ii++;
    }
}

void Network::solve_so(const DemandMap& demand, const LinkMapD& fftime,
                       const LinkMapD& capacity, int K0, double e0,
                       bool warm_start) {
    auto result = GP_SO(*this, demand, fftime, capacity, K0, e0, warm_start);
    pathflow = std::move(result.pathflow);
    flow = std::move(result.flow);
    path_number = result.numberofpath;
    pathset = std::move(result.pathset);
}

void Network::path_enumeration() {
    numberoflink = (int)Link_list.size();
    numberofod = (int)pathflow.size();
    numberofpath = 0;

    // Build link index
    std::map<Link, int> link_index;
    for (int i = 0; i < (int)Link_list.size(); ++i) {
        link_index[Link_list[i]] = i;
    }

    // Temporary storage for triplets
    std::vector<Eigen::Triplet<double>> edge_triplets;
    std::vector<Eigen::Triplet<double>> demand_triplets;

    path_id.clear();

    int kindex = 0;
    int odindex = 0;

    for (auto& [o, dests] : Odtree) {
        for (int d : dests) {
            OD od = {o, d};
            path_id[od].clear();

            if (pathflow.find(od) == pathflow.end()) {
                odindex++;
                continue;
            }

            for (auto& [path, f] : pathflow[od]) {
                numberofpath++;
                for (size_t i = 0; i + 1 < path.size(); ++i) {
                    Link lk = {path[i], path[i + 1]};
                    int lk_idx = link_index[lk];
                    edge_triplets.emplace_back(lk_idx, kindex, 1.0);
                }
                demand_triplets.emplace_back(odindex, kindex, 1.0);
                path_id[od][path] = kindex;
                kindex++;
            }
            odindex++;
        }
    }

    path_number = numberofpath;

    // Store triplets for generate_sparse_matrix
    // We generate the matrices right here to keep them in sync
    path_edge.resize(numberoflink, numberofpath);
    path_edge.setFromTriplets(edge_triplets.begin(), edge_triplets.end());

    path_demand.resize(numberofod, numberofpath);
    path_demand.setFromTriplets(demand_triplets.begin(), demand_triplets.end());
}

void Network::generate_sparse_matrix() {
    // Already generated in path_enumeration() for efficiency
    // This method exists for API compatibility with the Python version
    path_edge.makeCompressed();
    path_demand.makeCompressed();
}
