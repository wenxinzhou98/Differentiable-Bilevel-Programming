#include "instance.h"
#include <stdexcept>

Instance::Instance(const std::string& network, const std::string& base_dir) {
    // Determine which network file to read
    std::string net_name = network;
    if (network == "hearn1" || network == "hearn2") {
        net_name = "hearn";
    }

    NetworkData data = read_notoll(net_name, base_dir);

    Capacity = data.Capacity;
    Length = data.Length;
    Fftime = data.Fftime;
    Demand = data.Demand;
    Link_list = data.Link_list;
    Innode = data.Innode;
    Outnode = data.Outnode;
    Odtree = data.Odtree;

    int num_links = (int)Link_list.size();
    int num_ods = (int)Demand.size();

    // Convert to Eigen vectors (ordered by Link_list order)
    Eigen::VectorXd tfree_raw(num_links);
    Eigen::VectorXd cap_raw(num_links);
    for (int i = 0; i < num_links; ++i) {
        tfree_raw(i) = Fftime[Link_list[i]];
        cap_raw(i) = Capacity[Link_list[i]];
    }

    // Convert demand to vector (ordered by Demand map iteration order)
    demand_vec.resize(num_ods);
    int idx = 0;
    for (auto& [od, d] : Demand) {
        demand_vec(idx++) = d;
    }

    // Network-specific configuration
    if (network == "hearn1") {
        toll_link = {1, 10, 11};
        epsilon = 1e-6;
        xi = 1e-4;
        eps_eva = 1e-8;
        tfree = tfree_raw;
        cap = cap_raw;
    }
    else if (network == "hearn2") {
        toll_link = {10, 11};
        epsilon = 1e-6;
        xi = 1e-4;
        eps_eva = 1e-8;
        tfree = tfree_raw;
        cap = cap_raw;
    }
    else if (network == "sf") {
        // Scale demand and capacity
        demand_vec /= 1000.0;
        // Update Demand dict
        for (auto& [od, d] : Demand) d = demand_vec(0); // need index
        idx = 0;
        for (auto& [od, d] : Demand) {
            d = demand_vec(idx++);
        }

        tfree = tfree_raw / 0.1 * 0.06;
        // Update Fftime dict
        for (int i = 0; i < num_links; ++i) {
            Fftime[Link_list[i]] = tfree(i);
        }

        cap = cap_raw / 1000.0;
        // Update Capacity dict
        for (int i = 0; i < num_links; ++i) {
            Capacity[Link_list[i]] = cap(i);
        }

        toll_link = {11, 14, 21, 32, 35, 38, 41, 45, 46, 48, 51, 52, 57, 64, 66, 68, 70, 73};
        epsilon = 1e-6;
        xi = 1e-4;
        eps_eva = 1e-6;
    }
    else if (network == "bar") {
        demand_vec *= 2.0 / 1000.0;
        idx = 0;
        for (auto& [od, d] : Demand) {
            d = demand_vec(idx++);
        }

        tfree = tfree_raw;

        cap = cap_raw / 1000.0;
        for (int i = 0; i < num_links; ++i) {
            Capacity[Link_list[i]] = cap(i);
        }

        toll_link = {1514, 1464, 2265, 2517, 2274, 2062,  677,  704,  644,  674,  649, 2176,
                     1752, 1635, 2484, 1652, 1531, 1782, 2007, 1892};
        epsilon = 1e-4;
        xi = 1e-4;
        eps_eva = 1e-6;
    }
    else if (network == "cs") {
        demand_vec *= 2.0 / 1000.0;
        idx = 0;
        for (auto& [od, d] : Demand) {
            d = demand_vec(idx++);
        }

        tfree = tfree_raw;

        cap = cap_raw / 1000.0;
        for (int i = 0; i < num_links; ++i) {
            Capacity[Link_list[i]] = cap(i);
        }

        toll_link = {658,  566,  800, 1082,  810,  646,  654,  849, 1053,  807,  869,  857,
                     771,  789,  805,  847,  775,  650,  966,  801,  865,  796,  852,  621,
                     630,  853,  755,  861,  642,  638,  634,  782,  813,  856,  783,  561,
                     655,  662,  779,  860};
        epsilon = 1e-4;
        xi = 1e-4;
        eps_eva = 1e-6;
    }
    else {
        throw std::runtime_error("Unknown network: " + network);
    }
}
