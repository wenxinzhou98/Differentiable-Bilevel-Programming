#include "gp.h"
#include "network.h"
#include <cmath>
#include <algorithm>
#include <set>

Path sppconvert(const std::map<int, int>& predecessor, int o, int d) {
    if (d == o) return {o};
    Path convert;
    int p = d;
    convert.push_back(p);
    while (true) {
        p = predecessor.at(p);
        convert.insert(convert.begin(), p);
        if (p == o) break;
    }
    return convert;
}

double calc_pathcost(const Path& path, const LinkMapD& linkimpedance) {
    double cost = 0.0;
    for (size_t i = 0; i + 1 < path.size(); ++i) {
        Link l = {path[i], path[i + 1]};
        cost += linkimpedance.at(l);
    }
    return cost;
}

// Update link flow and time along a path
static void update_linkflowtime(const Path& spp, LinkMapD& flow, LinkMapD& time_map,
                                 const LinkMapD& fftime, const LinkMapD& capacity,
                                 double add_flow, const LinkMapD* toll) {
    for (size_t i = 0; i + 1 < spp.size(); ++i) {
        Link l = {spp[i], spp[i + 1]};
        flow[l] += add_flow;
        time_map[l] = BPR(fftime.at(l), capacity.at(l), flow[l]);
        if (toll) {
            time_map[l] += toll->at(l);
        }
    }
}

// Hessian diagonal for UE
static double calc_hkl_ue(const Path& path, const Path& spp,
                           const LinkMapD& flow, const LinkMapD& capacity,
                           const LinkMapD& fftime) {
    // Collect links in path and spp
    std::set<Link> path_links, spp_links;
    for (size_t i = 0; i + 1 < path.size(); ++i)
        path_links.insert({path[i], path[i + 1]});
    for (size_t i = 0; i + 1 < spp.size(); ++i)
        spp_links.insert({spp[i], spp[i + 1]});

    // Symmetric difference
    std::set<Link> sym_diff;
    for (auto& l : path_links)
        if (spp_links.find(l) == spp_links.end()) sym_diff.insert(l);
    for (auto& l : spp_links)
        if (path_links.find(l) == path_links.end()) sym_diff.insert(l);

    double hkl = 0.0;
    for (auto& l : sym_diff) {
        hkl += BPR_1_derivative(flow.at(l), capacity.at(l), fftime.at(l));
    }
    return hkl;
}

// Hessian diagonal for SO (uses marginal cost second derivative)
static double calc_hkl_so(const Path& path, const Path& spp,
                           const LinkMapD& flow, const LinkMapD& capacity,
                           const LinkMapD& fftime) {
    std::set<Link> path_links, spp_links;
    for (size_t i = 0; i + 1 < path.size(); ++i)
        path_links.insert({path[i], path[i + 1]});
    for (size_t i = 0; i + 1 < spp.size(); ++i)
        spp_links.insert({spp[i], spp[i + 1]});

    std::set<Link> sym_diff;
    for (auto& l : path_links)
        if (spp_links.find(l) == spp_links.end()) sym_diff.insert(l);
    for (auto& l : spp_links)
        if (path_links.find(l) == path_links.end()) sym_diff.insert(l);

    double hkl = 0.0;
    for (auto& l : sym_diff) {
        hkl += 2.0 * BPR_1_derivative(flow.at(l), capacity.at(l), fftime.at(l))
             + flow.at(l) * BPR_2_derivative(flow.at(l), capacity.at(l), fftime.at(l));
    }
    return hkl;
}

// Delete zero-flow paths
static void delete_zero_paths(PathFlowMap& pathflow, PathSetMap& pathset) {
    for (auto& [od, pf] : pathflow) {
        std::vector<Path> to_delete;
        for (auto& [path, f] : pf) {
            if (std::abs(f) < 1e-9) {
                to_delete.push_back(path);
            }
        }
        for (auto& dp : to_delete) {
            pf.erase(dp);
            pathset[od].erase(path_to_string(dp));
        }
    }
}

// Update link flow and time for SO (marginal cost)
static void update_linkflowtime_so(const Path& spp, LinkMapD& flow, LinkMapD& time_map,
                                    const LinkMapD& fftime, const LinkMapD& capacity,
                                    double add_flow) {
    for (size_t i = 0; i + 1 < spp.size(); ++i) {
        Link l = {spp[i], spp[i + 1]};
        flow[l] += add_flow;
        // Marginal cost: BPR + flow * BPR'
        time_map[l] = BPR(fftime.at(l), capacity.at(l), flow[l])
                     + flow[l] * BPR_1_derivative(flow[l], capacity.at(l), fftime.at(l));
    }
}

GPResult GP(Network& net, const DemandMap& demand, const LinkMapD& fftime,
            const LinkMapD& capacity, int K0, double e0,
            bool warm_start, const LinkMapD* toll, bool do_delete) {

    int maxInIter = 5;
    int k = 0;
    double RG = 1e30;

    PathFlowMap pathflow;
    PathSetMap pathset;
    LinkMapD flow;
    LinkMapD time_map;

    if (!warm_start) {
        // Cold start
        time_map = fftime;
        if (toll) {
            for (auto& lk : net.Link_list) {
                time_map[lk] += toll->at(lk);
            }
        }
        for (auto& lk : net.Link_list) flow[lk] = 0.0;
        for (auto& [od, d] : demand) {
            pathflow[od] = {};
            pathset[od] = {};
        }

        for (auto& [o, dests] : net.Odtree) {
            auto [impedance, predecessor] = net.LC(o, time_map);
            for (int d : dests) {
                Path spp = sppconvert(predecessor, o, d);
                std::string spp_str = path_to_string(spp);
                OD od = {o, d};
                pathset[od].insert(spp_str);
                pathflow[od][spp] = demand.at(od);
                update_linkflowtime(spp, flow, time_map, fftime, capacity, demand.at(od), toll);
            }
        }
    } else {
        // Warm start from existing solution
        pathflow = net.pathflow;
        pathset = net.pathset;
        flow = net.flow;
        for (auto& lk : net.Link_list) {
            time_map[lk] = BPR(fftime.at(lk), capacity.at(lk), flow[lk]);
        }
        if (toll) {
            for (auto& lk : net.Link_list) {
                time_map[lk] += toll->at(lk);
            }
        }
    }

    // Main loop
    while (true) {
        if (RG < 1e-6) maxInIter = 5;

        std::map<OD, Path> spps;
        double fenzi = 0.0, fenmu = 0.0;
        int numberofpath = 0;

        for (auto& [o, dests] : net.Odtree) {
            auto [impedance, predecessor] = net.LC(o, time_map);
            for (int d : dests) {
                OD od = {o, d};
                Path spp = sppconvert(predecessor, o, d);
                spps[od] = spp;
                std::string spp_str = path_to_string(spp);
                if (pathset[od].find(spp_str) == pathset[od].end()) {
                    pathset[od].insert(spp_str);
                    pathflow[od][spp] = 0.0;
                }
                fenzi += impedance.at(d) * demand.at(od);
                numberofpath += (int)pathset[od].size();
            }
        }
        for (auto& lk : net.Link_list) {
            fenmu += flow[lk] * time_map[lk];
        }
        RG = 1.0 - (fenzi / fenmu);

        if (k >= K0 || RG <= e0) {
            if (k >= K0) std::cerr << "does not converge" << std::endl;
            // Count final paths
            numberofpath = 0;
            for (auto& [od, ps] : pathset) numberofpath += (int)ps.size();
            return {pathflow, flow, numberofpath, pathset};
        }

        // Inner iterations: flow shifting
        for (int il = 0; il < maxInIter; ++il) {
            int total_shift_od_num = 0;
            for (auto& [o, dests] : net.Odtree) {
                for (int d : dests) {
                    OD od = {o, d};
                    const Path& spp = spps[od];
                    double maxCost = 0.0, minCost = 1e10;
                    for (auto& [path, f] : pathflow[od]) {
                        double cst = calc_pathcost(path, time_map);
                        maxCost = std::max(maxCost, cst);
                        minCost = std::min(minCost, cst);
                    }
                    double shiftFlow = maxCost - minCost;
                    if (shiftFlow > RG / 2.0) {
                        total_shift_od_num++;

                        // Copy keys to avoid iterator invalidation issues
                        std::vector<Path> paths;
                        for (auto& [p, f] : pathflow[od]) paths.push_back(p);

                        for (auto& path : paths) {
                            if (path != spp) {
                                double gkl = calc_pathcost(path, time_map) - calc_pathcost(spp, time_map);
                                double hkl = calc_hkl_ue(path, spp, flow, capacity, fftime);
                                double f_path = pathflow[od][path];
                                double f_spp = pathflow[od][spp];
                                double delta = std::max(f_spp + (f_path - std::max(0.0, f_path - gkl / hkl)), 0.0) - f_spp;
                                pathflow[od][path] -= delta;
                                pathflow[od][spp] += delta;
                                shiftFlow += std::abs(delta);
                                update_linkflowtime(path, flow, time_map, fftime, capacity, -delta, toll);
                                update_linkflowtime(spp, flow, time_map, fftime, capacity, delta, toll);
                            }
                        }
                    }
                }
            }
            if (total_shift_od_num < 3) break;
        }

        if (do_delete) {
            delete_zero_paths(pathflow, pathset);
        }
        k++;
    }
}

GPResult GP_SO(Network& net, const DemandMap& demand, const LinkMapD& fftime,
               const LinkMapD& capacity, int K0, double e0, bool warm_start) {

    int maxInIter = 5;
    int k = 0;
    double RG = 1e30;

    PathFlowMap pathflow;
    PathSetMap pathset;
    LinkMapD flow;
    LinkMapD time_map;

    if (!warm_start) {
        // Cold start with marginal cost
        for (auto& lk : net.Link_list) {
            time_map[lk] = fftime.at(lk);  // initial marginal cost = fftime
            flow[lk] = 0.0;
        }
        for (auto& [od, d] : demand) {
            pathflow[od] = {};
            pathset[od] = {};
        }

        for (auto& [o, dests] : net.Odtree) {
            auto [impedance, predecessor] = net.LC(o, time_map);
            for (int d : dests) {
                Path spp = sppconvert(predecessor, o, d);
                std::string spp_str = path_to_string(spp);
                OD od = {o, d};
                pathset[od].insert(spp_str);
                pathflow[od][spp] = demand.at(od);
                update_linkflowtime_so(spp, flow, time_map, fftime, capacity, demand.at(od));
            }
        }
    } else {
        pathflow = net.pathflow;
        pathset = net.pathset;
        flow = net.flow;
        for (auto& lk : net.Link_list) {
            time_map[lk] = BPR(fftime.at(lk), capacity.at(lk), flow[lk])
                         + flow[lk] * BPR_1_derivative(flow[lk], capacity.at(lk), fftime.at(lk));
        }
    }

    // Main loop
    while (true) {
        if (RG < 1e-6) maxInIter = 5;

        std::map<OD, Path> spps;
        double fenzi = 0.0, fenmu = 0.0;
        int numberofpath = 0;

        for (auto& [o, dests] : net.Odtree) {
            auto [impedance, predecessor] = net.LC(o, time_map);
            for (int d : dests) {
                OD od = {o, d};
                Path spp = sppconvert(predecessor, o, d);
                spps[od] = spp;
                std::string spp_str = path_to_string(spp);
                if (pathset[od].find(spp_str) == pathset[od].end()) {
                    pathset[od].insert(spp_str);
                    pathflow[od][spp] = 0.0;
                }
                fenzi += impedance.at(d) * demand.at(od);
                numberofpath += (int)pathset[od].size();
            }
        }
        for (auto& lk : net.Link_list) {
            fenmu += flow[lk] * time_map[lk];
        }
        RG = 1.0 - (fenzi / fenmu);

        if (k >= K0 || RG <= e0) {
            numberofpath = 0;
            for (auto& [od, ps] : pathset) numberofpath += (int)ps.size();
            return {pathflow, flow, numberofpath, pathset};
        }

        for (int il = 0; il < maxInIter; ++il) {
            int total_shift_od_num = 0;
            for (auto& [o, dests] : net.Odtree) {
                for (int d : dests) {
                    OD od = {o, d};
                    const Path& spp = spps[od];
                    double maxCost = 0.0, minCost = 1e10;
                    for (auto& [path, f] : pathflow[od]) {
                        double cst = calc_pathcost(path, time_map);
                        maxCost = std::max(maxCost, cst);
                        minCost = std::min(minCost, cst);
                    }
                    double shiftFlow = maxCost - minCost;
                    if (shiftFlow > RG / 2.0) {
                        total_shift_od_num++;

                        std::vector<Path> paths;
                        for (auto& [p, f] : pathflow[od]) paths.push_back(p);

                        for (auto& path : paths) {
                            if (path != spp) {
                                double gkl = calc_pathcost(path, time_map) - calc_pathcost(spp, time_map);
                                double hkl = calc_hkl_so(path, spp, flow, capacity, fftime);
                                double f_path = pathflow[od][path];
                                double f_spp = pathflow[od][spp];
                                double delta = std::max(f_spp + (f_path - std::max(0.0, f_path - gkl / hkl)), 0.0) - f_spp;
                                pathflow[od][path] -= delta;
                                pathflow[od][spp] += delta;
                                shiftFlow += std::abs(delta);
                                update_linkflowtime_so(path, flow, time_map, fftime, capacity, -delta);
                                update_linkflowtime_so(spp, flow, time_map, fftime, capacity, delta);
                            }
                        }
                    }
                }
            }
            if (total_shift_od_num < 3) break;
        }

        // Delete zero-flow paths
        delete_zero_paths(pathflow, pathset);
        k++;
    }
}
