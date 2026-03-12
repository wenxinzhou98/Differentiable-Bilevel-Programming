// sctp_sab.cpp
//
// C++ translation of sctp_sab.py
//
// Algorithm: Sensitivity Analysis Based (SAB) bilevel optimization for the
//            Single-level Congestion Toll Problem (SCTP).
//
// Matrix computations: Eigen3
// LP solver:           Gurobi
// UE solver:           Gradient Projection (GP) from GP.py / graph.py
//
// Original Python author: Jiayang Li

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <gurobi_c++.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace Eigen;
using namespace std;

// ============================================================
// Type aliases
// ============================================================

using Link = pair<int, int>;  // directed network edge (tail, head)
using OD   = pair<int, int>;  // origin-destination pair
using Path = vector<int>;     // sequence of node ids

// ============================================================
// BPR travel time functions  (from GP.py)
// ============================================================

// Bureau of Public Roads travel time: t0 * (1 + alpha*(x/c)^beta)
double BPR(double fftime, double capacity, double flow,
           double alpha = 0.15, double beta = 4.0) {
    return fftime * (1.0 + alpha * pow(flow / capacity, beta));
}

// First derivative of BPR w.r.t. flow:  4 * t0 * 0.15 * x^3 / c^4
double BPR_1_derivative(double x_a, double c_a, double t0) {
    return (4.0 * t0 * 0.15 / pow(c_a, 4.0)) * pow(x_a, 3.0);
}

// ============================================================
// Path utilities
// ============================================================

// Convert path to hyphen-separated string (e.g. "1-3-5") — used for pathset keys
string path_to_string(const Path& p) {
    string s;
    for (size_t i = 0; i < p.size(); ++i) {
        if (i > 0) s += '-';
        s += to_string(p[i]);
    }
    return s;
}

// Reconstruct a path from the predecessor map produced by LC/LS
// Mirrors sppconvert() in GP.py
Path sppconvert(const map<int, int>& predecessor, int o, int d) {
    if (d == o) return {o};
    Path path;
    int p = d;
    path.push_back(p);
    while (p != o) {
        p = predecessor.at(p);
        path.insert(path.begin(), p);
    }
    return path;
}

// ============================================================
// Network struct  (mirrors graph.py Network class)
// ============================================================

struct Network {
    string name;

    // Topology
    map<int, vector<int>> Innode;   // node -> predecessor nodes
    map<int, vector<int>> Outnode;  // node -> successor nodes
    vector<Link>           Link_list; // ordered list of all directed links

    // OD tree:  origin -> list of destinations
    map<int, vector<int>> Odtree;

    // UE solution (populated by GP / solve_ue)
    map<OD, map<Path, double>> pathflow;  // (o,d) -> { path -> flow }
    map<Link, double>           flow;     // link -> aggregate flow
    int                         path_number = 0;
    map<OD, set<string>>        pathset;  // (o,d) -> set of path-strings

    // OD index:  od pair -> row index used in path_demand
    map<OD, int> od2num;

    // Sparse path-edge / path-demand index lists (populated by path_enumeration)
    vector<int> path_edge_lk;        // row (link) indices for path_edge
    vector<int> path_edge_kindex;    // col (path) indices for path_edge
    vector<int> path_demand_odindex; // row (OD)   indices for path_demand
    vector<int> path_demand_kindex;  // col (path) indices for path_demand
    int numberoflink = 0, numberofod = 0, numberofpath = 0;

    map<OD, map<Path, int>> path_id; // (o,d) -> { path -> column index }
    vector<int>              path2od; // column k -> OD row index

    // Eigen sparse matrices (populated by generate_sparse_matrix)
    SparseMatrix<double> path_edge;   // [n_links x n_paths]
    SparseMatrix<double> path_demand; // [n_OD   x n_paths]

    // Demand dictionary (stored for warm-start access inside GP)
    map<OD, double> demand_dict;
};

// ============================================================
// Label Correcting (LC) shortest-path  (from graph.py)
// Returns {impedance map, predecessor map}
// ============================================================
pair<map<int, double>, map<int, int>>
LC(int o,
   const map<Link, double>& weight,
   const map<int, vector<int>>& Outnode)
{
    map<int, double> impedance;
    map<int, int>    predecessor;

    for (const auto& kv : Outnode) {
        int node = kv.first;
        impedance[node]  = (node == o) ? 0.0 : numeric_limits<double>::infinity();
        predecessor[node] = -1;
    }

    vector<int> Q = {o};
    while (!Q.empty()) {
        int i = Q.front();
        Q.erase(Q.begin());
        for (int j : Outnode.at(i)) {
            Link l = {i, j};
            double new_dist = impedance[i] + weight.at(l);
            if (impedance[j] > new_dist) {
                impedance[j]  = new_dist;
                predecessor[j] = i;
                if (find(Q.begin(), Q.end(), j) == Q.end())
                    Q.push_back(j);
            }
        }
    }
    return {impedance, predecessor};
}

// ============================================================
// update_linkflowtime  (from GP.py)
// Add add_flow along path spp, recompute BPR times.
// ============================================================
void update_linkflowtime(const Path& spp,
                         map<Link, double>& flow,
                         map<Link, double>& time_map,
                         const map<Link, double>& fftime,
                         const map<Link, double>& capacity,
                         double add_flow,
                         const map<Link, double>* toll = nullptr)
{
    for (size_t k = 0; k + 1 < spp.size(); ++k) {
        Link l = {spp[k], spp[k + 1]};
        flow[l]    += add_flow;
        time_map[l] = BPR(fftime.at(l), capacity.at(l), flow[l], 0.15, 4);
        if (toll) time_map[l] += toll->at(l);
    }
}

// ============================================================
// calc_pathcost  (from GP.py)
// ============================================================
double calc_pathcost(const Path& path, const map<Link, double>& linkimpedance) {
    double cost = 0.0;
    for (size_t k = 0; k + 1 < path.size(); ++k)
        cost += linkimpedance.at({path[k], path[k + 1]});
    return cost;
}

// ============================================================
// calc_hkl  (from GP.py)
// Sum BPR first derivatives over the symmetric difference of the
// link sets of `path` and `spp`.
// Note: the `isinspp` parameter exists in the Python signature but is unused
//       in the function body — it is omitted here.
// ============================================================
double calc_hkl(const Path& path,
                const Path& spp,
                const map<Link, double>& flow,
                const map<Link, double>& capacity,
                const map<Link, double>& fftime)
{
    set<Link> path_link, spp_link;
    for (size_t k = 0; k + 1 < path.size(); ++k)
        path_link.insert({path[k], path[k + 1]});
    for (size_t k = 0; k + 1 < spp.size(); ++k)
        spp_link.insert({spp[k], spp[k + 1]});

    // Symmetric difference: links in exactly one of the two paths
    set<Link> sym_diff;
    set_symmetric_difference(path_link.begin(), path_link.end(),
                             spp_link.begin(), spp_link.end(),
                             inserter(sym_diff, sym_diff.begin()));

    double hkl = 0.0;
    for (const Link& l : sym_diff)
        hkl += BPR_1_derivative(flow.at(l), capacity.at(l), fftime.at(l));
    return hkl;
}

// ============================================================
// calc_obj  (from GP.py)
// Compute BPR objective (used internally by GP for tracking convergence).
// ============================================================
double calc_obj(const map<Link, double>& flow,
                const map<Link, double>& fftime,
                const map<Link, double>& capacity,
                const map<Link, double>* toll = nullptr)
{
    double obj = 0.0;
    for (const auto& kv : flow) {
        const Link& l = kv.first;
        double fl = kv.second;
        obj += fftime.at(l) * fl
             + 0.15 * fftime.at(l) * pow(fl, 5.0) / (5.0 * pow(capacity.at(l), 4.0));
        if (toll) obj += fl * toll->at(l);
    }
    return obj;
}

// ============================================================
// GP: Gradient Projection User Equilibrium solver  (from GP.py)
// Populates net.pathflow, net.flow, net.path_number, net.pathset, net.od2num.
// ============================================================
void GP(Network& net,
        const map<OD, double>&   demand,
        const map<Link, double>& fftime,
        const map<Link, double>& capacity,
        int    K0,
        double e0,
        bool   warm_start    = false,
        const map<Link, double>* toll = nullptr,
        bool   delete_paths  = true)
{
    int    maxInIter = 5;
    int    k         = 0;
    double RG        = numeric_limits<double>::infinity();

    map<OD, map<Path, double>> pathflow;
    map<OD, set<string>>       pathset;
    map<Link, double>          flow, time_map;

    // ---- Initialization ----
    if (!warm_start) {
        // Cold start: initialise from all-shortest-path assignment
        for (const Link& l : net.Link_list) flow[l] = 0.0;
        for (const Link& l : net.Link_list)
            time_map[l] = fftime.at(l) + (toll ? toll->at(l) : 0.0);

        for (const auto& kv : demand) { pathflow[kv.first] = {}; pathset[kv.first] = {}; }

        for (const auto& kv1 : net.Odtree) {
            int o = kv1.first;
            auto [imp, pred] = LC(o, time_map, net.Outnode);
            for (int d : kv1.second) {
                OD   od      = {o, d};
                Path spp     = sppconvert(pred, o, d);
                string spp_s = path_to_string(spp);
                pathset[od].insert(spp_s);
                pathflow[od][spp] = demand.at(od);
                update_linkflowtime(spp, flow, time_map, fftime, capacity, demand.at(od), toll);
            }
        }
    } else {
        // Warm start: resume from existing net state
        pathflow = net.pathflow;
        pathset  = net.pathset;
        flow     = net.flow;
        for (const Link& l : net.Link_list)
            time_map[l] = BPR(fftime.at(l), capacity.at(l), flow[l], 0.15, 4)
                        + (toll ? toll->at(l) : 0.0);
    }

    // Ensure pathflow has entries for all OD pairs
    for (const auto& kv1 : net.Odtree)
        for (int d : kv1.second)
            if (!pathflow.count({kv1.first, d})) {
                pathflow[{kv1.first, d}] = {};
                pathset [{kv1.first, d}] = {};
            }

    // ---- Main GP loop ----
    while (true) {
        if (RG < 1e-6) maxInIter = 5;

        map<OD, Path> spps;
        double fenzi = 0.0, fenmu = 0.0;

        for (const auto& kv1 : net.Odtree) {
            int o = kv1.first;
            auto [imp, pred] = LC(o, time_map, net.Outnode);
            for (int d : kv1.second) {
                OD   od      = {o, d};
                Path spp     = sppconvert(pred, o, d);
                spps[od]     = spp;
                string spp_s = path_to_string(spp);
                if (!pathset[od].count(spp_s)) {
                    pathset[od].insert(spp_s);
                    pathflow[od][spp] = 0.0;
                }
                fenzi += imp.at(d) * demand.at(od);
            }
        }
        for (const Link& l : net.Link_list)
            fenmu += flow.at(l) * time_map.at(l);
        RG = 1.0 - fenzi / fenmu;

        // (calc_obj tracked internally; not exposed to caller in sctp_sab.py)
        double obj = calc_obj(flow, fftime, capacity, toll);
        (void)obj;

        if (k >= K0 || RG <= e0) {
            if (k >= K0) cout << "does not converge" << endl;
            break;
        }

        // ---- Inner gradient projection loop ----
        for (int il = 0; il < maxInIter; ++il) {
            int total_shift_od_num = 0;
            for (const auto& kv1 : net.Odtree) {
                int o = kv1.first;
                for (int d : kv1.second) {
                    OD          od  = {o, d};
                    const Path& spp = spps[od];

                    double maxCost = 0.0, minCost = 1e18;
                    for (const auto& pkv : pathflow[od]) {
                        double cst = calc_pathcost(pkv.first, time_map);
                        maxCost = max(maxCost, cst);
                        minCost = min(minCost, cst);
                    }
                    double shiftFlow = maxCost - minCost;

                    if (shiftFlow > RG / 2.0) {
                        total_shift_od_num++;

                        // Collect paths to avoid iterator invalidation
                        vector<Path> other_paths;
                        for (const auto& pkv : pathflow[od])
                            if (pkv.first != spp) other_paths.push_back(pkv.first);

                        for (const Path& path : other_paths) {
                            double gkl = calc_pathcost(path, time_map)
                                       - calc_pathcost(spp,  time_map);
                            double hkl = calc_hkl(path, spp, flow, capacity, fftime);
                            if (hkl == 0.0) continue;

                            // delta = flow shifted from path -> spp
                            // Python: delta = max(f_spp + (f_path - max(0, f_path - gkl/hkl)), 0) - f_spp
                            double f_path = pathflow[od][path];
                            double f_spp  = pathflow[od][spp];
                            double delta  = max(f_spp + (f_path - max(0.0, f_path - gkl / hkl)), 0.0) - f_spp;

                            pathflow[od][path] -= delta;
                            pathflow[od][spp]  += delta;
                            shiftFlow += abs(delta);
                            update_linkflowtime(path, flow, time_map, fftime, capacity, -delta, toll);
                            update_linkflowtime(spp,  flow, time_map, fftime, capacity,  delta, toll);
                        }
                    }
                }
            }
            if (total_shift_od_num < 3) break;
        }

        // Delete near-zero paths
        auto delete_zero_paths = [&]() {
            for (const auto& kv1 : net.Odtree) {
                int o = kv1.first;
                for (int d : kv1.second) {
                    OD od = {o, d};
                    vector<Path> to_delete;
                    for (const auto& pkv : pathflow[od])
                        if (abs(pkv.second) < 1e-9) to_delete.push_back(pkv.first);
                    for (const Path& p : to_delete) {
                        pathflow[od].erase(p);
                        pathset[od].erase(path_to_string(p));
                    }
                }
            }
        };
        if (delete_paths) delete_zero_paths();

        k++;
    }

    // Final deletion
    if (delete_paths) {
        for (const auto& kv1 : net.Odtree) {
            int o = kv1.first;
            for (int d : kv1.second) {
                OD od = {o, d};
                vector<Path> to_delete;
                for (const auto& pkv : pathflow[od])
                    if (abs(pkv.second) < 1e-9) to_delete.push_back(pkv.first);
                for (const Path& p : to_delete) {
                    pathflow[od].erase(p);
                    pathset[od].erase(path_to_string(p));
                }
            }
        }
    }

    // Count total path number
    int total_paths = 0;
    for (const auto& kv1 : net.Odtree)
        for (int d : kv1.second)
            total_paths += (int)pathset[{kv1.first, d}].size();

    // ---- Write results back to net ----
    net.pathflow     = pathflow;
    net.flow         = flow;
    net.path_number  = total_paths;
    net.pathset      = pathset;
    net.demand_dict  = demand;

    // od2num: assign OD row indices in the same iteration order as path_enumeration
    net.od2num.clear();
    int ii = 0;
    for (const auto& kv1 : net.Odtree)
        for (int d : kv1.second)
            net.od2num[{kv1.first, d}] = ii++;
}

// ============================================================
// path_enumeration  (from graph.py)
// Builds sparse index lists describing path_edge and path_demand.
// ============================================================
void path_enumeration(Network& net) {
    net.numberoflink = (int)net.Link_list.size();
    net.numberofod   = (int)net.od2num.size();
    net.numberofpath = 0;

    net.path_edge_lk.clear();
    net.path_edge_kindex.clear();
    net.path_demand_odindex.clear();
    net.path_demand_kindex.clear();
    net.path2od.clear();
    net.path_id.clear();

    // Map each link to its column index in path_edge
    map<Link, int> link_index;
    for (int i = 0; i < (int)net.Link_list.size(); ++i)
        link_index[net.Link_list[i]] = i;

    int kindex = 0, odindex = 0;
    for (const auto& kv1 : net.Odtree) {
        int o = kv1.first;
        for (int d : kv1.second) {
            OD od = {o, d};
            net.path_id[od] = {};

            if (net.pathflow.count(od)) {
                for (const auto& pkv : net.pathflow[od]) {
                    const Path& path = pkv.first;
                    net.numberofpath++;

                    // path-edge incidence: path passes through each consecutive link
                    for (size_t i = 0; i + 1 < path.size(); ++i) {
                        int lk = link_index.at({path[i], path[i + 1]});
                        net.path_edge_lk.push_back(lk);
                        net.path_edge_kindex.push_back(kindex);
                    }
                    // path-demand incidence: path belongs to this OD pair
                    net.path_demand_odindex.push_back(odindex);
                    net.path_demand_kindex.push_back(kindex);

                    net.path_id[od][path] = kindex;
                    net.path2od.push_back(net.od2num.at(od));
                    kindex++;
                }
            }
            odindex++;
        }
    }
}

// ============================================================
// generate_sparse_matrix  (from graph.py)
// Fills net.path_edge [n_links x n_paths] and
//       net.path_demand [n_OD   x n_paths] as Eigen sparse matrices.
// (Mirrors PyTorch sparse tensor construction from index lists in graph.py.)
// Uses Triplet-based construction for efficient batch insertion.
// ============================================================
void generate_sparse_matrix(Network& net) {
    using T = Triplet<double>;

    // path_edge: entry (link_row, path_col) = 1  for each link on each path
    vector<T> edge_triplets;
    edge_triplets.reserve(net.path_edge_lk.size());
    for (size_t i = 0; i < net.path_edge_lk.size(); ++i)
        edge_triplets.emplace_back(net.path_edge_lk[i], net.path_edge_kindex[i], 1.0);
    net.path_edge.resize(net.numberoflink, net.numberofpath);
    net.path_edge.setFromTriplets(edge_triplets.begin(), edge_triplets.end());

    // path_demand: entry (OD_row, path_col) = 1  for the OD pair each path serves
    vector<T> demand_triplets;
    demand_triplets.reserve(net.path_demand_odindex.size());
    for (size_t i = 0; i < net.path_demand_odindex.size(); ++i)
        demand_triplets.emplace_back(net.path_demand_odindex[i], net.path_demand_kindex[i], 1.0);
    net.path_demand.resize(net.numberofod, net.numberofpath);
    net.path_demand.setFromTriplets(demand_triplets.begin(), demand_triplets.end());
}

// ============================================================
// total_travel_time  (from sctp_sab.py)
// BPR total travel time: sum_a x_a * t0_a * (1 + 0.15*(x_a/c_a)^4)
// ============================================================
double total_travel_time(const VectorXd& x,
                         const VectorXd& tfree,
                         const VectorXd& cap)
{
    VectorXd ratio = x.cwiseQuotient(cap);
    VectorXd t     = tfree.cwiseProduct(
        VectorXd::Ones(x.size()) + 0.15 * ratio.array().pow(4).matrix());
    return x.dot(t);
}

// ============================================================
// find_linearly_independent_columns  (from sctp_sab.py)
// Uses QR decomposition with column pivoting to find a maximal set of
// linearly independent columns of A.
// Mirrors:  Q, R, P = scipy.linalg.qr(A, pivoting=True)
//           rank = sum(abs(diag(R)) > 1e-20)
//           return P[:rank]
// ============================================================
vector<int> find_linearly_independent_columns(const MatrixXd& A) {
    ColPivHouseholderQR<MatrixXd> qr(A);
    // Upper-triangular factor R
    MatrixXd R = qr.matrixR().triangularView<Upper>();
    const auto& perm = qr.colsPermutation();

    int rank = 0;
    int diag_size = (int)min(R.rows(), R.cols());
    for (int i = 0; i < diag_size; ++i)
        if (abs(R(i, i)) > 1e-20) rank++;

    vector<int> indices(rank);
    for (int i = 0; i < rank; ++i)
        indices[i] = perm.indices()(i);
    return indices;
}

// ============================================================
// solve_lp_active_paths  (Gurobi LP)
// Mirrors:  res = scipy.optimize.linprog(obj=0, A_eq=A_eq, b_eq=b_eq,
//                                        bounds=(0, None))
//           p = res.x;  inds = (p > 0)
//
// Solves:  min  0
//          s.t. [path_edge  ] * p == x           (link flow conservation)
//               [path_demand] * p == demand_vec  (OD demand satisfaction)
//               p >= 0
//
// Constraints are built by iterating sparse non-zeros via InnerIterator,
// avoiding materialisation of the full dense A_eq matrix.
// env is passed by reference and reused across calls (created once in
// run_sctp_sab) to avoid per-call licence-check overhead.
// Returns the optimal p vector (zeros and a warning if the LP fails).
// ============================================================
VectorXd solve_lp_active_paths(GRBEnv&                     env,
                                const SparseMatrix<double>& path_edge,
                                const SparseMatrix<double>& path_demand,
                                const VectorXd&             x,
                                const VectorXd&             demand_vec,
                                int                         n_paths)
{
    int n_links = (int)x.size();
    int n_OD    = (int)demand_vec.size();

    GRBModel model(env);

    // Loosen feasibility tolerance to absorb floating-point drift that
    // accumulates in GP link-flow updates (flow[l] += add_flow repeated many
    // times).  Gurobi's default 1e-6 is too tight for flows computed this way.
    model.set(GRB_DoubleParam_FeasibilityTol, 1e-4);

    // Decision variables: p[j] >= 0, objective coefficient = 0
    vector<GRBVar> p(n_paths);
    for (int j = 0; j < n_paths; ++j)
        p[j] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
    model.update();

    // Pre-allocate one linear expression per constraint row
    vector<GRBLinExpr> link_expr(n_links), od_expr(n_OD);

    // Accumulate coefficients by iterating sparse non-zeros column-by-column.
    // SparseMatrix<double> is column-major by default, so InnerIterator over
    // column j is maximally cache-friendly.
    for (int j = 0; j < n_paths; ++j) {
        for (SparseMatrix<double>::InnerIterator it(path_edge, j); it; ++it)
            link_expr[it.row()] += it.value() * p[j];
        for (SparseMatrix<double>::InnerIterator it(path_demand, j); it; ++it)
            od_expr[it.row()] += it.value() * p[j];
    }

    // Add equality constraints: path_edge * p == x
    for (int i = 0; i < n_links; ++i)
        model.addConstr(link_expr[i] == x(i));
    // Add equality constraints: path_demand * p == demand_vec
    for (int i = 0; i < n_OD; ++i)
        model.addConstr(od_expr[i] == demand_vec(i));

    // Minimise zero (feasibility LP)
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    model.optimize();

    int status = model.get(GRB_IntAttr_Status);
    VectorXd result = VectorXd::Zero(n_paths);
    if (status == GRB_OPTIMAL) {
        for (int j = 0; j < n_paths; ++j)
            result(j) = p[j].get(GRB_DoubleAttr_X);
    } else {
        cerr << "[solve_lp_active_paths] Gurobi status = " << status
             << " (not optimal); returning zero vector." << endl;
    }
    return result;
}

// ============================================================
// Helper: extract a subset of columns from a sparse matrix into a dense matrix.
// The result is MatrixXd because the extracted active-path submatrices (path_edge_i,
// path_demand_i) are small and used exclusively in dense operations downstream.
// InnerIterator visits only non-zeros, so zero-initialisation of result is correct.
// ============================================================
MatrixXd extract_columns(const SparseMatrix<double>& M, const vector<int>& cols) {
    MatrixXd result = MatrixXd::Zero(M.rows(), (int)cols.size());
    for (int i = 0; i < (int)cols.size(); ++i)
        for (SparseMatrix<double>::InnerIterator it(M, cols[i]); it; ++it)
            result(it.row(), i) = it.value();
    return result;
}

// ============================================================
// NetworkParams  (mirrors sctp_load.py Instance)
// Problem parameters for one network instance.
// ============================================================
struct NetworkParams {
    // Dicts used by GP (link/OD keyed)
    map<OD,   double> Demand;
    map<Link, double> Fftime;
    map<Link, double> Capacity;

    // Vectors aligned with net.Link_list / demand order (for Eigen ops)
    VectorXd demand;  // OD demand values
    VectorXd cap;     // link capacities
    VectorXd tfree;   // link free-flow travel times

    // Toll structure
    vector<int> toll_link;  // indices (into Link_list) of toll-able links
    vector<int> non_link;   // indices of links that must have zero toll

    // Solver parameters
    double alpha   = 0.05;   // base step size
    double beta    = 20.0;   // step decay parameter
    double epsilon = 1e-4;   // UE convergence tolerance
    double xi      = 1e-4;   // toll convergence threshold
    double eps_eva = 1e-8;   // final evaluation UE tolerance

    int n_links = 0;
    int n_OD    = 0;
};

// ============================================================
// run_sctp_sab  (main loop of sctp_sab.py)
//
// For each initialisation:
//   Iterates: solve UE → path enumeration → sensitivity analysis
//             → gradient computation → toll update → convergence check
// ============================================================
void run_sctp_sab(Network&                net,
                  const NetworkParams&    params,
                  const vector<VectorXd>& initializations,  // toll_add per init
                  const VectorXd&         x_ue,
                  const VectorXd&         x_so)
{
    const VectorXd& tfree      = params.tfree;
    const VectorXd& cap        = params.cap;
    const VectorXd& demand_vec = params.demand;
    int n_links = params.n_links;
    int n_OD    = params.n_OD;

    // Create Gurobi environment once for all LP solves in this run.
    // Avoids per-call licence-check overhead and prevents occasional
    // initialisation failures when GRBEnv is constructed inside a tight loop.
    GRBEnv grb_env = GRBEnv(true);
    grb_env.set(GRB_IntParam_OutputFlag, 0);
    grb_env.start();

    // Pre-compute UE and SO total travel times  (from sctp_sab.py top of loop)
    double TT_ue = total_travel_time(x_ue, tfree, cap);
    double TT_so = total_travel_time(x_so, tfree, cap);

    for (int n_ini = 0; n_ini < (int)initializations.size(); ++n_ini) {

        // ---- Initialise toll vector ----
        // Python: toll = zeros(len(tfree)); toll[toll_link] = toll_add
        VectorXd toll = VectorXd::Zero(n_links);
        const VectorXd& toll_add = initializations[n_ini];
        for (int i = 0; i < (int)params.toll_link.size(); ++i)
            toll(params.toll_link[i]) = toll_add(i);

        VectorXd toll_old = toll;

        auto tic = chrono::high_resolution_clock::now();
        double inverting_time = 0.0;
        int    iter_num       = 0;

        // ---- Main iteration (while True in Python) ----
        while (true) {

            // Build toll dict aligned with Link_list
            map<Link, double> Toll;
            for (int kk = 0; kk < n_links; ++kk)
                Toll[net.Link_list[kk]] = toll(kk);

            // 1. Solve UE with current tolls
            GP(net, params.Demand, params.Fftime, params.Capacity,
               1000, params.epsilon, /*warm_start=*/true, &Toll);

            // 2. Path enumeration and dense matrix generation
            path_enumeration(net);
            generate_sparse_matrix(net);

            int path_number = net.path_number;

            // path_edge   [n_links x n_paths]  — sparse
            // path_demand [n_OD   x n_paths]  — sparse
            const SparseMatrix<double>& path_edge   = net.path_edge;
            const SparseMatrix<double>& path_demand = net.path_demand;

            auto tic2 = chrono::high_resolution_clock::now();

            // 3. Extract link flows:  x = [net.flow[Link_list[i]]]
            VectorXd x(n_links);
            for (int i = 0; i < n_links; ++i)
                x(i) = net.flow.at(net.Link_list[i]);

            // 4. Solve LP to find active paths (p > 0)
            //    min 0  s.t.  path_edge*p == x,  path_demand*p == demand_vec,  p >= 0
            //    Constraints built from sparse non-zeros; no dense A_eq materialised.
            VectorXd p = solve_lp_active_paths(
                grb_env, path_edge, path_demand, x, demand_vec, path_number);

            // Collect active path column indices (p > 0)
            vector<int> inds;
            for (int j = 0; j < path_number; ++j)
                if (p(j) > 0.0) inds.push_back(j);
            int n_i = (int)inds.size();

            // 5. Restrict to active paths
            //    path_edge_i   [n_links x n_i]
            //    path_demand_i [n_OD   x n_i]
            MatrixXd path_edge_i   = extract_columns(path_edge,   inds);
            MatrixXd path_demand_i = extract_columns(path_demand, inds);

            // 6. Sensitivity analysis (KKT system)
            //
            //    u_x = 0.15 * 4 * tfree .* x^3 ./ cap^4
            //        = 0.6 * tfree .* x^3 ./ cap^4
            VectorXd u_x = 0.6
                * tfree.cwiseProduct(x.array().pow(3).matrix())
                         .cwiseQuotient(cap.array().pow(4).matrix());

            //    c_f = path_edge_i.T * diag(u_x) * path_edge_i   [n_i x n_i]
            MatrixXd c_f = path_edge_i.transpose()
                         * u_x.asDiagonal()
                         * path_edge_i;

            //    J = [ c_f              -path_demand_i.T ]   [(n_i+n_OD) x (n_i+n_OD)]
            //        [ path_demand_i    0                ]
            int n_total = n_i + n_OD;
            MatrixXd J(n_total, n_total);
            J.topLeftCorner(n_i, n_i)       = c_f;
            J.topRightCorner(n_i, n_OD)     = -path_demand_i.transpose();
            J.bottomLeftCorner(n_OD, n_i)   = path_demand_i;
            J.bottomRightCorner(n_OD, n_OD).setZero();

            //    c_toll = path_edge_i.T    [n_i x n_links]
            //    Analytical Jacobian of path costs w.r.t. toll vector:
            //      c_i = path_edge_i.T @ (f(x) + toll)  is linear in toll,
            //      so dc_i/d_toll = path_edge_i.T.
            //    Mirrors: c_toll = jacobian(func_c_i, toll.detach()) in Python.
            MatrixXd c_toll = path_edge_i.transpose();  // [n_i x n_links]

            //    Right = [ -c_toll              ]   [(n_i+n_OD) x n_links]
            //            [ zeros(n_OD, n_links) ]
            MatrixXd Right(n_total, n_links);
            Right.topRows(n_i).setZero();
            Right.topRows(n_i)    = -c_toll;
            Right.bottomRows(n_OD).setZero();

            //    Solve:  J * S = Right   ->  S [(n_i+n_OD) x n_links]
            MatrixXd S      = J.colPivHouseholderQr().solve(Right);
            MatrixXd f_toll = S.topRows(n_i);        // [n_i x n_links]

            //    x_toll = path_edge_i * f_toll   [n_links x n_links]
            //    Sensitivity of link flows to toll vector.
            MatrixXd x_toll = path_edge_i * f_toll;

            // 7. Gradient of TT w.r.t. toll  (manual backprop of objective.backward())
            //
            //    Python uses stop-gradient trick:
            //      x_vir = x_toll @ toll - stopgrad(x_toll @ toll - x)
            //    so x_vir == x numerically, but gradient flows through x_toll @ toll.
            //
            //    TT = dot(x_vir, t_vir)  with t_vir = tfree*(1+0.15*(x_vir/cap)^4)
            //    Chain rule gives:
            //      d(TT)/d(x_vir) = t_vir + x_vir .* u_x    (product rule)
            //      d(TT)/d(toll)  = x_toll.T @ (t_vir + x_vir .* u_x)
            //    (u_x evaluated at x_vir ≈ x)
            VectorXd ratio = x.cwiseQuotient(cap);
            VectorXd t_vir = tfree.cwiseProduct(
                VectorXd::Ones(n_links) + 0.15 * ratio.array().pow(4).matrix());
            VectorXd grad_x    = t_vir + x.cwiseProduct(u_x);   // [n_links]
            VectorXd grad_toll = x_toll.transpose() * grad_x;   // [n_links]

            auto toc2 = chrono::high_resolution_clock::now();
            inverting_time += chrono::duration<double>(toc2 - tic2).count();

            // 8. Toll update:  toll -= alpha * beta / (iter_num + beta) * grad_toll
            double step = params.alpha * params.beta / (iter_num + params.beta);
            toll -= step * grad_toll;

            // Enforce non-negativity (clamp) and zero out non-toll links
            toll = toll.cwiseMax(0.0);
            for (int idx : params.non_link) toll(idx) = 0.0;

            // 9. Convergence check: L-inf norm of toll change
            double change = (toll - toll_old).lpNorm<Infinity>();
            if (change < params.xi) break;

            toll_old = toll;
            iter_num++;
        }

        auto toc    = chrono::high_resolution_clock::now();
        double time = chrono::duration<double>(toc - tic).count();

        // ---- Solution evaluation ----
        // Re-solve UE with tighter tolerance, then compute normalised objective.
        {
            map<Link, double> Toll;
            for (int kk = 0; kk < n_links; ++kk)
                Toll[net.Link_list[kk]] = toll(kk);

            GP(net, params.Demand, params.Fftime, params.Capacity,
               1000, params.eps_eva, /*warm_start=*/true, &Toll);

            VectorXd x_eval(n_links);
            for (int i = 0; i < n_links; ++i)
                x_eval(i) = net.flow.at(net.Link_list[i]);

            double TT  = total_travel_time(x_eval, tfree, cap);
            double obj = (TT - TT_so) / (TT_ue - TT_so);

            cout << "objective: " << obj << endl;

            // ---- Store / print results ----
            // In Python, results are pickled to result_sctp/sab{n_ini}_{NETWORK}.pkl
            // Here we print; serialisation hooks can be added as needed.
            cout << "[init " << n_ini << "] "
                 << "obj=" << obj
                 << "  time=" << time << "s"
                 << "  inverting_time=" << inverting_time << "s"
                 << endl;

            cout << "toll at toll_link:";
            for (int idx : params.toll_link) cout << " " << toll(idx);
            cout << endl;
        }
    }
}

// ============================================================
// main  (mirrors the top-level script in sctp_sab.py)
//
// Network data (topology, demands, UE/SO flows, initialisations) must be
// loaded here.  In the Python code these come from pickle files written by
// sctp_pre.py; replace the TODO sections with your preferred I/O method.
// ============================================================
int main() {

    for (const string& NETWORK : {"hearn1", "hearn2"}) {

        // ---- Build / load Network object ----
        // Python: net = pickle.load(open('result_sctp/net_' + NETWORK + '.pkl'))
        Network net;
        net.name = NETWORK;

        // TODO: populate net.Innode, net.Outnode, net.Link_list, net.Odtree
        //       from the network text files (see read.py for file format)

        // ---- Build / load NetworkParams ----
        // Python: instance = Instance(NETWORK)  (see sctp_load.py)
        NetworkParams params;

        // TODO: populate params.Demand, params.Fftime, params.Capacity,
        //       params.demand, params.cap, params.tfree,
        //       params.toll_link, params.non_link, params.n_links, params.n_OD
        //
        // Example for hearn1:
        //   params.toll_link = {1, 10, 11};
        //   params.non_link  = all indices not in toll_link;
        //   params.epsilon   = 1e-4;
        //   params.xi        = 1e-4;
        //   params.eps_eva   = 1e-8;

        if (NETWORK == "hearn1" || NETWORK == "hearn2") {
            params.alpha   = 0.05;
            params.beta    = 20.0;
            params.epsilon = 1e-4;
        }

        // ---- Load UE and SO link flows ----
        // Python: x_ue = torch.load('result_sctp/ue_' + NETWORK + '.pt')
        //         x_so = torch.load('result_sctp/so_' + NETWORK + '.pt')
        VectorXd x_ue;  // TODO: load from file
        VectorXd x_so;  // TODO: load from file

        // ---- Load initialisations ----
        // Python: initializations = pickle.load(open(...))
        // Each element is a VectorXd of initial toll values at toll_link indices.
        vector<VectorXd> initializations;  // TODO: load from file

        // ---- Run SAB ----
        run_sctp_sab(net, params, initializations, x_ue, x_so);
    }

    return 0;
}
