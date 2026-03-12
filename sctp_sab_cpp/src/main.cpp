#include "types.h"
#include "read.h"
#include "network.h"
#include "instance.h"
#include "gp.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <gurobi_c++.h>

#include <iostream>
#include <chrono>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>

// Compute total travel time: TT = sum_i x_i * t_i where t_i = tfree_i * (1 + 0.15*(x_i/c_i)^4)
static double total_travel_time(const Eigen::VectorXd& x,
                                const Eigen::VectorXd& tfree,
                                const Eigen::VectorXd& cap) {
    int n = (int)x.size();
    double tt = 0.0;
    for (int i = 0; i < n; ++i) {
        double ratio = x(i) / cap(i);
        double r4 = ratio * ratio * ratio * ratio;
        double t_i = tfree(i) * (1.0 + 0.15 * r4);
        tt += x(i) * t_i;
    }
    return tt;
}

// Find linearly independent columns using Gurobi LP.
// Directly iterates sparse matrix non-zeros — no dense conversion needed.
// Solves: min 0^T f  s.t. [path_edge; path_demand] * f = [x; demand], f >= 0
// Returns boolean mask where f > 0 (basic variables = linearly independent columns).
static std::vector<bool> find_independent_paths_gurobi(
        const Eigen::SparseMatrix<double>& path_edge,
        const Eigen::SparseMatrix<double>& path_demand,
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& demand_vec,
        int num_paths) {

    GRBEnv env(true);
    env.set(GRB_IntParam_OutputFlag, 0);
    env.start();

    GRBModel model(env);

    int num_links = (int)path_edge.rows();
    int num_od = (int)path_demand.rows();
    int num_constraints = num_links + num_od;

    // Add variables f_k >= 0
    std::vector<GRBVar> f(num_paths);
    for (int k = 0; k < num_paths; ++k) {
        f[k] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
    }

    // Build constraint expressions by iterating sparse non-zeros
    std::vector<GRBLinExpr> constrs(num_constraints);
    for (int i = 0; i < num_constraints; ++i) {
        constrs[i] = 0;
    }

    // path_edge contributions (rows 0..num_links-1)
    // path_edge is column-major (Eigen default), iterate over columns
    for (int k = 0; k < path_edge.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(path_edge, k); it; ++it) {
            constrs[it.row()] += it.value() * f[it.col()];
        }
    }

    // path_demand contributions (rows num_links..num_links+num_od-1)
    for (int k = 0; k < path_demand.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(path_demand, k); it; ++it) {
            constrs[num_links + it.row()] += it.value() * f[it.col()];
        }
    }

    // Add equality constraints with RHS
    for (int i = 0; i < num_links; ++i) {
        model.addConstr(constrs[i] == x(i));
    }
    for (int i = 0; i < num_od; ++i) {
        model.addConstr(constrs[num_links + i] == demand_vec(i));
    }

    // Objective: min 0 (feasibility only)
    model.setObjective(GRBLinExpr(0), GRB_MINIMIZE);
    model.optimize();

    // Extract solution: f > 0 are independent paths
    std::vector<bool> inds(num_paths, false);
    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        for (int k = 0; k < num_paths; ++k) {
            if (f[k].get(GRB_DoubleAttr_X) > 1e-10) {
                inds[k] = true;
            }
        }
    } else {
        std::cerr << "Warning: Gurobi LP did not find optimal solution. Status: "
                  << model.get(GRB_IntAttr_Status) << std::endl;
    }

    return inds;
}

// Extract columns from sparse matrix by boolean mask, returning a sparse matrix
static Eigen::SparseMatrix<double> extract_columns_sparse(
        const Eigen::SparseMatrix<double>& sp,
        const std::vector<bool>& mask) {

    // Build column index mapping
    int n_cols = 0;
    std::vector<int> col_map(mask.size(), -1);
    for (int k = 0; k < (int)mask.size(); ++k) {
        if (mask[k]) {
            col_map[k] = n_cols++;
        }
    }

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(sp.nonZeros());  // upper bound

    for (int k = 0; k < sp.outerSize(); ++k) {
        if (!mask[k]) continue;
        int new_col = col_map[k];
        for (Eigen::SparseMatrix<double>::InnerIterator it(sp, k); it; ++it) {
            triplets.emplace_back(it.row(), new_col, it.value());
        }
    }

    Eigen::SparseMatrix<double> result(sp.rows(), n_cols);
    result.setFromTriplets(triplets.begin(), triplets.end());
    result.makeCompressed();
    return result;
}

int main(int argc, char* argv[]) {
    // Base directory (parent of sctp_sab_cpp/)
    std::string base_dir = "..";
    if (argc > 2) {
        base_dir = argv[2];
    }

    // Networks to process
    std::vector<std::string> networks = {"hearn1", "hearn2"};
    if (argc > 1) {
        networks = {argv[1]};
    }

    for (const auto& NETWORK : networks) {
        std::cout << "=== Processing network: " << NETWORK << " ===" << std::endl;

        // Load instance
        Instance instance(NETWORK, base_dir);

        const auto& Capacity = instance.Capacity;
        const auto& Fftime = instance.Fftime;
        const auto& Demand = instance.Demand;
        const auto& Link_list = instance.Link_list;

        const Eigen::VectorXd& cap = instance.cap;
        const Eigen::VectorXd& tfree = instance.tfree;
        const std::vector<int>& toll_link = instance.toll_link;
        double epsilon = instance.epsilon;
        double xi = instance.xi;
        double eps_eva = instance.eps_eva;

        int Numberoflink = (int)Link_list.size();

        // Build demand vector (ordered by Demand map)
        Eigen::VectorXd demand_vec = instance.demand_vec;

        // Non-toll links
        std::vector<int> non_link;
        std::set<int> toll_set(toll_link.begin(), toll_link.end());
        for (int i = 0; i < Numberoflink; ++i) {
            if (toll_set.find(i) == toll_set.end()) {
                non_link.push_back(i);
            }
        }

        // Hyperparameters
        double alpha, beta;
        if (NETWORK == "hearn1" || NETWORK == "hearn2") {
            alpha = 0.05;
            beta = 20.0;
            epsilon = 1e-4;
        } else if (NETWORK == "sf") {
            alpha = 0.001;
            beta = 100.0;
        } else if (NETWORK == "bar") {
            alpha = 0.02;
            beta = 10.0;
        } else if (NETWORK == "cs") {
            alpha = 0.05;
            beta = 20.0;
        } else {
            alpha = 0.05;
            beta = 20.0;
        }

        // Build Network and compute reference UE/SO
        Network net_ue(NETWORK, instance.Innode, instance.Outnode, Link_list, instance.Odtree);
        net_ue.solve_ue(Demand, Fftime, Capacity, 1000, 1e-8);

        Eigen::VectorXd x_ue(Numberoflink);
        for (int i = 0; i < Numberoflink; ++i) {
            x_ue(i) = net_ue.flow[Link_list[i]];
        }
        double TT_ue = total_travel_time(x_ue, tfree, cap);

        Network net_so(NETWORK, instance.Innode, instance.Outnode, Link_list, instance.Odtree);
        net_so.solve_so(Demand, Fftime, Capacity, 1000, 1e-8);

        Eigen::VectorXd x_so(Numberoflink);
        for (int i = 0; i < Numberoflink; ++i) {
            x_so(i) = net_so.flow[Link_list[i]];
        }
        double TT_so = total_travel_time(x_so, tfree, cap);

        std::cout << "TT_ue = " << TT_ue << ", TT_so = " << TT_so << std::endl;

        // --- Main SAB loop ---
        // Initialize toll (use zero initialization for simplicity)
        Eigen::VectorXd toll = Eigen::VectorXd::Zero(Numberoflink);
        Eigen::VectorXd toll_old = toll;

        // Build working network
        Network net(NETWORK, instance.Innode, instance.Outnode, Link_list, instance.Odtree);
        // Cold-start UE first
        net.solve_ue(Demand, Fftime, Capacity, 1000, epsilon);

        auto tic = std::chrono::high_resolution_clock::now();
        double inverting_time = 0.0;
        int iter_num = 0;

        while (true) {
            // Convert toll vector to dict
            LinkMapD Toll;
            for (int k = 0; k < Numberoflink; ++k) {
                Toll[Link_list[k]] = toll(k);
            }

            // Solve UE with current toll
            net.solve_ue(Demand, Fftime, Capacity, 1000, epsilon, true, &Toll);

            // Path enumeration and sparse matrix generation
            net.path_enumeration();
            net.generate_sparse_matrix();

            int path_number = net.path_number;

            auto tic2 = std::chrono::high_resolution_clock::now();

            // Build x (link flow vector)
            Eigen::VectorXd x(Numberoflink);
            for (int i = 0; i < Numberoflink; ++i) {
                x(i) = net.flow[Link_list[i]];
            }

            // ============================================================
            // Find linearly independent paths via Gurobi LP
            // Directly uses sparse matrices — no dense conversion
            // ============================================================
            auto inds = find_independent_paths_gurobi(
                net.path_edge, net.path_demand, x, demand_vec, path_number);
            int n_i = 0;
            for (bool b : inds) if (b) n_i++;

            // ============================================================
            // Extract independent path columns as sparse matrices
            // ============================================================
            Eigen::SparseMatrix<double> path_edge_i = extract_columns_sparse(net.path_edge, inds);
            Eigen::SparseMatrix<double> path_demand_i = extract_columns_sparse(net.path_demand, inds);

            // ============================================================
            // Sensitivity analysis (all sparse where possible)
            // ============================================================

            // u_x = dt/dx = 0.15 * tfree * 4 * (x/cap)^3 / cap
            Eigen::VectorXd u_x(Numberoflink);
            for (int i = 0; i < Numberoflink; ++i) {
                double ratio = x(i) / cap(i);
                u_x(i) = 0.15 * tfree(i) * 4.0 * ratio * ratio * ratio / cap(i);
            }

            // c_f = path_edge_i^T * diag(u_x) * path_edge_i   (n_i x n_i, dense OK)
            // Build sparse diagonal matrix from u_x
            Eigen::SparseMatrix<double> diag_ux(Numberoflink, Numberoflink);
            {
                std::vector<Eigen::Triplet<double>> diag_triplets;
                diag_triplets.reserve(Numberoflink);
                for (int i = 0; i < Numberoflink; ++i) {
                    if (std::abs(u_x(i)) > 1e-30) {
                        diag_triplets.emplace_back(i, i, u_x(i));
                    }
                }
                diag_ux.setFromTriplets(diag_triplets.begin(), diag_triplets.end());
            }

            // Sparse multiply: path_edge_i^T * diag(u_x) * path_edge_i → dense n_i x n_i
            Eigen::SparseMatrix<double> temp_sp = path_edge_i.transpose() * diag_ux * path_edge_i;
            Eigen::MatrixXd c_f = Eigen::MatrixXd(temp_sp);

            int num_od = (int)demand_vec.size();

            // Build KKT Jacobian J  (dense, size (n_i + num_od)^2, moderate)
            // J = [c_f, -path_demand_i^T; path_demand_i, 0]
            int J_size = n_i + num_od;
            Eigen::MatrixXd J = Eigen::MatrixXd::Zero(J_size, J_size);
            J.topLeftCorner(n_i, n_i) = c_f;

            // Fill path_demand_i into J (sparse → dense block, size num_od x n_i, acceptable)
            Eigen::MatrixXd pd_i_dense = Eigen::MatrixXd(path_demand_i);
            J.topRightCorner(n_i, num_od) = -pd_i_dense.transpose();
            J.bottomLeftCorner(num_od, n_i) = pd_i_dense;

            // ============================================================
            // Adjoint method: compute gradient with a SINGLE linear solve
            // instead of |links| solves.
            //
            // Original: S = J^{-1} * Right  (Right has |links| columns)
            //           x_toll = path_edge_i * S[:n_i,:]  (|links| x |links|)
            //           grad = x_toll^T * dTT/dx
            //
            // Adjoint:  v = path_edge_i^T * dTT/dx         (sparse x vec)
            //           Solve J^T * y = [v; 0]             (1 solve!)
            //           grad = -path_edge_i * y[:n_i]       (sparse x vec)
            //
            // Proof:  grad = x_toll^T * dTT/dx
            //       = (Δ_I · f_toll)^T · dTT/dx
            //       = f_toll^T · (Δ_I^T · dTT/dx)
            //       = S[:n_i,:]^T · v
            //       = (J^{-1} · Right)[:n_i,:]^T · v
            //     Let w = [v; 0], then S^T · w = grad (for all columns)
            //     S^T = Right^T · J^{-T}
            //     grad = Right^T · J^{-T} · w = Right^T · y  where J^T y = w
            //          = [-c_toll^T, 0] · y = -Δ_I · y[:n_i]
            // ============================================================

            // dTT/dx_i = t_i(x_i) + x_i * dt_i/dx_i = tfree_i * (1 + 0.75*(x_i/cap_i)^4)
            Eigen::VectorXd dTT_dx(Numberoflink);
            for (int i = 0; i < Numberoflink; ++i) {
                double ratio = x(i) / cap(i);
                double r4 = ratio * ratio * ratio * ratio;
                dTT_dx(i) = tfree(i) * (1.0 + 0.75 * r4);
            }

            // v = path_edge_i^T * dTT_dx   (sparse x vector → n_i vector)
            Eigen::VectorXd v = Eigen::VectorXd(path_edge_i.transpose() * dTT_dx);

            // w = [v; 0]
            Eigen::VectorXd w = Eigen::VectorXd::Zero(J_size);
            w.head(n_i) = v;

            // Solve J^T * y = w   (single linear solve, size n_i + num_od)
            Eigen::VectorXd y = J.transpose().partialPivLu().solve(w);

            // grad = -path_edge_i * y[:n_i]   (sparse x vector → |links| vector)
            Eigen::VectorXd grad_toll = -(Eigen::VectorXd(path_edge_i * y.head(n_i)));

            auto toc2 = std::chrono::high_resolution_clock::now();
            inverting_time += std::chrono::duration<double>(toc2 - tic2).count();

            // Update toll with projected gradient descent
            double lr = alpha * beta / (iter_num + beta);
            toll -= lr * grad_toll;

            // Project: clamp to non-negative
            toll = toll.cwiseMax(0.0);

            // Zero out non-toll links
            for (int i : non_link) {
                toll(i) = 0.0;
            }

            // Check convergence
            double change = (toll - toll_old).lpNorm<Eigen::Infinity>();

            if (change < xi) break;

            toll_old = toll;
            iter_num++;
        }

        auto toc = std::chrono::high_resolution_clock::now();
        double total_time = std::chrono::duration<double>(toc - tic).count();

        // Final evaluation: solve UE with high precision
        LinkMapD Toll_final;
        for (int k = 0; k < Numberoflink; ++k) {
            Toll_final[Link_list[k]] = toll(k);
        }
        net.solve_ue(Demand, Fftime, Capacity, 1000, eps_eva, true, &Toll_final);

        Eigen::VectorXd x_final(Numberoflink);
        for (int i = 0; i < Numberoflink; ++i) {
            x_final(i) = net.flow[Link_list[i]];
        }
        double TT_final = total_travel_time(x_final, tfree, cap);
        double obj = (TT_final - TT_so) / (TT_ue - TT_so);

        std::cout << "objective: " << obj << std::endl;
        std::cout << "iterations: " << iter_num << std::endl;
        std::cout << "total time: " << total_time << "s" << std::endl;
        std::cout << "inverting time: " << inverting_time << "s" << std::endl;

        // Print toll values on toll links
        std::cout << "toll values on toll links:" << std::endl;
        for (int idx : toll_link) {
            std::cout << "  link " << idx << " (" << Link_list[idx].first
                      << "," << Link_list[idx].second << "): " << toll(idx) << std::endl;
        }
        std::cout << std::endl;
    }

    return 0;
}
