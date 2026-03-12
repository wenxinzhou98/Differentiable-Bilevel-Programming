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

// Find linearly independent columns using Gurobi LP
// Solve: min 0^T f  s.t. A_eq * f = b_eq, f >= 0
// Return indices where f > 0 (basic variables = linearly independent columns)
static std::vector<bool> find_independent_paths_gurobi(
        const Eigen::MatrixXd& A_eq,
        const Eigen::VectorXd& b_eq,
        int num_paths) {

    GRBEnv env(true);
    env.set(GRB_IntParam_OutputFlag, 0);  // Suppress output
    env.start();

    GRBModel model(env);

    // Add variables f_k >= 0
    std::vector<GRBVar> f(num_paths);
    for (int k = 0; k < num_paths; ++k) {
        f[k] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "f_" + std::to_string(k));
    }

    // Add equality constraints: A_eq * f = b_eq
    int num_constraints = (int)A_eq.rows();
    for (int i = 0; i < num_constraints; ++i) {
        GRBLinExpr expr = 0;
        for (int k = 0; k < num_paths; ++k) {
            if (std::abs(A_eq(i, k)) > 1e-15) {
                expr += A_eq(i, k) * f[k];
            }
        }
        model.addConstr(expr == b_eq(i), "eq_" + std::to_string(i));
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

// Extract columns from sparse matrix by boolean mask
static Eigen::MatrixXd extract_columns(const Eigen::SparseMatrix<double>& sp, const std::vector<bool>& mask) {
    int n_cols = 0;
    for (bool b : mask) if (b) n_cols++;

    Eigen::MatrixXd result(sp.rows(), n_cols);
    int col = 0;
    for (int k = 0; k < (int)mask.size(); ++k) {
        if (mask[k]) {
            result.col(col) = Eigen::VectorXd(sp.col(k));
            col++;
        }
    }
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

            // Build x (link flow vector) and demand vector
            Eigen::VectorXd x(Numberoflink);
            for (int i = 0; i < Numberoflink; ++i) {
                x(i) = net.flow[Link_list[i]];
            }

            // Convert sparse matrices to dense for LP and sensitivity analysis
            Eigen::MatrixXd path_edge_dense = Eigen::MatrixXd(net.path_edge);
            Eigen::MatrixXd path_demand_dense = Eigen::MatrixXd(net.path_demand);

            // Build A_eq = [path_edge; path_demand] and b_eq = [x; demand]
            int num_rows = Numberoflink + (int)demand_vec.size();
            Eigen::MatrixXd A_eq(num_rows, path_number);
            A_eq.topRows(Numberoflink) = path_edge_dense;
            A_eq.bottomRows(demand_vec.size()) = path_demand_dense;

            Eigen::VectorXd b_eq(num_rows);
            b_eq.head(Numberoflink) = x;
            b_eq.tail(demand_vec.size()) = demand_vec;

            // Find linearly independent paths via Gurobi LP
            auto inds = find_independent_paths_gurobi(A_eq, b_eq, path_number);
            int n_i = 0;
            for (bool b : inds) if (b) n_i++;

            // Extract independent path columns
            Eigen::MatrixXd path_edge_i = extract_columns(net.path_edge, inds);      // |links| x n_i
            Eigen::MatrixXd path_demand_i = extract_columns(net.path_demand, inds);   // |ODs| x n_i

            // Sensitivity analysis
            // u_x = 0.15 * tfree * 4 * (x / cap)^3 / cap = 0.6 * tfree * x^3 / cap^4
            Eigen::VectorXd u_x(Numberoflink);
            for (int i = 0; i < Numberoflink; ++i) {
                double ratio = x(i) / cap(i);
                u_x(i) = 0.15 * tfree(i) * 4.0 * ratio * ratio * ratio / cap(i);
            }

            // c_f = path_edge_i^T * diag(u_x) * path_edge_i   (n_i x n_i)
            Eigen::MatrixXd diag_ux = u_x.asDiagonal();
            Eigen::MatrixXd c_f = path_edge_i.transpose() * diag_ux * path_edge_i;

            int num_od = (int)demand_vec.size();

            // Build KKT Jacobian J
            // J = [c_f, -path_demand_i^T; path_demand_i, 0]
            int J_size = n_i + num_od;
            Eigen::MatrixXd J = Eigen::MatrixXd::Zero(J_size, J_size);
            J.topLeftCorner(n_i, n_i) = c_f;
            J.topRightCorner(n_i, num_od) = -path_demand_i.transpose();
            J.bottomLeftCorner(num_od, n_i) = path_demand_i;
            // bottom-right is zeros

            // c_toll = jacobian of func_c_i w.r.t. toll
            // func_c_i(toll) = path_edge_i^T * (tfree*(1+0.15*(x/cap)^4) + toll)
            // dc_i/dtoll = path_edge_i^T * I = path_edge_i^T
            // This is an n_i x Numberoflink matrix
            Eigen::MatrixXd c_toll = path_edge_i.transpose();  // n_i x Numberoflink

            // Right-hand side: [-c_toll; 0]
            Eigen::MatrixXd Right = Eigen::MatrixXd::Zero(J_size, Numberoflink);
            Right.topRows(n_i) = -c_toll;

            // Solve: J * S = Right  =>  S = J^{-1} * Right
            Eigen::MatrixXd S = J.partialPivLu().solve(Right);

            // Extract path flow sensitivity: f_toll = S[:n_i, :]
            Eigen::MatrixXd f_toll = S.topRows(n_i);

            // Link flow sensitivity: x_toll = path_edge_i * f_toll   (Numberoflink x Numberoflink)
            Eigen::MatrixXd x_toll = path_edge_i * f_toll;

            auto toc2 = std::chrono::high_resolution_clock::now();
            inverting_time += std::chrono::duration<double>(toc2 - tic2).count();

            // Compute gradient analytically (replacing PyTorch autograd)
            // TT = sum_i x_vir_i * t_i(x_vir_i)
            // where x_vir has forward value = x, gradient dx_vir/dtoll = x_toll
            //
            // dTT/dx_i = t_i(x_i) + x_i * dt_i/dx_i
            //          = tfree_i * (1 + 0.15*(x_i/cap_i)^4) + x_i * 0.6*tfree_i*x_i^3/cap_i^4
            //          = tfree_i * (1 + 0.15*(x_i/cap_i)^4 + 0.6*(x_i/cap_i)^4)
            //          = tfree_i * (1 + 0.75*(x_i/cap_i)^4)
            //
            // dTT/dtoll = x_toll^T * dTT/dx
            Eigen::VectorXd dTT_dx(Numberoflink);
            for (int i = 0; i < Numberoflink; ++i) {
                double ratio = x(i) / cap(i);
                double r4 = ratio * ratio * ratio * ratio;
                dTT_dx(i) = tfree(i) * (1.0 + 0.75 * r4);
            }

            Eigen::VectorXd grad_toll = x_toll.transpose() * dTT_dx;

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
