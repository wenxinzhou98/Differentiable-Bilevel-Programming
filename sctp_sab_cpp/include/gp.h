#pragma once

#include "types.h"

// BPR travel time function: t = t0 * (1 + 0.15 * (x/c)^4)
inline double BPR(double fftime, double capacity, double flow) {
    double ratio = flow / capacity;
    return fftime * (1.0 + 0.15 * ratio * ratio * ratio * ratio);
}

// First derivative of BPR: dt/dx = 0.6 * t0 * x^3 / c^4
inline double BPR_1_derivative(double x_a, double c_a, double t0) {
    return 4.0 * t0 * 0.15 / (c_a * c_a * c_a * c_a) * (x_a * x_a * x_a);
}

// Second derivative of BPR (used in SO): d2t/dx2 = 1.8 * t0 * x^2 / c^4
inline double BPR_2_derivative(double x_a, double c_a, double t0) {
    return 3.0 * 4.0 * t0 * 0.15 / (c_a * c_a * c_a * c_a) * (x_a * x_a);
}

// Convert predecessor array to path
Path sppconvert(const std::map<int, int>& predecessor, int o, int d);

// Calculate path cost
double calc_pathcost(const Path& path, const LinkMapD& linkimpedance);

// Result of Frank-Wolfe algorithm
struct GPResult {
    PathFlowMap pathflow;
    LinkMapD flow;
    int numberofpath;
    PathSetMap pathset;
};

// Forward declaration of Network class
class Network;

// Frank-Wolfe Gradient Projection for User Equilibrium
GPResult GP(Network& net, const DemandMap& demand, const LinkMapD& fftime,
            const LinkMapD& capacity, int K0, double e0,
            bool warm_start = false, const LinkMapD* toll = nullptr,
            bool do_delete = true);

// Frank-Wolfe Gradient Projection for System Optimal
GPResult GP_SO(Network& net, const DemandMap& demand, const LinkMapD& fftime,
               const LinkMapD& capacity, int K0, double e0,
               bool warm_start = false);
