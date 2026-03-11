# CLAUDE.md — Codebase Guide for AI Assistants

## Project Overview

This repository implements **Differentiable Bilevel Programming** methods for solving transportation network optimization problems. The core contribution is making lower-level user equilibrium (UE) solvers differentiable via PyTorch autograd, enabling gradient-based optimization of bilevel problems at scale.

Two main problem types are studied:

- **CNDP** (Capacity Network Design Problem): Augment link capacities to minimize total travel time plus augmentation cost.
- **SCTP** (Spatial Toll Collection Problem): Set link tolls to maximize revenue subject to user equilibrium.

Both problems share the bilevel structure:

```
Upper level (leader):     min_y  obj_leader(x(y), y)
Lower level (followers):  x(y) = argmin_x  obj_followers(x, y)
```

---

## Repository Structure

```
.
├── graph.py              # Core Network class: shortest paths, UE/SO solvers, path enumeration
├── read.py               # I/O utilities: parse .net, .trp, .nod, .tol, .tb files
├── GP.py                 # Gradient Projection user equilibrium solver
├── GP_SO.py              # System Optimal solver (marginal cost BPR)
│
├── cndp_pre.py           # Preprocessing: run once to enumerate paths, saves .pkl
├── cndp_load.py          # Instance class: loads/configures a CNDP instance
├── cndp_setting.py       # Loads existing results, selects network/parameters
├── cndp_so.py            # System Optimal baseline for CNDP
├── cndp_da.py            # Dual Annealing baseline for CNDP
├── cndp_sa.py            # Simulated Annealing baseline for CNDP
├── cndp_braess.py        # Standalone Braess paradox experiment (executable)
├── cndp_dolmd.py         # DolMD algorithm for CNDP
├── cndp_sab.py           # SAB algorithm for CNDP
├── cndp_silmd.py         # SilMD algorithm for CNDP
│
├── sctp_pre.py           # Preprocessing for SCTP
├── sctp_init.py          # Initialize toll solutions
├── sctp_load.py          # Instance loader for SCTP
├── sctp_setting.py       # Configuration and result loading for SCTP
├── sctp_cp.py            # Cross-derivative Partial solver
├── sctp_ct.py            # Cost-gradient Toll solver
├── sctp_da.py            # Dual Annealing baseline for SCTP
├── sctp_mct.py           # Modified Cost-gradient Toll
├── sctp_emcdt.py         # Enhanced Modified CDT
├── sctp_dolmd.py         # DolMD algorithm for SCTP
├── sctp_sab.py           # SAB algorithm for SCTP
├── sctp_silmd.py         # SilMD algorithm for SCTP
├── sctp_hearn.py         # Hearn network experiments
├── sctp_hearn2.py        # Hearn network variant
│
├── network/              # Transportation network data
│   ├── bar/              # Barnett network
│   ├── cs/               # Corridor (Sioux Falls) network
│   ├── sf/               # San Francisco network
│   └── hearn/            # Hearn network
│
├── result_cndp/          # Serialized results from CNDP experiments (.pkl files)
├── result_sctp/          # Serialized results from SCTP experiments (.pkl files)
└── LICENSE               # MIT License
```

---

## Key Algorithms

### 1. GP — Gradient Projection (UE Solver)
**File:** `GP.py`

Classic iterative user equilibrium solver. Not differentiable. Used as ground-truth lower-level solver in baselines.

- `BPR(x, t0, c)` — Bureau of Public Roads congestion: `t = t0 * (1 + 0.15*(x/c)^4)`
- `GP(network, demand, ...)` — Main solver; alternates shortest-path assignment and flow rebalancing. Accepts optional tolls for SCTP.

### 2. DolMD — Differentiation On-Level MultiDimension
**Files:** `cndp_dolmd.py`, `sctp_dolmd.py`

Solves lower-level equilibrium approximately using a **logit-based path choice model**. Path probabilities update as `p *= exp(-r * c)` and are differentiable w.r.t. upper-level variables. Gradient of upper-level objective flows through PyTorch autograd.

### 3. SAB — Semi-Analytic Bilevel
**Files:** `cndp_sab.py`, `sctp_sab.py`

Computes **implicit derivatives** via Jacobian inversion of equilibrium conditions. Uses `torch.autograd.functional.jacobian`. Most mathematically rigorous but computationally expensive for large networks.

### 4. SilMD — Single-Internal-Level MultiDimension
**Files:** `cndp_silmd.py`, `sctp_silmd.py`

Runs `T` inner iterations of equilibrium solver then differentiates through those `T` steps. Interpolates between DolMD (T=1, fast/approximate) and SAB (T→∞, exact).

### 5. Baselines
- **IOA** — Iterative Optimization Assignment: alternates optimizing upper-level and solving lower-level exactly. No gradients.
- **Dual Annealing** (`cndp_da.py`, `sctp_da.py`) — scipy `dual_annealing` black-box baseline.
- **Simulated Annealing** (`cndp_sa.py`) — scipy `dual_annealing` with SA settings.
- **SO** — System Optimal (lower bound for CNDP, not a true bilevel solution).

---

## Network Model

### BPR Function
```python
t(x) = t0 * (1 + alpha * (x / c) ** beta)
# Standard parameters: alpha=0.15, beta=4
```

### Network File Formats
- `.net` — Link data: from_node, to_node, capacity, length, free-flow-time, ...
- `.trp` — OD demand matrices
- `.nod` — Node coordinates
- `.tol` — Toll values per link
- `.tb` — Toll bounds

All parsed by `read.py` helpers: `read_net()`, `read_trp()`, `read_nod()`, `read_tol()`, `read_tb()`.

---

## Dependencies

No `requirements.txt` exists. Install these packages:

```
torch          # Differentiable computation, autograd, sparse tensors
numpy          # Array operations
scipy          # Optimization (minimize, dual_annealing, root), linear algebra (qr)
numba          # JIT compilation for tight inner loops
networkx       # Graph utilities (used in sctp_setting.py)
pandas         # Data frames (used in GP_SO.py)
matplotlib     # Plotting (used in cndp_braess.py)
```

Example install:
```bash
pip install torch numpy scipy numba networkx pandas matplotlib
```

---

## Development Workflow

### 1. Preprocessing (run once per network)

Before running any algorithm, preprocess the network to enumerate paths:

```bash
python cndp_pre.py    # Saves preprocessed CNDP network to result_cndp/*.pkl
python sctp_pre.py    # Saves preprocessed SCTP network to result_sctp/*.pkl
```

### 2. Running Experiments

Each algorithm file is a standalone script. Select the network and parameters inside the file via the configuration variables at the top (e.g., `NETWORK`, `alpha`, `r`), then run:

```bash
# CNDP experiments
python cndp_braess.py    # Braess paradox study (self-contained, no preprocessing needed)
python cndp_dolmd.py     # Run DolMD on CNDP
python cndp_sab.py       # Run SAB on CNDP
python cndp_silmd.py     # Run SilMD on CNDP
python cndp_da.py        # Run Dual Annealing baseline

# SCTP experiments
python sctp_dolmd.py     # Run DolMD on SCTP
python sctp_sab.py       # Run SAB on SCTP
python sctp_silmd.py     # Run SilMD on SCTP
python sctp_da.py        # Run Dual Annealing baseline
```

### 3. Loading Results

Results are serialized to pickle files in `result_cndp/` or `result_sctp/`. Use `cndp_setting.py` / `sctp_setting.py` to load and compare stored results.

---

## Test Networks

| Name | Code | Nodes | Description |
|------|------|-------|-------------|
| San Francisco | `sf` | ~76 | Medium-scale urban network |
| Barnett | `bar` | ~24 | Small benchmark network |
| Corridor | `cs` | ~13 | Small corridor network |
| Hearn | `hearn` | ~4 | Toy network for SCTP |
| Braess | (inline) | 4 | Toy network in `cndp_braess.py` |

---

## Key Classes and Functions

### `Network` (graph.py)
Central data structure for a transportation network.

```python
net = Network(...)
net.LC(od)              # List-based shortest path
net.LS(od)              # Dijkstra shortest path
net.LC_2factors(od)     # Two-criteria shortest path
net.solve_ue()          # Solve user equilibrium
net.solve_so()          # Solve system optimal
```

### `Instance` (cndp_load.py / sctp_load.py)
Problem instance loader. Configures network, enhancement links, scaling parameters.

```python
inst = Instance('sf')   # Load San Francisco CNDP instance
```

### `GP(network, demand, ...)` (GP.py)
Returns: `(path_flows, link_flows, travel_times, gap_history)`

### `GP_SO(network, demand, ...)` (GP_SO.py)
Returns system-optimal flows.

---

## Code Conventions

- **No formal testing framework** — experiments are executable scripts; correctness is verified by comparing against known baselines (SO lower bound, GP exact UE).
- **PyTorch tensors** are used throughout for differentiable operations. CPU/GPU device is typically selected at the top of each script.
- **Pickle serialization** — intermediate results (paths, equilibrium solutions) are cached as `.pkl` files to avoid recomputation.
- **Configuration at top of script** — each experiment file has a block at the top (network name, algorithm parameters, random seeds) that serves as configuration.
- **Sparse tensors** — path-link incidence matrices and flow propagation use `torch.sparse` for memory efficiency on large networks.
- **Numba JIT** — tight loops (e.g., in GP projection steps) use `@numba.jit` for performance.
- **Result naming convention** — output files follow `result_{problem}/{network}_{algorithm}_{params}.pkl`.

---

## Mathematical Background

### CNDP Objective
```
min_{y}  TT(x(y), y) + eta * sum_i(cash_i * delta_c_i^r)
```
- `TT`: Total travel time = `sum(flow * time)`
- `delta_c`: Capacity augmentation on each link
- `cash`: Unit augmentation cost
- `eta`, `r`: Scaling parameters

### SCTP Objective
```
max_{tau}  sum_i(tau_i * x_i(tau))
```
- `tau`: Toll vector (upper-level decision)
- `x(tau)`: Link flows at UE with tolls

### Logit Path Choice (DolMD)
```
p_k  <-  p_k * exp(-r * c_k) / Z
```
- `c_k`: Cost of path k
- `r`: Smoothing parameter (larger = sharper, closer to shortest-path)
- `Z`: Normalization

---

## Common Pitfalls

1. **Skip preprocessing**: Always run `cndp_pre.py` / `sctp_pre.py` before algorithm scripts; they load the serialized network.
2. **Device mismatch**: Ensure all tensors are on the same device (CPU or CUDA). Check `device` variable at top of each script.
3. **Large `r` parameter**: Very large `r` in DolMD can cause numerical overflow in `exp(-r * c)`. Use log-sum-exp tricks or moderate `r`.
4. **Memory on large networks**: San Francisco network with full path enumeration can be memory-intensive. SAB's Jacobian computation scales poorly; prefer DolMD or SilMD for large instances.
5. **Pickle version compatibility**: `.pkl` files saved with one Python/PyTorch version may not load with another. Re-run preprocessing if loading fails.

---

## License

MIT License — Copyright 2024 Jiayang Li
