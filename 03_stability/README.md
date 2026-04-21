# 03_stability: Linearised stability analysis

Linearises the full six-degree-of-freedom dynamics about each trimmed state and computes the eigenvalues and stability derivatives, nominal and Monte Carlo.

## Contents

```
03_stability/
├── stability_analysis_mc.m              Trim + linearise + MC eigenvalue clouds
├── plot_stability_results.m             Generates the dissertation stability figures
└── results/
    └── stability_results_QAV250.mat
```

## Workflow

1. **Run `stability_analysis_mc.m`** (~45 minutes). This trims internally (no dependency on the trim folder) across V = 0 to 26 m/s at 1 m/s resolution, computes the 8×8 Jacobian about each trim, and extracts eigenvalues and four nominal stability derivatives (M_u, M_q, M_w, X_u). It then runs 500 MC samples at six target speeds (5, 10, 15, 17, 20, 24 m/s), storing the full eigenvalue cloud and per-sample M_u and M_q for each target.

2. **Run `plot_stability_results.m`** (seconds). Reads the `.mat` and produces four dissertation figures plus a console summary of mean ± σ values for max Re(λ), M_u, and M_q at each target speed.

## How to adapt

### Airspeed grid

Top of `stability_analysis_mc.m`, `V_sweep` (default `0:1:26`). The Jacobian is linearised about each trimmed point, so extending the range only produces useful results where trim is feasible.

### MC target speeds

`V_mc_targets` (default `[5, 10, 15, 17, 20, 24]`). The full 500-sample cloud is expensive; the target-speed approach gives a clear view of the eigenvalue spread at key operating points without running MC at every airspeed.

### Number of MC samples

`N_MC` (default 500). Same CRN seed as the trim stage (42) so the MC samples index the same physical population.

### Finite-difference step sizes for the Jacobian

Inside `compute_jacobian`, the per-state perturbation vector `delta` controls the accuracy of the central-difference derivative. The values here (`[0.1; 0.1; 0.1; 0.01; 0.01; 0.01; 0.005; 0.005]`) were verified using a plateau check. Change only if the equations of motion or state normalisation change.

### Extending to controllability

The current stage returns A matrices only. Adding a linearised controllability analysis (extension) with uncertainty propagation would require computing the B matrix (control Jacobian) alongside A, tracking the rank of [B, AB, A²B, ...] under MC variation. 
