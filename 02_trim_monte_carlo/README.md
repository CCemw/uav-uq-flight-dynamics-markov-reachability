# 02_trim_monte_carlo: Trim sweep and P(trim) envelope

Produces the probability-of-trim envelope across forward airspeed. This is the central result the downstream Markov chain consumes.

## Contents

```
02_trim_monte_carlo/
├── run_nominal_trim_range.m     Deterministic trim sweep (nominal polars only)
├── run_mc_trim_range.m          Monte Carlo trim sweep with CRN
├── plot_trim_results.m          Generates the dissertation trim figures
└── results/
    ├── nominal_trim_results_QAV250.mat
    └── mc_trim_results_QAV250.mat
```

## Workflow

1. **Run `run_nominal_trim_range.m`** (~1 minute). This sweeps forward airspeed 0 to 30 m/s at 0.25 m/s resolution and saves the trim states, per-rotor thrust, torque, and power into `results/nominal_trim_results_QAV250.mat`.

2. **Run `run_mc_trim_range.m`** (multiple hours!). This runs 500 Monte Carlo samples at each of an airspeed grid. Each sample perturbs the blade polars and twist by a fresh draw from the σ envelope derived in Stage 1, using common random numbers (the same draw indexes across every airspeed for a given sample). The output is `results/mc_trim_results_QAV250.mat`.

3. **Run `plot_trim_results.m`** (seconds). This reads both `.mat` files and produces the three dissertation trim figures plus a V_max confidence table.

The committed `.mat` files are the results used in the dissertation. An examiner can go straight to the plotting script without re-running anything.

## How to adapt

### Airspeed range

Top of each run script, `V_sweep`. **IT MAY BE BETTER TO USE A LARGER STEP SIZE TO REDUCE RUNTIME!!**  >> RESULTS PRODUCED USE 0.25 M/S STEPS.

### Number of MC samples

Top of `run_mc_trim_range.m`, `N_MC` (default 500). For a tighter estimate of the P(trim) tails, increase it. Runtime scales linearly.

### RNG seed

`RNG_SEED` (default 42). Change for a different realisation of the CRN draws. The dissertation uses 42 for reproducibility.

### Twist uncertainty standard deviation

`sigma_twist` (default 1°). Adjust for a different manufacturing process.

### Convergence threshold

The MC classifier accepts residual norms below 10⁻³ (the relaxed threshold); the nominal script runs to 10⁻⁸ with damped Newton-Raphson. The relaxed MC threshold is discussed in dissertation Appendix D; changing it slightly affects the steep transition region of the P(trim) curve but not the qualitative shape.

### Integration method

The solver is damped Newton-Raphson with backtracking line search. The damping factor α = 0.6 and the three-halving backtrack sequence can be tuned inside `trim_solve`, but the values here are load-bearing and were tuned to converge the full MC.

## Solver diagnostics

Both run scripts print one line per trim solve: airspeed, trimmed pitch angle, rotor speeds. The MC script additionally prints a per-speed P(trim) and ETA.
