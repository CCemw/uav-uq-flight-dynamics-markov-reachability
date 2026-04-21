# Uncertainty-Embedded Flight Dynamics Modelling and Markov Chain Reachability Analysis for a UAV

BEng Aerospace Engineering dissertation
University of Glasgow, 2026

**Author:** Edward Mills-Wierda
**Supervisors:** Dr Ye Yuan (first), Dr Chongfeng Wei (second)

---

## What this repository is

A MATLAB implementation of a general methodology for propagating blade-section manufacturing variability through the complete flight-dynamics pipeline of a small quadrotor UAV, demonstrated on the QAV250 Holybro UAV with APC 5×4.5" three-bladed propellers.

The pipeline runs in five stages. Manufacturing variation in the blade polars and twist is characterised from repeated wind-tunnel measurements. This resulting aerodynamic uncertainty is then propagated through a blade element momentum theory rotor model, a six-degree-of-freedom trim solver, a linearised stability analysis about each trimmed state, a closed-loop pitch manoeuvre, and finally a discrete-time Markov chain that rolls the trim credibility up into a per-flight-hour crash probability over a phased urban mission.

Everything is reproducible from the committed Monte Carlo data files. Every stage can also be re-run from scratch, and the scripts are self-contained except for the DTMC, which reads the trimmed Monte Carlo output.

The methodology is deliberately not specific to the QAV250. The platform parameters and blade polars are isolated in `00_inputs/`, and every stage script exposes its assumptions at the top of the file. Swap in a different rotorcraft, a different mission profile, or a different controller and the pipeline runs the same way.

---

## Repository layout
uav-uq-flight-dynamics-markov-reachability/
├── 00_inputs/                     Platform parameters and blade-section polars
├── 01_uncertainty_quantification/ Polar uncertainty characterisation (plot only)
├── 02_trim_monte_carlo/           Trim sweep: nominal + 500-sample Monte Carlo
├── 03_stability/                  Eigenvalue analysis, nominal + Monte Carlo
├── 04_closed_loop_manoeuvre/      PID-stabilised pulse response, nominal + Monte Carlo
├── 05_mission_dtmc/               Discrete-time Markov chain mission risk model
├── dissertation.pdf               Full written dissertation (methodology + results)
├── LICENSE                        MIT
└── README.md                      This file

Each folder has its own `README.md` describing its scripts, data files, and the settings a user would most likely want to change. The per-folder READMEs are the practical reference (for extensibility etc.), and this README guides the full structure and usage.

---

## How the five stages fit together

**Stage 0 — Inputs.** `00_inputs/platform/platform_parameters.m` defines mass, geometry, inertia, blade planform, and fuselage coefficients. `00_inputs/polars/afl{13,14,15}.c81.txt` are the three measured blade-section polars that together bracket the manufacturing variation. All downstream stages read from here.

**Stage 1 — Uncertainty quantification.** `01_uncertainty_quantification/plot_polar_uncertainty.m` reads the three polars and produces the mean ± σ envelope used throughout the rest of the pipeline. σ is derived from the range of the three measurements divided by the Tippett constant d₂.

**Stage 2 — Trim.** `02_trim_monte_carlo/` contains two runnable scripts. `run_nominal_trim_range.m` sweeps forward airspeed with the nominal polars only, producing the baseline trim envelope in about one minute. `run_mc_trim_range.m` runs the same sweep 500 times with common random numbers applied to the blade polars and twist. The resulting probability-of-trim envelope is the central input to the downstream Markov chain.

**Stage 3 — Stability.** `03_stability/stability_analysis_mc.m` linearises the six-degree-of-freedom dynamics about each trimmed state and computes the eigenvalues, nominal and Monte Carlo. 

**Stage 4 — Closed-loop manoeuvre.** `04_closed_loop_manoeuvre/closed_loop_manoeuvre_mc.m` simulates the time-domain pitch response to a rear-rotor RPM pulse under a cascaded PID controller. It is run four times: nominal and Monte Carlo, each at V = 5 and V = 15. (can be changed of course)

**Stage 5 — Mission DTMC.** `05_mission_dtmc/mission_dtmc.m` evolves a four-state chain (Normal Flight, Degraded Flight, Abort, Crash) through a user-defined mission. The per-phase transition probabilities are derived by convolving the P(trim) envelope from Stage 2 with a low-altitude Dryden gust distribution at the phase altitude. The final state distribution gives the probability of each mission outcome. Crash probability divided by mission duration gives the per-flight-hour crash rate which is compared against the JARUS SORA SAIL thresholds.

---

## How to run

### Prerequisites
- MATLAB R2024b or later
- Base MATLAB only. No additional toolboxes required. In particular, the DTMC is *purposefully* implemented directly by matrix multiplication rather than via MATLAB's `dtmc` class, to avoid a dependency on the Econometrics Toolbox. The phased structure of the mission, with a separate transition matrix per phase and an abort sub-chain entered between phases, also does not fit the single-chain assumption of the built-in class cleanly, so the raw matrix-multiplication implementation is both lighter and more natural. This can, of course, be adjusted.

### Running a stage
Navigate into the stage folder and run the script. For example:

```matlab
cd 02_trim_monte_carlo
run_nominal_trim_range
```

Each script writes its output to a `results/` subfolder within the same stage folder. Figures are written to a `figs/` subfolder, which is gitignored.

### Cross-folder dependencies
The only cross-folder read in the repository is the DTMC, which loads the trim Monte Carlo result from `02_trim_monte_carlo/results/`. All other stages re-trim internally when needed so they are independently runnable.

### Expected runtimes
| Script                               | Runtime     
-------------------------------------------------------------
| `run_nominal_trim_range`             | ~ seconds - minutes   
| `run_mc_trim_range`                  | ~ can be multiple hours (longest)
| `stability_analysis_mc`              | ~ can be multiple hours
| `closed_loop_manoeuvre_mc` (V=15 MC) | ~ can be a few hours
| `mission_dtmc`                       | ~ few seconds

All plotting scripts run in seconds and read from the committed `.mat` files, so anyone can regenerate every figure in the dissertation without re-running any Monte Carlo.

---

## Adjusting the model

Each stage folder's `README.md` lists the adjustments that can be made to each stage. In brief:

- **Platform.** Edit `00_inputs/platform/platform_parameters.m` to use a different airframe, mass, inertia, blade planform, or fuselage reference model.
- **Polar uncertainty.** Replace the three `.c81.txt` files in `00_inputs/polars/` with a different set. The current σ grid formulation (range over d₂) assumes exactly three samples, and hence the d₂ constant would need updating for different n.
- **Twist uncertainty.** Adjusted at the top of each Monte Carlo script (`sigma_twist` parameter).
- **Number of MC samples and RNG seed.** Top of each Monte Carlo script (`N_MC`, `RNG_SEED`).
- **Manoeuvre disturbance and PID gains.** Top of `closed_loop_manoeuvre_mc.m` — pulse magnitude, pulse duration, trim speed, and all four controller gains are exposed.
- **Mission profile.** Top of `mission_dtmc.m` — the mission is a cell array where each row is one phase. Add or remove rows, change phase durations, altitudes, or speeds.
- **Contingency routing.** Top of `mission_dtmc.m` — the emergency-landing success probability η, the degraded-state landing probability η_DF, and the severity-routing parameter s are exposed as scalars.
- **Atmospheric model.** Top of `mission_dtmc.m` — baseline wind, urban intensity scaling table, and correlation-length factor are editable.
- **SAIL thresholds.** The JARUS SORA table at the top of `mission_dtmc.m` can be replaced with any other regulatory threshold set.

## Extending the methodology/framework

Beyond the tunable parameters listed per stage, the repository is structured so a re-user can build on top of it. The pipeline is open to extensions. The following are examples: the stability stage could be augmented with linearised controllability analysis and uncertainty propagation; the DTMC could be replaced with a continuous-time Markov chain or a semi-Markov process to model phase-length variability; additional failure modes (battery, motor, sensor drop-out) could be added as further states without changing the trim-derived transition probability construction.

More broadly, the dissertation identifies the following directions for future work:

> Future work and extensions can strengthen the framework through improved empirical grounding, broader model fidelity, and wider application. Additional stand, wind-tunnel, flight-test, or operational data would allow the multiple stages of the methodology to be compared more directly against measured behaviour. In particular, the uncertainty model could be improved if the specimen dataset became larger, allowing assumed distributions to be replaced with ones that are grounded in data, and a Bayesian updating or calibration approach would be a natural next step, allowing prior uncertainty models to be refined as new evidence becomes available.
>
> In terms of the modelling scope, this could be extended through the inclusion of blade flapping, rotor-to-rotor interaction, rotor-to-body interaction, improved fuselage aerodynamics, and higher fidelity atmospheric or urban disturbance models. The same methodology could then be applied to different airframe configurations, broader operating environments, and alternative mission profiles, and in addition, future work could investigate controller design explicitly for robustness to the propagated uncertainty and assess the resulting disturbance rejection and mission-level outcomes. Finally, the present aerodynamic uncertainty contribution demonstrated here could also be further integrated with contributors such as propulsion faults, sensor failures, navigation uncertainty, traffic encounters, and contingency-management data, allowing the DTMC structure to form one component within a broader operational risk model.

None of the above are implemented. A user who wants any of these extensions would need to add them and the pipeline has been structured so that doing so does not require rewriting the existing stages.

---

## Citation

If you use this code in academic work, please cite the dissertation:

> Mills-Wierda, E. (2026). *Uncertainty-Embedded Flight Dynamics Modelling and Markov Chain Reachability Analysis for a UAV.* BEng dissertation, University of Glasgow.

---

## Licence

MIT. See `LICENSE`.
