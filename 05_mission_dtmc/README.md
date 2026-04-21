# 05_mission_dtmc: Mission-level Markov chain risk model

Evolves a four-state discrete-time Markov chain (Normal Flight, Degraded Flight, Abort, Crash) through a phased mission. Per-phase transition probabilities are derived by convolving the trim probability envelope from Stage 2 with a low-altitude Dryden gust distribution at the phase altitude. The final state distribution gives the mission outcome probabilities, and the crash probability divided by mission duration gives the per-flight-hour crash rate for comparison against JARUS SORA SAIL thresholds.

## Contents

```
05_mission_dtmc/
├── mission_dtmc.m              Full DTMC mission model, all figures and tables
└── README.md
```

There are no `.mat` files in this folder. The DTMC derives **EVERYTHING** at runtime from the trim Monte Carlo output loaded from `../02_trim_monte_carlo/results/`.

## Workflow

1. Ensure Stage 2 has been run and `mc_trim_results_QAV250.mat` exists.
2. Run `mission_dtmc.m` from this folder. The whole run takes a few seconds.

Produces the four mission-evolution figures and two console tables (terminal outcomes at the two reference cruise speeds and V_max compliance at each SAIL level).

## How to adapt

### Mission profile

The mission is a cell array at the top of the script:

```matlab
mission = {
    'Takeoff', 3,  10, 30;
    'Cruise',  12, 60, 600;
    'Descent', 12, 40, 120;
    'Loiter',  5,  30, 180;
    'Landing', 3,  10, 30
};
```

Each row is one phase: name, nominal airspeed (m/s), altitude (m), duration (s). Add rows, remove rows, or edit values freely. The "Cruise" phase (second row) receives a special substitution: its airspeed is overridden by the cruise speed being analysed, so the user can vary cruise speeds without editing the profile.

### Contingency routing parameters

Three scalars at the top of the script:

- `eta`: probability that a sample in the Abort state successfully lands
- `eta_DF`: probability that a sample in the Degraded Flight state lands safely when routed to the abort sub-phase
- `s`: severity-routing parameter, the fraction of trim-loss events that go straight to the terminal branch rather than to Degraded Flight

These are the dominant free parameters in the DTMC. The current values are chosen for a plausible urban mission with moderate assurance. In my actual dissertation, Appendix F.4 presents the sensitivity across a factorial grid and shows that the qualitative result (manufacturing uncertainty produces a meaningful crash rate well below the deterministic envelope) is robust to the choice.

### Atmospheric model

`W20_kt` (wind at 20 ft, in knots), `h_urban_table` and `f_urban_table` (altitude-to-urban-intensity mapping), and `f_Lu` (urban correlation-length factor) are exposed at the top. The low-altitude Dryden form is MIL-F-8785C; the urban correction factors come from the literature survey in dissertation Chapter 2.

### Cruise speeds analysed

- `V_pair` (default `[12, 17]`): the two cruise speeds shown side-by-side in the state-evolution figures.
- `V_sweep` (default `8:1:22`): the grid used for the crash-rate vs. cruise-speed figure.

### SAIL thresholds

The current JARUS SORA thresholds can be replaced with any other regulatory set (for example FAA ACs or EASA certification targets) by editing `sail_targets` and `sail_names`.

### Abort sub-phase

The abort sub-phase is defined as a single struct at the top. The abort is entered from the Degraded Flight state at the end of any non-final phase and terminates in either a successful emergency landing (with probability η) or a crash.

## Extensions the structure supports

- **Data-driven transition probabilities.** The current per-phase probabilities come from convolving the trim envelope with a Dryden gust distribution. A user with measured or simulated encounter data (flight-test logs, high-fidelity gust simulations, operational incident records) could replace `compute_phase_probabilities` with a function that returns `p_N`, `p_DF`, `p_L` from the data directly, leaving the rest of the chain evolution untouched.

- **Additional risk factors.** The transition into the trim-loss branch is currently driven by aerodynamic uncertainty only. Other contributors, such as propulsion faults, sensor failures, navigation errors, or traffic encounters, could be included by adding their per-phase probability into the existing `p_L` term, or by introducing new terminal branches. The DTMC framework does not depend on the physical origin of the probabilities, only that they sum correctly across each row.

- **Further contingency states.** Additional states such as "battery critical", "motor fault", "sensor drop-out", or "communications lost" can be added by extending the 4×4 transition matrix and defining the corresponding rows. The trim-derived probability construction for the existing states carries over unchanged.

- **Richer contingency routing.** The between-phase routing currently sends the full Degraded Flight residual to a single abort sub-phase. This block can be replaced with any user-defined routing logic (mission-specific emergency landing sites, conditional branching based on altitude or remaining mission time, or multiple competing contingency options, each with its own sub-phase).

- **Multiple atmospheric scenarios.** The current script evaluates at a single `W20`. A sweep over wind strength or urban intensity produces an additional axis of risk variation and is a natural starting point for envelope analysis.

## Numerical method

The chain evolution uses raw matrix multiplication rather than MATLAB's built-in `dtmc` class. This keeps the repository free of toolbox dependencies (the `dtmc` class requires the Econometrics Toolbox), and is the natural form for this problem because each mission phase has its own phase-specific transition matrix rather than a single chain.
