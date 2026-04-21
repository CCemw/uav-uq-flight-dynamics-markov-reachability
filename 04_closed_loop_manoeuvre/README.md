# 04_closed_loop_manoeuvre: PID-stabilised manoeuvre response

Time-domain simulation of the closed-loop pitch response to a symmetric rear-rotor RPM pulse, under a cascaded angle/rate PID controller, nominal and Monte Carlo.

## Contents

```
04_closed_loop_manoeuvre/
├── closed_loop_manoeuvre_mc.m           Trim + closed-loop simulation
├── plot_manoeuvre_response.m            Plots the dissertation manoeuvre figures
└── results/
    ├── manoeuvre_V5_nominal_QAV250.mat
    ├── manoeuvre_V5_mc_QAV250.mat
    ├── manoeuvre_V15_nominal_QAV250.mat
    └── manoeuvre_V15_mc_QAV250.mat
```

## Workflow

The core script runs **one case at a time** based on two settings at the top of the file: `V_trim` (5 or 15) and `RUN_MODE` (`'nominal'` or `'mc'`). Editing these and re-running produces one of the four `.mat` files listed above. The filename is auto-generated from these two settings so there's no risk of overwriting.

1. **Run `closed_loop_manoeuvre_mc.m`** four times with the four combinations. Each run takes seconds (nominal) to about 25 minutes (MC at V=15).

2. **Run `plot_manoeuvre_response.m`** for each airspeed. The plotting script's `V_plot` setting selects which pair of `.mat` files to load. Produces `figs/Fig10_manoeuvre_V5.pdf` and `figs/Fig10_manoeuvre_V15.pdf`.

## How to adapt

### Trim airspeed

`V_trim` at top of core script. The script trims internally up to V_trim + 2 m/s to have guess material for the solver, then simulates at V_trim.

### Disturbance pulse

`dRPM_pulse` (magnitude) and `t_pulse_dur` (duration) at top of core script. The pulse is applied symmetrically to both rear rotors, producing a pure nose-up pitch disturbance.

### Controller gains

All four PID gains (`Kp_angle`, `Kp_rate`, `Ki_rate`, `Kd_rate`) and the derivative low-pass cutoff `f_D_cutoff` are explicitly at the top of the script. The current values are tuned for the QAV250; a different UAV will require re-tuning.

### Divergence threshold

`theta_diverge` (default 60°). User-set abort criterion: if the pitch angle exceeds this magnitude during the simulation, the run is logged as diverged and the time history is held at the last valid value for the remainder of the window. (Arbitrary "divergence", can be adjusted or removed completely)

### Mission-style extensions

The current simulation is a single-axis pulse response. The RK4 integration loop and controller block are written so a full trajectory controller (position loops, waypoint following) could work with the same 8-state EOM and trim engine without changing the uncertainty propagation structure.

### Gain scheduling

The controller uses fixed gains regardless of airspeed. Adding a simple gain schedule (for example, looking up `Kp_rate` from a table indexed on `V_trim`) would go inside the simulation loop without changing the trim or MC structure.
