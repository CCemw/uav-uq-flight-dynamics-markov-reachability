# 01_uncertainty_quantification: Polar uncertainty characterisation

Visualises the blade-section polar measurements and derives the pointwise standard deviation envelope that feeds every downstream Monte Carlo.

## Contents

```
01_uncertainty_quantification/
├── plot_polar_uncertainty.m     Produces the polar envelope figures
└── README.md
```
## What the script does

`plot_polar_uncertainty.m` loads the three C81 tables from `../00_inputs/polars/`, extracts the lift and drag coefficients at Mach 0.30, and produces two figures:

- `figs/Fig05_Cl_Cd_polars.pdf`: the three measured polars plotted together
- `figs/Fig05b_Cl_Cd_sigma.pdf`: the mean polar with a pointwise standard deviation band

The σ grid is computed as `(max - min) / d2_n3` with `d2_n3 = 1.693` for n = 3 samples. This same grid is rebuilt inside every downstream Monte Carlo script, so editing the polars or the d₂ constant in `00_inputs/` is sufficient to change the uncertainty model everywhere.

## How to adapt

- **Mach number.** Edit `M_plot` at the top of the script to visualise the polars at a different slice through the (α, M) grid.
- **Number of samples.** If the number of polar files in `00_inputs/polars/` changes, update `d2_n3` to the appropriate Tippett constant.
- **Alpha range.** Edit the two `xlim([-25 25])` calls for a different axis range.

## How to run

```
matlabcd 01_uncertainty_quantification
plot_polar_uncertainty
```
Produces two PDFs in `figs/`.
