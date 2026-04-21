# 00_inputs: Platform parameters and blade-section polars

Shared inputs consumed by every downstream stage. Nothing is computed here and this folder just holds the data.

## Contents

```
00_inputs/
├── platform/
│   └── platform_parameters.m      Airframe parameters as a MATLAB function
└── polars/
    ├── afl13.c81.txt              Variant 1 blade-section polar
    ├── afl14.c81.txt              Nominal blade-section polar
    └── afl15.c81.txt              Variant 2 blade-section polar
```

## platform_parameters.m

Returns a struct `p` containing everything the BEMT and 6-DOF models need:

- Environmental constants (`p.g`, `p.rho`, `p.a_sound`)
- Masses (`p.m_base`, `p.m_pay`, `p.m`)
- Full 3×3 inertia tensor `p.I` and its inverse `p.Iinv` 
- Rotor hub positions `p.r_cg` and spin directions `p.spin`
- Fuselage aerodynamic reference (`p.Vref_fuse`, `p.Dref`, `p.Lref`, `p.Mref`)
- Rotor geometry (`p.R`, `p.Nb`, `p.r_blade`, `p.c_blade`, `p.theta_blade`)
- Pre-computed derived quantities: blade annular widths `p.dr_blade`, advance-ratio switch `p.mu_switch`, azimuth discretisation matrices

This is the **only** file in the repository that contains platform-specific numeric data. Every downstream script calls `p = platform_parameters()` and uses `p.field_name` accessors.

### How to adapt for a different UAV

Open the file and edit the values in place. The function name and all field names stay the same. Everything propagates automatically through the rest of the pipeline.

Items to edit:
- `p.m_base`, `p.m_pay`: dry mass and payload (kg)
- `Ixx` through `Iyz` (and the assembled `p.I`): inertia tensor components (kg·m²)
- `p.r_cg`: 3×4 matrix of rotor hub positions in body axes (m)
- `p.spin`: signs of each rotor's spin direction
- `p.Dref`, `p.Lref`, `p.Mref`, `p.Vref_fuse`: fuselage quadratic aerodynamic model
- `p.R`, `p.Nb`: rotor radius and blade count
- `p.r_blade`, `p.c_blade`, `p.theta_blade`: blade planform (radial stations, chord, twist)
- `p.mu_switch`: advance ratio threshold separating hover from forward flight (default 0.01)
- `M_azi`: azimuth discretisation count for forward-flight BEMT (default 12)

The derived fields (`p.dr_blade`, the `*_mat` matrices, `p.Iinv`) are computed from the primary values and shouldn't be edited directly.

## C81 polar files

Three standard C81-format tables containing lift and drag coefficients as functions of angle of attack and Mach number. The three samples are used to characterise blade-section manufacturing variation.

### How to adapt for different polars

Replace the three files with measurements from a different airfoil or a different manufacturing sample set. The repository convention is three samples, which matches the Tippett d₂ constant of 1.693 used throughout the Monte Carlo. If using a different number of samples, update `d2_n3` in every Monte Carlo script to the corresponding tabulated value.

The C81 row-column layout expected is:
- Row 1: Mach breakpoints (columns 2–11)
- Rows 2–100: CL table (column 1 = alpha, columns 2–11 = CL at each Mach)
- Row 101: Mach breakpoints repeated
- Rows 102–200: CD table (same layout)
