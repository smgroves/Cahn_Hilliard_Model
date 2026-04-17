# Scalar Auxiliary Variable (SAV) Solver

This directory contains the **Scalar Auxiliary Variable (SAV)** method solver for the Cahn-Hilliard equation, implemented in MATLAB (with an early Julia prototype). SAV is an alternative to NMG that guarantees **unconditional energy stability** and allows larger time steps than explicit methods.

## Method

SAV reformulates the CH equation by introducing a scalar auxiliary variable:

```
r(t) = sqrt(E₁(φ) + C₀)
```

where `E₁(φ)` is the nonlinear part of the free energy and `C₀ > 0` is a constant ensuring the square root is always positive. This transforms the nonlinear problem into a sequence of linear solves. The energy stability is guaranteed regardless of time step size.

**Parameters specific to SAV:**
- `C0` — Constant ensuring `E₁ + C0 > 0`
- `Beta` — Regularization parameter
- `gamma0` — Auxiliary energy parameter

## Development History

The directory structure traces the SAV solver development timeline:

```
SAV/
├── old_code/                          # Feb 2024 — initial development
│   ├── SAV_filter_clean_up_02282024.m # Early MATLAB SAV implementation
│   ├── SAV_filter.jl                  # Julia SAV prototype (FFT-based)
│   ├── CH1d_SAV_1st_order.mlx         # 1D SAV verification
│   └── compare_SAV_multigrid/         # CPU timing: SAV vs. NMG
├── old_code_doesn't_work_11-2024/     # Nov 2024 — abandoned intermediate version
├── CH_SAV_20241213/                   # Dec 2024 — improved version
│   └── run_SAV_droplets.m             # Two-droplet IC test
├── CH_SAV_20250106/                   # Jan 2025 — current version
│   ├── run_SAV_droplets.m             # Two-droplet IC simulation
│   └── run_SAV_spinodal.m             # Spinodal decomposition simulation
└── output/                            # Simulation results
```

**The `CH_SAV_20250106/` directory contains the current working code.**

## Current Version (`CH_SAV_20250106/`)

### `run_SAV_spinodal.m`
Runs SAV solver on a spinodal decomposition initial condition.

**IC source:** Reads from `nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/IC/`  
**Output:** Writes to `SAV/output/spinodal_smooth_relax_function/`

### `run_SAV_droplets.m`
Runs SAV solver on a two-droplet initial condition with parameters `r1`, `r2`, `spacing`.

**Output:** Writes to `SAV/output/large_and_small_droplets/`

## Output

`output/` contains:
- `large_and_small_droplets/` — Two-droplet simulations
- `spinodal_smooth_relax_function/` — Spinodal decomposition simulations

Each contains MP4 movies and PNG snapshots.

## Julia Prototype (`old_code/SAV_filter.jl`)
An early Julia implementation of SAV using FFT-based solvers. Not maintained; the MATLAB version is the production implementation.

## Comparison to NMG

SAV results are compared to NMG and MATLAB results in:
- `plotting/plot_MATLAB_vs_Julia_SAV.py` — MATLAB vs. Julia SAV comparison
- `plotting/comparison_sim_results_across_methods/` — FD, NMG, SAV comparison figures

The canonical MATLAB SAV solver is in `CHsolvers_package_copy_0206/CahnHilliard_SAV.m` and `CahnHilliard_SAV_SMG.m`.
