# Optimization Testing

This directory contains benchmarking, performance testing, and numerical validation work for the CH solvers — primarily comparing the Julia and C implementations.

## Subdirectories

### `julia_early_tesing/`

Early Julia solver development and iteration. Contains:

| File/Directory | Purpose |
|----------------|---------|
| `solver_type_stable_v4.jl` | Julia solver iteration (type-stable version 4) |
| `solver_type_stable_v5.jl` | Julia solver iteration (type-stable version 5) |
| `run_solver.jl` | Run script for early versions |
| `Project.toml` / `Manifest.toml` | Julia package environment (dependencies) |
| `test_256/` | Intermediate numerical outputs for debugging at 256×256 |

These files document the development process before the solver was moved to `nonlinear_multigrid/julia_multigrid/`.

### `julia_c_256_error/`

Detailed **numerical comparison between the Julia and C solvers** at 256×256 resolution. This is the primary validation evidence that the Julia solver correctly reimplements the C algorithm.

Contents:
- PNG images of phi and mu fields at each V-cycle level, per multigrid level
- `plot_error_mu_julia_c.py` — Python script generating the comparison plots
- `logging_plot_error_mu_julia_c.ipynb` — Jupyter notebook version with logged outputs
- `relax/` — Per-level comparison after relaxation step
- `tan_IC/` — Comparisons with tanh initial condition
- `vcycle/` — Per-level comparison through full V-cycle
- `MinJhe_Mac/` — Comparison runs on Min-Jhe Lu's machine (cross-platform validation)
- `output_C_on_mac_MinJhe/` — Raw text outputs from C solver for comparison

The comparison checks phi and mu at each grid level of the multigrid hierarchy, verifying that Julia and C produce identical results (within floating-point tolerance).
