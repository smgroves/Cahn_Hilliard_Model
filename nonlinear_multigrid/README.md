# Nonlinear Multigrid (NMG) Solvers

This directory contains NMG solver implementations in three languages. The **Julia implementation** is the primary production solver used for all manuscript results.

## Subdirectories

### `julia_multigrid/` — Primary Production Solver
The main solver used for all manuscript-quality results. Fast, feature-complete, and validated against the C reference. See [julia_multigrid/README.md](julia_multigrid/README.md).

### `python_multigrid/` — Reference (Slow)
Direct Python port of the C solver. Extremely slow; for demonstration and educational use only. See [python_multigrid/README.md](python_multigrid/README.md).

### `matlab_multigrid/`
MATLAB NMG solver. The canonical packaged MATLAB solvers (including NMG) are in `CHsolvers_package_copy_0206/`.

## Algorithm Overview

All three implementations use the **Full Approximation Scheme (FAS)** nonlinear multigrid method:

1. **Pre-smoothing**: Apply Gauss-Seidel relaxation on the fine grid
2. **Restriction**: Transfer residual and approximation to the coarse grid
3. **Coarse grid correction**: Recursively solve (or relax) on the coarser grid
4. **Prolongation**: Transfer correction back to the fine grid
5. **Post-smoothing**: Apply Gauss-Seidel relaxation again

The grid hierarchy goes from nx×nx down to 2×2, with `log2(nx)` levels.

## Cross-validation

The Julia and C solvers have been numerically compared at 256×256 resolution. Per-level phi and mu field comparisons are in `optimization_testing/julia_c_256_error/`.
