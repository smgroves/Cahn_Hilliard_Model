# Original C Nonlinear Multigrid Solver

This directory contains the **reference C implementation** of the nonlinear multigrid (NMG) solver for the Cahn-Hilliard equation, based on Lee et al. (2014). All other solver implementations in this repository (Julia, Python, MATLAB) are derived from or validated against this code.

## Algorithm

The solver uses the **Full Approximation Scheme (FAS)** multigrid method:
- **Relaxation**: Gauss-Seidel smoothing, solving a 2×2 system at each grid point via Cramer's rule
- **Restriction**: Full-weighting operator (coarsen by 2×)
- **Prolongation**: Constant interpolation (expand by 2×)
- **Boundary conditions**: Neumann (zero-flux)
- **Time stepping**: Implicit, converged via V-cycles until residual < tolerance

## Source Files

| File | Description |
|------|-------------|
| `CHsolver.c` | Base solver, 64×64 grid |
| `CHsolver_CPC.c` | CPC geometry variant |
| `CHsolver_CPC_big.c` | CPC variant for 256×256 grid |
| `CHsolver_CPC_big_from_file.c` | CPC variant reading IC from CSV file |
| `CHsolver_CPC_circle_IC.c` | CPC variant with circular initial condition |
| `CHsolver_CPC_compare_to_py.c` | Comparison variant against Python output |
| `CHsolver_convergence.c` | Convergence testing variant |
| `CHsolver_CPC_big128.c` | CPC variant for 128×128 grid |

## Core Functions (in each solver)

| Function | Purpose |
|----------|---------|
| `laplace()` | Discrete Laplacian with Neumann BCs |
| `relax()` | Gauss-Seidel smoother (2×2 system per grid point) |
| `restrictCH()` | Full-weighting restriction (fine → coarse grid) |
| `prolong_ch()` | Constant prolongation (coarse → fine grid) |
| `vcycle()` | Recursive FAS V-cycle |
| `defect()` | Defect computation for FAS |
| `error2()` | Frobenius norm of residual |
| `cahn()` | Time-stepping loop (calls vcycle until convergence) |

## Pre-compiled Executables

| Executable | Grid size |
|-----------|-----------|
| `solver_16` | 16×16 |
| `solver_32` | 32×32 |
| `solver_64` | 64×64 |
| `solver_128` | 128×128 |
| `solver_256` | 256×256 |
| `CHsolver` | 64×64 (base) |
| `CHsolver_CPC` | CPC geometry |

## Scripts

- `runtime_solver.sh` — Batch benchmark script, runs compiled executables with varying `max_it` and `max_it_CH` parameters
- `src/` — Alternative organized source directory

## Output

Solvers write the phi field as text files to `outputs/` or `out/`, with one NxN block per saved timestep appended consecutively.

## Relationship to Other Solvers

This C code is the numerical ground truth for:
- **Julia NMG** (`nonlinear_multigrid/julia_multigrid/`) — primary production solver, validated against C
- **Python NMG** (`nonlinear_multigrid/python_multigrid/`) — direct Python port for demonstration
- **MATLAB NMG** (`CHsolvers_package_copy_0206/CahnHilliard_NMG.m`) — MATLAB port

Detailed numerical comparisons between Julia and C at 256×256 are in `optimization_testing/julia_c_256_error/`.

## Reference

Lee, H.G. et al. (2014). A second-order accurate non-linear difference scheme for the N-component Cahn–Hilliard system. *Physica A*.

