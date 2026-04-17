# MATLAB CH Solvers Package (February 2025)

This directory is a **curated package of all MATLAB Cahn-Hilliard solver implementations** as of February 6, 2025. It serves as the canonical MATLAB solver collection used by other parts of the repository (particularly `image_to_CH/`).

## Directory: `CahnHilliard_MATLAB_solvers/`

### Solver Functions

| File | Method | Description |
|------|--------|-------------|
| `CahnHilliard_FD.m` | Forward Euler (explicit) | Finite difference solver; small `dt` required for stability |
| `CahnHilliard_FD_SMG.m` | Forward Euler + SMG | FD with Smooth Multigrid variant |
| `CahnHilliard_NMG.m` | Nonlinear Multigrid (FAS) | Implicit solver; converges via V-cycles |
| `CahnHilliard_NMG_SMG.m` | NMG + SMG | NMG with Smooth Multigrid variant |
| `CahnHilliard_SAV.m` | Scalar Auxiliary Variable | Energy-stable implicit solver; returns phi in memory |
| `CahnHilliard_SAV_SMG.m` | SAV + SMG | SAV variant that streams phi to file (for large simulations) |

### Solver Interface (common parameters)

All solvers share a similar interface:

**Inputs:**
- `phi0` — Initial condition matrix (nx × ny)
- `t_iter` — Number of timesteps
- `dt` — Time step size
- `dt_out` — Save output every `dt_out` timesteps
- `m` — Interface mesh points (controls `epsilon`)
- `epsilon2` — Cahn number squared (can be provided directly)
- `boundary` — Boundary condition type (`'neumann'` or `'periodic'`)
- `c_relax` — Smoothing steps per V-cycle (NMG only)
- `solver_iter` — Max V-cycles per timestep (NMG only)
- `tol` — Convergence tolerance (NMG only)
- `domain` — Physical domain size
- `printres` — Print residuals (bool)
- `printphi` — Write phi to file (bool)
- `pathname` — Output directory path

**Outputs:**
- `t_out` — Vector of saved timestep indices
- `phi_t` — 3D array of phi (nx × ny × n_saved) or written to file
- `delta_mass_t` — Mass change at each saved timestep
- `E_t` — Discrete energy at each saved timestep

**SAV-specific parameters:** `C0`, `Beta`, `gamma0` (modify the auxiliary energy)

### Utility Functions

| File | Purpose |
|------|---------|
| `ch_initialization.m` | IC generator: `spinodal`, `unirand`, `droplet` (tanh), `geometric` (CPC rectangle+circle), `file` (from CSV) |
| `ch_movie.m` | Create MP4 animation from phi_t 3D array |
| `ch_movie_from_file.m` | Create animation from saved phi CSV (for large sims) |
| `ch_discrete_energy.m` | Compute discrete Cahn-Hilliard free energy |
| `ch_laplace.m` | Discrete Laplacian (Neumann BCs) |
| `Lap_SAV.m` | Laplacian with periodic BCs via FFT (for SAV) |
| `ch_smooth.m` | Gauss-Seidel smoothing operator |
| `nmg_solver.m` / `nmg_vcycle.m` | NMG solver components |
| `sav_solver.m` | SAV solver core |
| `fd_solver.m` | FD solver core |
| `f.m` / `df.m` | Free energy `f(φ)` and its derivative `f'(φ)` |
| `b_fun.m`, `r_fun.m`, `r0_fun.m` | SAV auxiliary variable functions |
| `g_fun_CN.m`, `A_inv_CN.m` | SAV Crank-Nicolson helper functions |
| `fft2_filtered.m` | FFT-based filtered Laplacian for SAV periodic BCs |
| `ch_error2.m` | Residual (Frobenius norm) computation |
| `read_phi_from_file.m` | Read phi output files from disk |
| `run_spinodal_decomp.m` | Example run script |

## Usage

This package is called from other parts of the repository:
- `image_to_CH/run_CH_on_image.m` calls `CahnHilliard_SAV_SMG` from this package
- Add this directory to your MATLAB path before calling any function

## Job Logging

`Job_specs.txt` records MATLAB simulation runs with parameters and runtimes. The top-level `Job_specs.csv` records runs across all languages.
