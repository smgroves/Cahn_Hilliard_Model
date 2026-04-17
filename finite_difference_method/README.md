# Finite Difference Method (Forward Euler)

This directory contains a standalone MATLAB **explicit Forward Euler finite difference** solver for the Cahn-Hilliard equation. This is the simplest solver in the repository — straightforward but requires very small time steps for stability.

## Algorithm

**Method:** Forward Euler (explicit time integration) with second-order finite differences in space.

```
φⁿ⁺¹ = φⁿ + dt · Δ((φⁿ)³ - φⁿ - ε² Δφⁿ)
```

Since this is explicit, stability requires `dt ≤ h⁴ / (4ε²)` approximately. In practice, much smaller `dt` than the NMG solver is needed.

## Files

### `spinodal_decomp.m` — Main Solver Function

**Signature:** `spinodal_decomp(D, gamma, options)`

**Inputs:**
- `D` — Diffusion coefficient (scales the time derivative)
- `gamma` — Ratio `ε²/h²` (encodes Cahn number relative to grid spacing)
- `options` struct:
  - `dt` — Time step
  - `GridSize` — nx (grid is nx × nx)
  - `NumIterations` — Total timesteps
  - `FrameSpacing` — Save every N timesteps
  - `CaptureMode` — `'all'` or `'final'`
  - `ImgStyle` — Colormap style
  - `FileName` — Output filename base
  - `InputMatrix` — Path to CSV file for initial condition
  - `InputType` — IC type (`'csv'`, `'random'`, etc.)
  - `ConstantColorbar` — Fix colorbar range (bool)
  - `write_phi` — Write phi to CSV (bool)
  - `write_residual` — Write residual to CSV (bool)

**Outputs:**
- `phi_t` — 3D array (nx × ny × n_frames)
- `mass_t` — Mass at each saved frame
- `E_t` — Discrete energy at each saved frame
- MP4 video of phi evolution

### `laplacian.m`
Discrete Laplacian operator with **Neumann (zero-flux) boundary conditions** using finite differences. This is the corrected version; see `output/obsolete_output/obsolete_laplace/` for the old (incorrect) implementation.

### `ch_movie.m`
Generates MP4 animation from a 3D phi array.

## Run Scripts

All run scripts read initial conditions from the Julia NMG solver outputs and write results to `output/`.

| Script | IC Source | Scenario |
|--------|-----------|----------|
| `run_fd_method.m` | Various | General FD runs |
| `run_fd_spinodal_MG_timepoint.m` | `nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_MG_timepoint_IC/` | Spinodal from NMG timepoint |
| `run_FD_sd_smooth.m` | `nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/` | Spinodal from smooth IC |
| `run_fd_droplets.m` | `nonlinear_multigrid/julia_multigrid/manuscript_output/large_and_small_droplets/` | Two-droplet IC |
| `run_fd_checkerboard.m` | Checkerboard pattern | Checkerboard IC |
| `run_fd_flat_interface.m` | `nonlinear_multigrid/julia_multigrid/manuscript_output/flat_interface/` | Flat interface |

## Output

Results are written to `output/` with standardized filenames encoding simulation parameters:
```
{solver}_{iterations}_dt{dt}_N{nx}_gamma{gamma}_D{D}_{IC_type}
```

Subdirectories:
- `checkerboard_IC/`
- `flat_interface/`
- `large_and_small_droplets/`
- `spinodal_+1_-1/`
- `spinodal_MG_timepoint_IC/`
- `spinodal_smooth_relax_function/`

Each subdirectory may contain:
- `obsolete_laplace/` — Results from old (incorrect) Laplacian
- `periodic_BC_laplace/` — Results with periodic BC Laplacian for comparison

## Connection to Other Solvers

FD solver results are compared against NMG results in `plotting/comparison_sim_results_across_methods/` and `plotting/plot_energy_mass_across_methods.py`. The FD solver uses the **same initial conditions** as the NMG and SAV solvers for fair comparison.
