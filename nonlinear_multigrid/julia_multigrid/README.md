# Julia NMG Solver

The primary production solver for the Cahn-Hilliard equation. Implements the Full Approximation Scheme (FAS) nonlinear multigrid method in Julia. All manuscript figures are generated using this solver (or its Rivanna HPC version).

## Solver Files

Multiple solver variants exist for different use cases:

| File | Description |
|------|-------------|
| `CH_multigrid_solver.jl` | Base solver, Neumann (zero-flux) BCs |
| `CH_multigrid_solver_periodic_bc.jl` | Periodic boundary condition variant |
| `CH_multigrid_solver_with_alpha.jl` | Adds asymmetry parameter `α` to free energy for CPC physics |
| `CH_multigrid_solver_with_alpha_v2.jl` | Updated `α` version |
| `CH_multigrid_solver_with_alpha_change_domain.jl` | `α` version with variable domain size |
| `CH_multigrid_solver_large_domain.jl` | Large domain variant |
| `CH_multigrid_solver_optimized.jl` | Performance-optimized version |

### The `α` Parameter

The `_with_alpha` variants add an asymmetry term to the double-well free energy:

```
f(φ) = φ⁴/4 - φ²/2 + α·φ
```

This breaks the symmetry between the two phases and is used to model CPC condensate formation, where one phase is energetically preferred.

## Main Entry Point: `multigrid_solver()`

All solver variants expose a `multigrid_solver()` function with the following interface:

**Required arguments:**
- `oc` — Initial phi matrix (nx × nx)
- `nx` — Grid size
- `tol` — Convergence tolerance (e.g., `1e-6`)
- `outdir` — Output directory path

**Named parameters:**
- `max_it` — Max V-cycles per timestep
- `max_it_CH` — Max CH timesteps
- `dt` — Time step size
- `m` — Interface mesh points (controls epsilon; typically `8`)
- `epsilon` — Interface width (if provided directly, overrides `m`)
- `ns` — Save output every `ns` timesteps
- `print_phi` — Write phi field to file (bool)
- `print_mass` — Write mass to file (bool)
- `print_e` — Write energy to file (bool)
- `print_r` — Write residual to file (bool)
- `overwrite` — Overwrite existing output (bool)
- `check_dir` — Check if output directory exists (bool)

## Output Files

All outputs are written to `outdir` as plain text files. For each saved timestep, an nx×nx phi block is appended:

| File | Contents |
|------|----------|
| `*_phi.txt` | Phase field φ at each saved timestep (NxN blocks, concatenated) |
| `*_mass.txt` | Total mass at each saved timestep |
| `*_energy.txt` | Discrete free energy at each saved timestep |
| `*_residual.txt` | NMG residual at each saved timestep |

Use `get_last_n_lines.py` to extract the last N timestep snapshots from a `phi.txt` file.

## Initialization Functions

| Function | Description |
|----------|-------------|
| `initialization_from_function()` | Analytical IC (e.g., tanh droplet profile) |
| `initialization_random()` | Uniform random noise |
| `initialization_spinodal()` | Small random perturbation around mean for spinodal decomposition |
| `initialization_from_file()` | Read phi from a CSV file |
| `initialize_geometric_CPC()` | Rectangle + circle geometry for CPC/chromosome modeling |
| `initialize_round_CPC()` | Circular CPC geometry |
| `initialization()` | Dispatcher — selects IC based on keyword argument |

## Diagnostics

- `calculate_mass()` — Total mass (integral of phi)
- `calculate_discrete_energy()` — Discrete Cahn-Hilliard free energy
- `calculate_discrete_norm_energy()` — Normalized energy

## Run Scripts

### `run_solver_example.jl`
Minimal usage example showing two scenarios: spinodal decomposition and initialization from file. Good starting point for new users.

### `run_solver_manuscript.jl`
Full parameter sweeps for all manuscript figures. Includes:
- Droplet simulations (critical radius, large/small droplets)
- Spinodal decomposition (multiple IC types)
- CPC geometry simulations (varying `CPC_width`, `cohesin_width`, `alpha`)
- Some parameter sweeps were run as SLURM job arrays on the **Rivanna HPC cluster** (see comments in file)

## Utility

### `get_last_n_lines.py`
Reads the last N lines from a `phi.txt` output file, corresponding to the last N saved timestep snapshots. Useful when only the final state is needed.

## Output Data (`manuscript_output/`)

Subdirectories contain simulation results for each scenario:

| Directory | Contents |
|-----------|----------|
| `CPC_geometry_m=8/` | CPC/cohesin geometry sims at varying widths |
| `critical_radius/` | Droplet critical radius vs epsilon; `from_Rivanna/` has Rivanna results |
| `flat_interface/` | Flat interface benchmark |
| `large_and_small_droplets/` | Two-droplet simulations |
| `spinodal_+1_-1_IC/` | Spinodal IC from ±1 initial condition |
| `spinodal_MG_timepoint_IC/` | Spinodal IC from NMG timepoint |
| `spinodal_smooth_relax_function/` | Spinodal IC from smoothed relaxation; `IC/` subdirectory has ICs used by FD and SAV solvers |
| `spinodal_normal_IC/` | Normally distributed IC |
| `stripe_IC/` | Stripe pattern IC |
| `optimized_testing/` | Performance benchmarking outputs |

### Connection to Other Solvers

The `spinodal_smooth_relax_function/IC/` subdirectory serves as a shared IC source:
- `finite_difference_method/` run scripts read ICs from here
- `SAV/CH_SAV_20250106/run_SAV_spinodal.m` reads ICs from here

### `critical_radius/initialize_alt_visual.jl`
Script used to test alternative ICs for critical radius simulations on Rivanna.

### `critical_radius/critical_radius_vs_epsilon.py`
Python plotting script that reads from `from_Rivanna/` and `from_Rivanna_old/` subdirectories to generate critical radius vs. epsilon figures.

## Rivanna HPC

Several large-scale parameter sweeps were submitted as SLURM job arrays on Rivanna. The companion repository **Cahn_Hilliard_Model_HPC** contains the job submission scripts. Results are stored locally in `from_Rivanna/` subdirectories.
