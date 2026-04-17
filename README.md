# Cahn-Hilliard Model

This repository implements and studies numerical solvers for the **Cahn-Hilliard (CH) equation**, a fourth-order PDE that models phase separation in binary mixtures (spinodal decomposition, droplet dynamics). The primary scientific application is modeling **Chromosomal Passenger Complex (CPC) condensate formation** on chromosomes during cell division in HeLa and MCF10A cell lines.

## The Cahn-Hilliard Equation

The CH equation solved throughout this codebase is:

```
∂φ/∂t = Δμ
μ = φ³ - φ - ε² Δφ
```

where:
- `φ` is the order parameter (phase field), ranging from -1 to +1
- `ε` (Cahn number) controls the interface width
- `μ` is the chemical potential
- The double-well free energy is `f(φ) = (1/4)(φ² - 1)²`

An asymmetry parameter `α` can be added to the free energy for CPC condensate physics.

## Repository Structure

```
Cahn_Hilliard_Model/
├── original_C_NMG/          # Reference C implementation (Lee et al. 2014)
├── nonlinear_multigrid/     # NMG solvers in Julia, Python, and MATLAB
│   ├── julia_multigrid/     # Primary production Julia solver
│   ├── python_multigrid/    # Python port (slow, for reference only)
│   └── matlab_multigrid/    # MATLAB NMG solver
├── CHsolvers_package_copy_0206/  # Packaged MATLAB solvers (canonical Feb 2025)
├── finite_difference_method/     # Forward Euler explicit FD solver (MATLAB)
├── SAV/                     # Scalar Auxiliary Variable solver (MATLAB/Julia)
├── finite_element_method/   # FEniCS FEM prototype (Python, slow)
├── image_to_CH/             # Microscopy image → CH initial condition pipeline
├── vcell_to_chsolver/       # VCell simulation → CH initial condition bridge
├── plotting/                # Post-processing, analysis, and visualization
├── optimization_testing/    # Benchmarking and solver comparison
├── utils/                   # General utilities (file merging, etc.)
└── Job_specs.csv            # Centralized runtime log for all simulations
```

## Data Flow

The overall pipeline connects these components:

```
VCell exports ──────────────────────────────────────────────┐
Microscopy images (JPEGs) ──────────────────────────────────┤
                                                            ▼
                                              IC generation (CSV files)
                                                            │
                              ┌─────────────────────────────┤
                              ▼                             ▼
              Julia NMG solver               MATLAB solvers (NMG/FD/SAV)
     (nonlinear_multigrid/julia_multigrid)   (CHsolvers_package_copy_0206)
                              │                             │
                              └──────────┬──────────────────┘
                                         ▼
                                 phi.txt / phi.csv outputs
                                         │
                                         ▼
                              plotting/ (kymographs, energy,
                              droplet analysis, figures)
```

**Key connections:**
- `finite_difference_method/` run scripts read ICs from `nonlinear_multigrid/julia_multigrid/manuscript_output/`
- `SAV/` run scripts read ICs from `nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/IC/`
- `image_to_CH/` calls solvers from `CHsolvers_package_copy_0206/`
- `vcell_to_chsolver/` writes ICs to `data/` for use in other solvers
- `plotting/` reads outputs from all solvers

## Rivanna HPC

Several large-scale simulations were run on **Rivanna**, the University of Virginia HPC cluster, submitted as SLURM job arrays. Results are stored locally in subdirectories named `from_rivanna/` or `from_Rivanna/`. The companion repository **Cahn_Hilliard_Model_HPC** contains the job submission scripts and Julia solver code used on Rivanna.

Key Rivanna-generated data:
- `nonlinear_multigrid/julia_multigrid/manuscript_output/critical_radius/from_Rivanna/` — critical radius vs. epsilon parameter sweeps
- `plotting/` references `domain_0_2_from_rivanna_kymographs_e_0.0075` — kymograph data from job arrays
- `utils/combine_all_folders.py` — merges two batches of Rivanna kymograph output

## Solvers Overview

| Solver | Language | Method | Speed | Use Case |
|--------|----------|--------|-------|----------|
| Julia NMG | Julia | Nonlinear Multigrid (FAS) | Fast | Primary production solver |
| C NMG | C | Nonlinear Multigrid (FAS) | Fastest | Reference implementation |
| MATLAB NMG | MATLAB | Nonlinear Multigrid (FAS) | Moderate | Cross-validation |
| MATLAB FD | MATLAB | Forward Euler explicit | Slow (small dt) | Simple reference |
| MATLAB SAV | MATLAB | Scalar Auxiliary Variable | Fast, energy-stable | Alternative implicit solver |
| Python NMG | Python | Nonlinear Multigrid (FAS) | Very slow | Demonstration only |
| Python FEM | Python | Finite Element (FEniCS) | Very slow | Prototype |

## Key Parameters

| Parameter | Meaning | Typical Value |
|-----------|---------|---------------|
| `nx`/`ny` | Grid points per side | 128, 256 |
| `h` | Spatial step size = 1/nx | 1/128 ≈ 0.0078 |
| `dt` | Time step | 1e-5 to 2.5e-5 |
| `epsilon` | Interface width (Cahn number) | derived: `m*h/(2√2·atanh(0.9))` |
| `m` | Interface mesh points | 8 |
| `tol` | NMG convergence tolerance | 1e-6 |
| `alpha` | Free energy asymmetry for CPC | -0.5, -0.2, 0.0, 0.2 |

## Contributors

- **C code**: Lee et al. (2014) — original NMG reference implementation
- **Julia code**: Sarah Groves (2024)
- **Python code**: Sarah Groves (2024)
- **MATLAB SAV code**: Min-Jhe Lu (2024)
- **MATLAB NMG/FD code**: Sarah Groves (2024)

## Job Logging

All simulation runs are logged in `Job_specs.csv` (top level) with columns: date, simulation name, language, grid size, epsilon, dt, tolerance, timesteps, max multigrid iterations, and wall-clock time. MATLAB runs also log to `CHsolvers_package_copy_0206/Job_specs.txt`.
