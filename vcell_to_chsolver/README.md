# VCell to CH Solver

This directory bridges **mechanistic signaling models** (from VCell, the Virtual Cell computational biology platform) to the Cahn-Hilliard phase-field simulations in this repository.

## Purpose

VCell simulations model the detailed biochemical reaction-diffusion dynamics of CPC (Chromosomal Passenger Complex) proteins during mitosis. The output of VCell (spatial protein concentration maps) can serve as physically-motivated initial conditions for the CH model, replacing synthetic ICs like tanh droplets or random noise.

## Pipeline

```
VCell simulation exports
(~/Box/CPC_Model_Project/VCell_Exports/)
          │
          ▼
   normalize_CPC.py
   - Sum all CPC species at a timepoint
   - Rescale by factor (default: 10)
   - Pad to square grid
   - Optionally rescale to φ ∈ [-1, 1]
   - Generate 1×, 2×, 4× resolutions
          │
          ▼
   data/*.csv    (IC files at multiple resolutions)
          │
          ▼
   CH solvers (julia_multigrid, CHsolvers_package_copy_0206, etc.)
```

## File

### `normalize_CPC.py`

**Inputs:**
- VCell export CSVs from `~/Box/CPC_Model_Project/VCell_Exports/` (not tracked in this repo)
- `timepoint` — Which simulation timepoint to use (100–400 s available)
- `scale_factor` — Rescaling factor for concentrations (default: 10)

**Processing:**
1. Sum all CPC species: `CPCi`, `CPCa`, `pH2A_Sgo1_CPCa`, and other variants
2. Rescale by `scale_factor`
3. Pad to square (chromosome shape is rectangular)
4. Optionally rescale to φ ∈ [-1, 1]
5. Apply `prolong()` (nearest-neighbor upsampling) to generate higher-resolution ICs

**Outputs (to `data/`):**
- `*_1x.csv` — Original resolution
- `*_2x.csv` — 2× upsampled
- `*_4x.csv` — 4× upsampled

Multiple VCell simulation IDs are referenced in commented-out code, corresponding to different CPC model conditions (relaxed/tensed chromosome configurations) at timepoints from 100–400 s.

## VCell Models Referenced

The commented code references several VCell simulation IDs corresponding to different biological conditions:
- **Relaxed CPC model** — chromosome in relaxed conformation
- **Tensed CPC model** — chromosome under tension (as during chromosome biorientation)

These simulations are stored in `~/Box/CPC_Model_Project/VCell_Exports/` (external to this repository).

## Connection to Rest of Repository

The CSV files output to `data/` are read by:
- Julia NMG solver: `initialization_from_file()` in `nonlinear_multigrid/julia_multigrid/`
- MATLAB solvers: `ch_initialization.m` with `InputType='file'` in `CHsolvers_package_copy_0206/`
