# Finite Element Method (FEM) Solver

This directory contains a prototype **finite element** solver for the Cahn-Hilliard equation using FEniCS/dolfinx.

> **Note**: This implementation is very slow and is not integrated into the rest of the repository pipeline. It exists as a prototype/proof-of-concept.

## File

### `fe_solver.py`

A Python FEM solver based on the Garth N. Wells FEniCS demo (2009), adapted for this project in November 2024.

**Method details:**
- Framework: FEniCS / dolfinx
- Time integration: Crank-Nicolson (second-order)
- Spatial discretization: Finite elements, mixed function spaces
- Domain: Unit square, 96×96 mesh
- Weak form formulation of the 4th-order CH equation split into two second-order equations

**Why it's slow:** FEniCS overhead and the 96×96 mesh make this impractical for the parameter sweeps done elsewhere in the repository.

## Comparison to Other Solvers

| Solver | Grid | Method | Speed |
|--------|------|--------|-------|
| FEM (this) | 96×96 | FEniCS CG elements | Very slow |
| Julia NMG | 128–256 | Multigrid FAS | Fast |
| MATLAB FD | 128–256 | Forward Euler | Slow |

For production simulations, use the Julia NMG solver in `nonlinear_multigrid/julia_multigrid/`.
