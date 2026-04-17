# pyCHSolver — Python NMG Solver

Python Implementation of the Nonlinear Multigrid Method for the Cahn-Hilliard Equation.

> **Warning**: This implementation is extremely slow and should not be used for production simulations. It exists for demonstration and educational purposes only.

## Purpose

This is a direct Python port of the C NMG solver (`original_C_NMG/`), implementing the same Full Approximation Scheme (FAS) multigrid algorithm. It is useful for:
- Understanding the algorithm step-by-step in a readable language
- Debugging and unit testing against known C outputs
- Educational demonstrations of the NMG method

## Files

### `pyCHSolver/solver.py`
Complete Python implementation of the NMG solver. Uses global variables with a fixed grid size of 64×64 (`nx=64, ny=64`). Implements the same functions as the C solver:
- `laplace()`, `relax()`, `restrict_ch()`, `prolong_ch()`
- `vcycle()`, `defect()`, `error2()`, `cahn()`

### `pyCHSolver/relax.py`
Unit test comparing the Python `relax()` output against known reference values from the C solver. Used for cross-validation.

### `tests/test_.py`
Test suite for the Python solver.

## Relationship to Other Implementations

This solver is derived from `original_C_NMG/CHsolver.c`. For production use, use the Julia solver in `nonlinear_multigrid/julia_multigrid/` instead.

