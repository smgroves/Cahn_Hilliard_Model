# Plotting and Analysis

This directory contains all post-processing, analysis, and visualization scripts. Most scripts read simulation outputs from `nonlinear_multigrid/julia_multigrid/manuscript_output/` and `finite_difference_method/output/`.

## Python Scripts

### Energy and Solver Comparison
| Script | Purpose | Reads From |
|--------|---------|------------|
| `plot_energy_mass_across_methods.py` | Compare energy and mass trajectories across FD and NMG solvers for identical ICs | `*_energy.csv`, `*_mass.csv` from solver outputs |
| `plot_MATLAB_vs_Julia_SAV.py` | Compare MATLAB and Julia SAV solver outputs | SAV output directories |
| `plotting_runtimes.py` | Compare wall-clock runtime across all methods | `comp_time_comparison/` CSVs |

### Droplet Analysis
| Script | Purpose | Reads From |
|--------|---------|------------|
| `count_droplets_plotting.py` | Frequency of droplet counts vs. CPC/cohesin parameters | `simulated_droplet_distributions_e_0.0075_noisy_cohesin.csv` |
| `distance_between_droplets_across_params.py` | Statistical analysis of inter-droplet distances (KS tests, KDE) vs. CPC/cohesin params | `distance_between_droplets.csv` |
| `chromosome_lengths.py` | Chromosome arm length statistics for HeLa and MCF10A; violin plots | `count_peaks_imageN_.csv` files in `image_analysis/` |

### Figure Assembly
| Script | Purpose |
|--------|---------|
| `make_CPC_radii_pdf.py` | Assembles simulation snapshot images into PDF |
| `make_kymograph_pdf.py` | Assembles kymograph images into PDF |
| `plot_figure_2.py` | Figure 2 visualization (multigrid snapshots, V-cycle debugging) |
| `plot_CPC_snapshot.py` | CPC geometry snapshots |

### Other
| Script | Purpose |
|--------|---------|
| `plot_ave_mass.py` | Average mass over time for different tolerances |
| `plot_resids.py` | Solver residuals over time |
| `slice_1d_plot.py` | 1D slice plots through phi field |

## MATLAB Scripts

### Kymographs
| Script | Purpose |
|--------|---------|
| `kymograph_central_droplets.m` | Creates kymograph (space-time plot) by slicing phi along central axis. **Inputs:** `indir`, `outdir`, `name`, `dt`, `dtout`, `transposed`, `cutoff`. **Outputs:** PDF kymograph |
| `kymograph_central_droplets_domain.m` | Same but with physical domain scaling |
| `make_kymographs.sh` | Batch kymograph generation across many CPC/cohesin/alpha/epsilon parameter combinations; references Rivanna output directories |

### Droplet Tracking
| Script | Purpose |
|--------|---------|
| `level_set_radius.m` | Track single droplet radius over time via phi=0 contour |
| `level_set_radius_multiple_droplets.m` | Track multiple droplet radii over time |
| `run_level_set_radius_multiple_droplets.sh` | Batch runner for level set radius analysis |
| `find_inflections.m` | Find inflection points in radius-vs-time curves (for critical radius) |
| `plot_inflection_pt.m` | Plot inflection points |
| `count_droplets_over_time.m` | Count number of droplets per timeframe |

### Inter-droplet Distance
| Script | Purpose |
|--------|---------|
| `distance_between_droplets_sim.m` | Compute inter-droplet distances from simulation phi fields |
| `distance_between_droplets_sim_chr_lengths.m` | Same, with chromosome length data; reads from `domain_0_2_from_rivanna_kymographs_e_0.0075` |
| `distance_between_droplets_sim_chr_lengths_256.m` | 256×256 grid version |

### Visualization Helpers
| Script | Purpose |
|--------|---------|
| `CHplotting_function.m` | Standard CH field visualization |
| `CHplotting_gif.m` | GIF output version |
| `CHplotting_legacy.m` | Older plotting function |
| `redbluecmap.csv` | Custom red-blue colormap data |

## Subdirectories

| Directory | Contents |
|-----------|----------|
| `comparison_sim_results_across_methods/` | PDF/PNG comparison figures: FD vs. NMG for spinodal, checkerboard, two-droplet |
| `comp_time_comparison/` | Computation time CSVs and PDFs comparing Python, C, Julia |
| `image_analysis/` | `April_28.ipynb` (chromosome image analysis notebook), chromosome length CSVs, histogram PDFs |
| `simulated_droplet_distributions/` | Droplet count frequency data vs. CPC/cohesin parameters |
| `distance_between_droplets/` | Inter-droplet distance analysis data |
| `kymograph_pdfs/` | Generated kymograph PDFs |
| `radii_over_time_pdfs/` | Droplet radius vs. time PDFs |
| `radii_lineplots_kymographs/` | Combined radius line plots and kymographs |
| `critical_radius_timepoint_plots/` | Snapshots at critical radius timepoints |
| `manuscript/` | Final manuscript figures |

## Rivanna Data

Several analysis scripts specifically read from directories containing Rivanna HPC results:
- `distance_between_droplets_sim_chr_lengths.m` reads from `domain_0_2_from_rivanna_kymographs_e_0.0075`
- `make_kymographs.sh` references multiple Rivanna job array outputs
- `image_analysis/April_28.ipynb` references `from_rivanna_kymographs_e_0.0075`

These directories were populated using `utils/combine_all_folders.py` to merge batches of Rivanna output.
