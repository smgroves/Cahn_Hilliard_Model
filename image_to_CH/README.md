# Image to Cahn-Hilliard Pipeline

This directory converts **real fluorescence microscopy images of chromosomes** into initial conditions for Cahn-Hilliard simulations. This allows the CH model to be initialized with experimentally observed CPC/protein distributions.

## Pipeline Overview

```
Fluorescence microscopy JPEG
          │
          ▼
   image_to_ch.py          (Python: resize, grayscale, rescale to φ ∈ [-1,1])
   intensity_to_conc.m     (MATLAB equivalent)
          │
          ▼
   input/*.csv             (phi initial condition, 128×128 or 256×256)
          │
          ▼
   run_CH_on_image.m       (MATLAB: calls CahnHilliard_SAV_SMG)
          │
          ▼
   output/*.mp4 / *.png    (simulation results)
```

## Files

### `image_to_ch.py`
Converts a JPEG microscopy image to a CH initial condition CSV.

- **Input:** `T6_I8_Lng-cropped-v2_CPC.jpg` (fluorescence image of chromosome with CPC labeling)
- **Processing:**
  1. Load JPEG
  2. Resize to 256×128 (matching chromosome aspect ratio)
  3. Convert to grayscale
  4. Rescale pixel intensities to φ ∈ [-1, 1]
  5. Pad to square grid
- **Output:** CSV file saved to `input/`

### `intensity_to_conc.m`
MATLAB equivalent of `image_to_ch.py` for the same preprocessing pipeline.

### `run_CH_on_image.m`
Runs the SAV solver on an image-derived initial condition.

- **Input:** Reads phi CSV from `input/`
- **Solver:** Calls `CahnHilliard_SAV_SMG` from `CHsolvers_package_copy_0206/`
- **Domain scaling:** Applies physical domain dimensions:
  - MCF10A cells: 6.6 μm × 3.2 μm
  - HeLa cells: 4.4 μm × 3.2 μm
- **Output:** MP4 movies and PNG snapshots in `output/`

## Input Files

| File | Description |
|------|-------------|
| `T6_I8_Lng-cropped-v2_CPC.jpg` | Source fluorescence microscopy image |
| `MCF10A_T6I8_chromosome_phi_IC_128.csv` | MCF10A IC at 128×128 resolution |
| `MCF10A_T6I8_chromosome_phi_IC_256.csv` | MCF10A IC at 256×256 resolution |
| `MCF10A_T6I8_chromosome_phi_IC_256.csv_C_um_AurB.csv` | CPC Aurora B kinase concentration map |
| `MCF10A_T6I8_chromosome_phi_IC_256.csv_C_um_Borealin.csv` | CPC Borealin concentration map |
| `figure_5_chromosome_phi_IC.csv` | IC used for manuscript Figure 5 |
| `figure_5_chromosome_phi_IC_128.csv` | Figure 5 IC at 128×128 |

## Output

`output/` contains MP4 simulation movies:
- `HeLa_*` — Simulations using HeLa cell parameters
- `MCF10A_*` — Simulations using MCF10A cell parameters
- `MCF10A_Borealin_*` — MCF10A simulations with Borealin as the CPC marker

Each filename encodes: cell line, IC source, solver parameters (epsilon, dt, alpha, domain).

## Biology Context

The CPC (Chromosomal Passenger Complex) includes Aurora B kinase and Borealin, which localize to the inner centromere during mitosis. The fluorescence images capture the spatial distribution of these proteins along chromosome arms. CH simulations model how condensate-like CPC domains form and evolve given the chromosomal geometry.
