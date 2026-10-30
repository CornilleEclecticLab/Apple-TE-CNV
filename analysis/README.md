# Downstream analysis — Figures 1 to 3 and the statistics

The `scripts/` folder at the root of this repository covers the upstream
pipeline: building and curating the TE library, annotating and classifying
repeats, mapping the reads and counting TE copies.

This folder covers everything after that: the population-genomic analyses, the
reviewer-response analyses, and the code that draws Figures 1, 2 and 3 of the
manuscript.

## Running order

The scripts are numbered in dependency order and each one writes its tables and
figures under `outputs/<script name>/`. Run them from the folder that holds
`scripts/`, `data/` and `outputs/`, so that the relative paths resolve:

```r
source("scripts/00_setup.R")   # shared paths, seed, session log
```

`00_setup.R` fixes the random seed (20260824) and must be sourced first; every
other script does so itself.

| script | what it produces |
|---|---|
| `01_snp_pca_axes.R` | SNP PCA axes used as covariates |
| `02_candidate_TEs.R` | the candidate TE list behind Figure 3 |
| `03_spatial_vs.R` | the spatial model (Model B) and the vector-strength table |
| `04_elevation_control.R` | altitude control |
| `05`, `12`, `13`, `14` | transcription-factor binding-site enrichment |
| `06`, `11` | genomic context and annotation of the candidates |
| `07_map_figure3A.R` | the arrow geometry of the Figure 3 map |
| `08`, `09`, `10` | per-TE copy-number surfaces and the focal candidate |
| `15`, `19`, `28` | sequence composition, CpG and its deamination products |
| `16`, `27` | windowed genome mosaic and ancestry PCA |
| `17`, `21`, `23`, `26` | successive versions of the assembled Figure 3 |
| `18`, `20`, `22`, `24`, `25` | climate-model diagnostics and permutations |
| `29`, `30` | short-read vs assembly copy number; ComBat sensitivity |
| **`31_figure3_cvd.R`** | **Figure 3 as published, colour-blind-safe palette** |
| **`32_figure3_layers_cvd.R`** | Figure 3 panel A as separate layers, for manual assembly |
| **`33_fig1_heatmaps_cvd.R`** | **Figures 1E and 1F, the TE heatmaps** |
| **`34_palette_reference.R`** | the reference swatch sheet for the whole palette |
| **`35_figS5_snp_pca.R`** | **Supplementary Figure S5, the SNP PCA, as vector output** |
| **`36_figS8_te_pca.R`** | **Supplementary Figure S8, the TE PCA in three dimensions** |
| **`37_figS6_snp_pca_3d.R`** | **Supplementary Figure S6, the SNP PCA in three dimensions** |
| **`38_figS9_te_composition.R`** | **Supplementary Figure S9, structural TE composition** |

## Inputs

These scripts read the TE count matrix, the TE consensus database and the
sample metadata. The two large files are not in this repository; they are
archived on Zenodo at <https://doi.org/10.5281/zenodo.19354151>:

- `TE_counts_clean_wide_with_metadata.csv` — the TE copy-number matrix
- `Final_V1_DUAN_DACCORD_ATIP_DATABASE_non-redundant_150sup.fasta` — the TE
  consensus library
- `250923_TE_structure.txt`, `250923_TE_residuals.txt` — the two TE lists

Place them under `data/` before running, or point the scripts at them.

## Site-specific settings

Two environment variables, both with portable defaults:

- `FIG_OUT` — where the figures are written (default: `../revision_finale/figures_cvd`)
- `TE_DATA_DIR` — where `33_fig1_heatmaps_cvd.R` looks for the count matrix
  (default: `../250923_TE-2`)

## Figure colours

Figures 1, 2 and 3 share a single palette, one colour per category, designed to
stay readable under red-green colour blindness. It contains no green, because
green folds onto yellow under deuteranopia and protanopia and then merges with
orange; red is kept. Every legend was validated with CIEDE2000 computed on
simulated colours. `34_palette_reference.R` draws the full swatch sheet.
