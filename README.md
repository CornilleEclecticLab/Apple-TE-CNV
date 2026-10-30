# Apple-TE-CNV

## Overview

This repository contains the scripts used for the analysis of transposable
elements (TEs) in "Transposable element copy number variation reveals
demographic and environment-associated dynamics during apple domestication".

The pipeline is designed to:

- generate and curate TE libraries
- annotate and classify repeats
- quantify TE abundance
- perform downstream statistical analyses (normalisation, PCA, K2P divergence)
- produce publication-ready figures

All scripts are organised as a step-by-step workflow corresponding to the
analyses presented in the manuscript.

## Repository structure

```
config/                                 # Configuration templates (PASTEC, REPdenovo)
data/metadata/                          # Sample-level metadata
analysis/                               # Downstream analyses and the figures (see analysis/README.md)
scripts/
├── s01_RepDeNovo/                      # TE discovery with REPdenovo (+ HTCondor runner)
├── s02_deduplicate_te_fasta.py         # Remove redundant TE sequences
├── s03_Pastec_annotation.sh            # TE annotation with PASTEC
├── s04_TESorter_annotation.sh          # TE classification with TEsorter
├── s05_stampy.pl                       # Read mapping (SortMeRNA + Stampy)
├── s06_normalisation_PCA.R             # Normalisation and PCA
├── s07_classify_repeats.R              # TE classification and summary statistics
├── s08_heatmaps.R                      # Heatmaps and graphical summaries
├── s09_1_launch_panTEannot/            # panREPET: TE library filtering
├── s09_2_panREPET_analysis/            # panREPET: file merging, plots, GO enrichment
├── s09_3_cal_k2p/                      # K2P divergence from BAM pileups
└── s10_climate_adaptation.R            # Downstream ecological analysis
```

`scripts/` is the upstream pipeline, from raw reads to TE copy numbers.
`analysis/` holds everything downstream of that: the population-genomic and
reviewer-response analyses, and the code that draws Figures 1, 2 and 3. See
[analysis/README.md](analysis/README.md).

## Pipeline description

**Step 1 — TE discovery.** `s01_RepDeNovo/`, tool: REPdenovo. Builds a de novo
TE library from the sequencing data.

**Step 2 — TE library curation.** `s02_deduplicate_te_fasta.py`. Removes
redundant TE sequences using identity and coverage thresholds.

**Step 3 — TE annotation.** `s03_Pastec_annotation.sh`, tool: PASTEC.
Annotates TE sequences from structural and homology-based features.

**Step 4 — TE classification.** `s04_TESorter_annotation.sh`, tool: TEsorter.
Classifies TE families from conserved domains.

**Step 5 — Read mapping.** `s05_stampy.pl`, tools: SortMeRNA + Stampy. Maps
the sequencing reads onto the TE references.

**Step 6 — Normalisation and PCA.** `s06_normalisation_PCA.R`. Normalises TE
counts and runs principal component analysis on SNP dosages and TE counts.

**Step 7 — Repeat classification analysis.** `s07_classify_repeats.R`.
Aggregates TE classes and computes summary statistics.

**Step 8 — Visualisation.** `s08_heatmaps.R`. Produces the heatmaps and
graphical summaries of TE distributions.

**Step 9 — panREPET analysis.** Three sub-modules, kept separate because they
are more involved than the rest of the workflow:

- `s09_1_launch_panTEannot/` — TE library filtering (`1_filter_TElibrary.py`)
- `s09_2_panREPET_analysis/` — file merging (`merge_files.py`), plotting and GO
  enrichment (`plot_panREPET.py`)
- `s09_3_cal_k2p/` — FASTA reformatting (`s01.*`), BAM filtering and pileups
  (`s02`, `s03`), K2P divergence (`s04_calc_k2p_multi.py`) and its plots
  (`s05_plot_k2p.py`)

**Step 10 — Climate adaptation analysis.** `s10_climate_adaptation.R`.
Downstream ecological and statistical analyses linking TE variation to
environmental variables.

## Paths and remote hosts

Scripts read their inputs from paths relative to the repository root
(`../../data/...`). Two settings are site-specific and are declared as
variables at the top of the script that uses them, so they can be changed
without touching the rest of the code:

- `REMOTE` and `REMOTE_STAGE` in `scripts/s09_3_cal_k2p/s00_transfer.sh`, the
  host and directory the mapped BAM files are copied from
- `DIR_VS_LISTS` in `scripts/s09_2_panREPET_analysis/plot_panREPET.py`, the
  directory holding the `TE_list_VS_gt_<VS>.txt` files

Several scripts carry SLURM or HTCondor headers; adjust the partition, time
and memory requests to your own scheduler.

## License

This repository is released under the MIT License, an OSI-approved open source
license (https://opensource.org/licenses/MIT). See the LICENSE file for the
full text.

## How to cite

If you use this code, please cite the archived release on Zenodo and the
accompanying publication.
