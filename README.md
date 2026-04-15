# Apple-TE-CNV

Overview

This repository contains the scripts used for the analysis of transposable elements (TEs) in the context of "Variation in transposable element copy number reveals neutral and environment-associated dynamics across the apple domestication continuum".

The pipeline is designed to:
	•	generate and curate TE libraries
	•	annotate and classify repeats
	•	quantify TE abundance
	•	perform downstream statistical analyses (normalization, PCA, K2P divergence)
	•	produce publication-ready figures

All scripts are organized as a step-by-step workflow corresponding to the analyses presented in the manuscript.

Repository Structure
scripts/
├── s01_RepDeNovo/                  # TE discovery using REPdenovo
├── s02_deduplicate_te_fasta.py    # Remove redundant TE sequences
├── s03_Pastec_annotation/         # TE annotation with PASTEC
├── s04_TESorter_annotation/       # TE classification with TEsorter
├── s05_SortmeRNA/                 # rRNA filtering preparation
├── s06_stampy/                    # Mapping pipeline (SortMeRNA + Stampy)
├── s07_normalisation+PCA/         # Normalization and PCA analysis (R)
├── s08_classify_repeats.R         # TE classification and summarization
├── s09_heatmaps/                  # Heatmap generation
├── s10_panREPET_analysis/         # Advanced TE analysis (separate module)
│   ├── 1_filter_TElibrary.py
│   ├── merge_files.py
│   ├── plot_panREPET.py
│   ├── s01.*_reformat_TE_fasta_db.py
│   ├── s04_calc_k2p_multi.py
│   └── s05_plot_k2p.py
└── s11_climate_adaptation/        # Downstream ecological analysis


Pipeline Description

Step 1 – TE Discovery
	•	Tool: REPdenovo
	•	Script: s01_RepDeNovo
	•	Generates a de novo TE library from sequencing data.


Step 2 – TE Library Curation
	•	Script: s02_deduplicate_te_fasta.py
	•	Removes redundant TE sequences based on similarity thresholds.

Step 3 – TE Annotation
	•	Tool: PASTEC
	•	Script: s03_Pastec_annotation
	•	Annotates TE sequences using structural and homology-based features.


Step 4 – TE Classification
	•	Tool: TEsorter
	•	Script: s04_TESorter_annotation
	•	Classifies TE families based on conserved domains.



Step 5 – Read Mapping
	•	Tools: SortMeRNA + Stampy
	•	Script: s05_stampy
	•	Maps sequencing reads to TE references.


Step 6 – Normalization and PCA
	•	Language: R
	•	Script: s06_normalisation+PCA
	•	Performs normalization of TE counts and principal component analysis.


Step 7 – Repeat Classification Analysis
	•	Script: s07_classify_repeats.R
	•	Aggregates TE classes and computes summary statistics.

Step 8 – Visualization
	•	Script: s08_heatmaps
	•	Generates heatmaps and graphical summaries of TE distributions.


Step 09 – panREPET Analysis (Advanced Module)

This module contains additional scripts for TE library refinement and evolutionary analysis:
	•	TE library filtering (1_filter_TElibrary.py)
	•	File merging (merge_files.py)
	•	Visualization (plot_panREPET.py)
	•	FASTA reformatting (s01.*)
	•	K2P divergence calculation (s04_calc_k2p_multi.py)
	•	K2P plotting (s05_plot_k2p.py)

This module is more complex and is separated for clarity.


Step 10 – Climate Adaptation Analysis
	•	Script: s11_climate_adaptation
	•	Performs downstream ecological/statistical analyses linking TE variation to environmental variables.


