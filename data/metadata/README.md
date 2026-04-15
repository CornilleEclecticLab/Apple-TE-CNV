# Sample metadata

This folder contains sample-level metadata used to annotate and subset individuals across analyses (e.g., PCA on SNP dosages and TE counts).

## File

- `sample_metadata.csv`

Each row corresponds to one individual (sample). The primary key is `Sample_ID`.

## Columns

- `Sample_ID`  
  Unique sample identifier (string). This ID is used to match samples across:
  - TE count matrices
  - SNP dosage matrices
  - downstream result tables (e.g., PCA coordinates)

- `Species`  
  Taxonomic label (e.g., *M. domestica*, *M. sylvestris*, *M. sieversii*, *M. orientalis*, *M. baccata*).

- `Geographic_Origin`  
  Country/region of origin. May be `NA` for cultivated accessions.

- `Usage`  
  Biological/agronomic category (e.g., `wild`, `cider`, `dessert`, `rootstock`).

- `Site_Cultivars`  
  Sampling site or cultivar name (free-text). For wild samples, may represent the collection site code.

- `latitude`, `longitude`  
  Sampling coordinates in decimal degrees. Missing values are recorded as `NA`.

- `Code_group_geo_residual`  
  Geographic grouping code used for TE PCA legend/shape-color mapping in the analysis script.
  Examples include: `DOMC`, `DOMD`, `DOMR`, `SYL_DK`, `SYL_FR_EAST`, `SYL_RO`, `SIE_KZ`, `OR_AM`, `BAC_CN`, etc.

- `group_structure_CORE_admix_K6`  
  Genetic cluster assignment (K=6) used to annotate individuals for SNP PCA.
  Values include `Q1`–`Q6` and `Q*_Admix` for admixed assignments.

- `Dataset_source`  
  Source dataset / origin label (e.g., `Chen`, `Duan`). This helps track provenance when combining cohorts.

- `Excluded for TE_structure_residuals analyses`  
  Exclusion flag for specific analyses. This column is intended as a simple marker:
  - blank / `NA`: sample included
  - `x`: sample excluded  
  (Scripts should treat this column as optional and robustly handle missing values.)

## Notes and conventions

- Missing data are represented by `NA`.
- `Sample_ID` must be consistent across files. If other matrices use prefixes (e.g., `A_` or `B_`) or suffixes, cleaning rules may be applied in scripts to match IDs.
- If you modify group codes (`Code_group_geo_residual`) or cluster labels (`group_structure_CORE_admix_K6`), you may need to update the mapping rules used for plotting.

## Recommended usage in scripts

In R:

```r
meta <- readr::read_csv("data/metadata/sample_metadata.csv", show_col_types = FALSE)