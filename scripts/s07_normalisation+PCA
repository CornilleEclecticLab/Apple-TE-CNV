#!/usr/bin/env Rscript

# ============================================================
# PCA on SNP dosages and TE counts
# - Cleans sample IDs across sources
# - Intersects samples between TE counts / sample metadata / SNP dosages
# - Runs PCA (3D) and exports plots + coordinates
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(tibble)
  library(Biostrings)
  library(scatterplot3d)
})

# -----------------------------
# User parameters (edit these)
# -----------------------------
INPUT_SAMPLE_INFO <- "250923_sample_list_TE.csv"
INPUT_TE_COUNTS   <- "data_frame_all_TE 2.csv"
INPUT_DOSAGE_RDA  <- "dosage_data.rda"  # loads object 'dosage_data'
INPUT_TE_FASTA    <- "Final_V1_DUAN_DACCORD_ATIP_DATABASE_non-redundant_150sup.fasta"

OUTDIR <- "results/pca"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helpers
# -----------------------------
message_info <- function(...) message("[INFO] ", ...)
message_ok   <- function(...) message("[OK] ", ...)
message_warn <- function(...) warning("[WARN] ", ..., call. = FALSE)

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing input file: ", path)
}

export_scatterplot3d <- function(x, groups, pch, col, variance_explained, title, pdf_file, svg_file) {
  plot_fun <- function() {
    scatterplot3d(
      x[, 1:3],
      angle = 50,
      pch = pch,
      color = col,
      cex.symbols = 1.5,
      main = title,
      xlab = paste0("PC1 (", variance_explained[1], "%)"),
      ylab = paste0("PC2 (", variance_explained[2], "%)"),
      zlab = paste0("PC3 (", variance_explained[3], "%)")
    )
    legend("topright", legend = groups$legend_labels, pch = groups$legend_pch, col = groups$legend_col, cex = 0.85)
  }

  pdf(pdf_file, width = 7, height = 7)
  plot_fun()
  dev.off()

  svg(svg_file, width = 7, height = 7)
  plot_fun()
  dev.off()
}

remove_zero_variance_cols <- function(mat, out_txt = NULL) {
  variances <- apply(mat, 2, var, na.rm = TRUE)
  zero_var_cols <- names(variances[variances == 0])

  if (length(zero_var_cols) > 0) {
    if (!is.null(out_txt)) writeLines(zero_var_cols, out_txt)
    mat <- mat[, !(colnames(mat) %in% zero_var_cols), drop = FALSE]
    message_info("Removed zero-variance columns (n=", length(zero_var_cols), ").")
  } else {
    message_info("No zero-variance columns found.")
  }
  mat
}

# -----------------------------
# 0) Input checks
# -----------------------------
stop_if_missing(INPUT_SAMPLE_INFO)
stop_if_missing(INPUT_TE_COUNTS)
stop_if_missing(INPUT_DOSAGE_RDA)
stop_if_missing(INPUT_TE_FASTA)

# -----------------------------
# 1) Load inputs
# -----------------------------
message_info("Loading inputs...")
sample_info  <- readr::read_csv(INPUT_SAMPLE_INFO, show_col_types = FALSE)
count_te_raw <- readr::read_csv(INPUT_TE_COUNTS, show_col_types = FALSE)
load(INPUT_DOSAGE_RDA)  # expects dosage_data
fasta_te <- Biostrings::readDNAStringSet(INPUT_TE_FASTA)

if (!exists("dosage_data")) {
  stop("Object 'dosage_data' not found in ", INPUT_DOSAGE_RDA)
}

# -----------------------------
# 2) Clean TE count sample IDs
# -----------------------------
count_te_clean <- count_te_raw %>%
  as_tibble() %>%
  rename(Sample_ID_raw = 1) %>%
  mutate(
    Sample_ID = Sample_ID_raw %>%
      gsub("^A_|^B_", "", .) %>%
      gsub("_?list_sorting\\.txt$", "", .)
  ) %>%
  select(-Sample_ID_raw) %>%
  select(Sample_ID, everything()) %>%
  column_to_rownames("Sample_ID")

message_info("Example TE sample IDs after cleaning:")
print(head(rownames(count_te_clean), 10))

# -----------------------------
# 3) Clean sample_info
# -----------------------------
if (!"Sample_ID" %in% names(sample_info)) {
  sample_info <- sample_info %>% rename(Sample_ID = 1)
}
message_info("Example sample_info IDs:")
print(head(sample_info$Sample_ID, 10))

# -----------------------------
# 4) Intersect samples across sources
# -----------------------------
samples_te   <- rownames(count_te_clean)
samples_info <- sample_info$Sample_ID
samples_snp  <- rownames(dosage_data)

common_samples <- Reduce(intersect, list(samples_te, samples_info, samples_snp))
message_ok("Number of samples shared by TE counts, sample_info and SNP dosages: ", length(common_samples))

if (length(common_samples) == 0) {
  stop("No common samples found. Check ID cleaning rules and input files.")
}

# -----------------------------
# 5) Filter objects to common samples
# -----------------------------
count_te_final <- count_te_clean[common_samples, , drop = FALSE]
sample_info_final <- sample_info %>% filter(Sample_ID %in% common_samples)
dosage_final <- dosage_data[common_samples, , drop = FALSE]

# Keep sample order consistent
sample_info_final <- sample_info_final %>% arrange(match(Sample_ID, common_samples))
count_te_final <- count_te_final[match(sample_info_final$Sample_ID, rownames(count_te_final)), , drop = FALSE]
dosage_final   <- dosage_final[match(sample_info_final$Sample_ID, rownames(dosage_final)), , drop = FALSE]

stopifnot(all(rownames(count_te_final) == sample_info_final$Sample_ID))
stopifnot(all(rownames(dosage_final) == sample_info_final$Sample_ID))

# -----------------------------
# 6) Assign TE names from FASTA
# -----------------------------
sequence_names <- names(fasta_te)
if (length(sequence_names) == ncol(count_te_final)) {
  colnames(count_te_final) <- sequence_names
  message_ok("TE matrix column names updated from FASTA sequence headers.")
} else {
  message_warn("TE count columns (", ncol(count_te_final), ") != FASTA sequences (", length(sequence_names), "). Keeping original column names.")
}

# ============================================================
# PART A — PCA on SNP dosages
# ============================================================
message_info("Running PCA on SNP dosages...")

# Recode group variable (requires this column to exist)
if (!"group_structure_CORE_admix_K6" %in% names(sample_info_final)) {
  stop("Missing column 'group_structure_CORE_admix_K6' in sample_info.")
}

group_recoded <- dplyr::recode(sample_info_final$group_structure_CORE_admix_K6,
  "Q1"       = "SYLW",
  "Q2"       = "DOMD",
  "Q3"       = "BACC",
  "Q4"       = "SIOR",
  "Q5"       = "SYLE",
  "Q6"       = "DOMC",
  "Q1_Admix" = "SYLW_Admix",
  "Q2_Admix" = "DOMD_Admix",
  "Q3_Admix" = "BACC_Admix",
  "Q4_Admix" = "SIOR_Admix",
  "Q5_Admix" = "SYLE_Admix",
  "Q6_Admix" = "DOMC_Admix"
)

group_levels <- c("SYLW", "DOMD", "BACC", "SIOR", "SYLE", "DOMC",
                  "SYLW_Admix", "DOMD_Admix", "BACC_Admix", "SIOR_Admix", "SYLE_Admix", "DOMC_Admix")

group_colors <- c("#3171ad", "#dca327", "#6eb2e4", "#c66526", "#eee461", "#c17da5",
                  "#3171ad", "#dca327", "#6eb2e4", "#c66526", "#eee461", "#c17da5")

col_snp <- group_colors[match(group_recoded, group_levels)]
pch_snp <- ifelse(grepl("Admix", group_recoded), 1, 16)

# Remove zero-variance SNPs
dosage_filt <- remove_zero_variance_cols(dosage_final)

pca_snp <- prcomp(dosage_filt, scale. = TRUE)
var_snp <- round((pca_snp$sdev^2 / sum(pca_snp$sdev^2)) * 100, 1)

unique_groups <- sort(unique(group_recoded))
legend_pch  <- ifelse(grepl("Admix", unique_groups), 1, 16)
legend_cols <- group_colors[match(unique_groups, group_levels)]

groups_snp <- list(
  legend_labels = unique_groups,
  legend_pch = legend_pch,
  legend_col = legend_cols
)

export_scatterplot3d(
  x = pca_snp$x,
  groups = groups_snp,
  pch = pch_snp,
  col = col_snp,
  variance_explained = var_snp,
  title = "SNP PCA (3D)",
  pdf_file = file.path(OUTDIR, "PCA_SNP_plot.pdf"),
  svg_file = file.path(OUTDIR, "PCA_SNP_plot.svg")
)

pca_snp_coords <- as.data.frame(pca_snp$x[, 1:3])
pca_snp_coords$Sample_ID <- rownames(pca_snp_coords)
pca_snp_coords <- pca_snp_coords %>%
  left_join(sample_info_final %>% select(Sample_ID, group_structure_CORE_admix_K6), by = "Sample_ID") %>%
  select(Sample_ID, group_structure_CORE_admix_K6, PC1, PC2, PC3)

readr::write_csv(pca_snp_coords, file.path(OUTDIR, "PCA_SNP_coordinates.csv"))
readr::write_delim(pca_snp_coords, file.path(OUTDIR, "PCA_SNP_coordinates.tsv"), delim = "\t")
message_ok("SNP PCA outputs written to: ", OUTDIR)

# ============================================================
# PART B — PCA on TE counts
# ============================================================
message_info("Running PCA on TE counts...")

if (!"Code_group_geo_residual" %in% names(sample_info_final)) {
  stop("Missing column 'Code_group_geo_residual' in sample_info.")
}

get_color <- function(code) {
  if (code == "DOMD") return("#8b0000")
  else if (code == "DOMC") return("#ff0000")
  else if (code == "DOMR") return("#f08080")
  else if (code == "BAC_CN") return("#000000")
  else if (code == "BAC_FSU") return("#bebebe")
  else if (code == "SIE_CN") return("#a9746e")
  else if (code == "SIE_KZ") return("#5d4037")
  else if (code == "OR_AM") return("#7e57c2")
  else if (code == "OR_TK") return("#512da8")
  else if (code == "SYL_DK") return("#1d8348")
  else if (code == "SYL_AT") return("#28b463")
  else if (code == "SYL_DE") return("#2ecc71")
  else if (code == "SYL_FR_EAST") return("#abebc6")
  else if (code == "SYL_FR_MID") return("#d5f5e3")
  else if (code == "SYL_FR_WEST") return("#eafbf2")
  else if (code == "SYL_IT") return("#82e0aa")
  else if (code == "SYL_MKD") return("#58d68d")
  else if (code == "SYL_RO") return("#239b56")
  else if (code == "SYL_NW") return("#3b6440")
  else return("black")
}

get_shape <- function(code) {
  if (grepl("^DOM", code)) return(11)
  prefix <- ifelse(grepl("_", code), sub("_.*", "", code), code)
  if (prefix == "BAC") return(17)
  else if (prefix == "SIE") return(15)
  else if (prefix == "OR")  return(16)
  else if (prefix == "SYL") return(18)
  else return(1)
}

colors_te <- sapply(sample_info_final$Code_group_geo_residual, get_color)
shapes_te <- sapply(sample_info_final$Code_group_geo_residual, get_shape)

count_te_filtered <- remove_zero_variance_cols(
  count_te_final,
  out_txt = file.path(OUTDIR, "TE_zero_variance_removed.txt")
)

pca_te <- prcomp(count_te_filtered, scale. = TRUE)
var_te <- round((pca_te$sdev^2 / sum(pca_te$sdev^2)) * 100, 1)

unique_codes <- sort(unique(sample_info_final$Code_group_geo_residual))
legend_colors <- sapply(unique_codes, get_color)
legend_shapes <- sapply(unique_codes, get_shape)

groups_te <- list(
  legend_labels = unique_codes,
  legend_pch = legend_shapes,
  legend_col = legend_colors
)

export_scatterplot3d(
  x = pca_te$x,
  groups = groups_te,
  pch = shapes_te,
  col = colors_te,
  variance_explained = var_te,
  title = "TE PCA (3D, geography-coded)",
  pdf_file = file.path(OUTDIR, "PCA_TE_plot.pdf"),
  svg_file = file.path(OUTDIR, "PCA_TE_plot.svg")
)

pca_te_coords <- as.data.frame(pca_te$x[, 1:3])
pca_te_coords$Sample_ID <- rownames(pca_te_coords)
pca_te_coords <- pca_te_coords %>%
  left_join(sample_info_final %>% select(Sample_ID, Code_group_geo_residual), by = "Sample_ID") %>%
  select(Sample_ID, Code_group_geo_residual, PC1, PC2, PC3)

readr::write_csv(pca_te_coords, file.path(OUTDIR, "PCA_TE_coordinates.csv"))
readr::write_delim(pca_te_coords, file.path(OUTDIR, "PCA_TE_coordinates.tsv"), delim = "\t")
message_ok("TE PCA outputs written to: ", OUTDIR)

message_ok("Done.")