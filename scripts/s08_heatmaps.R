#!/usr/bin/env Rscript

# ============================================================
# TE Heatmaps (group means, log1p, scaled per TE)
# - Builds two heatmaps: STRUCTURE and RESIDUALS
# - Uses TE lists (one per heatmap)
# - Annotates TE columns with: Class (I/II), Family, Superfamily
#   * Superfamily is only Gypsy/Copia for LTR, otherwise Unknown
# - Output: PDF + PNG per heatmap
#
# Inputs:
#   1) TE counts matrix: data_frame_countTE_cleaned.csv
#      - First column must be Sample IDs (default: "V1")
#   2) Metadata: 250923_sample_list_TE.csv
#      - Must include Sample_ID
#      - Must include structure cluster column (default name below)
#      - Must include residual group column (default name below)
#   3) TE list files: 250923_TE_structure.txt / 250923_TE_residuals.txt
#
# Example:
#   Rscript scripts/te_heatmaps.R \
#     --counts data_frame_countTE_cleaned.csv \
#     --meta 250923_sample_list_TE.csv \
#     --te_structure 250923_TE_structure.txt \
#     --te_residuals 250923_TE_residuals.txt \
#     --outdir results_heatmaps_TE
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(pheatmap)
  library(optparse)
})

options(stringsAsFactors = FALSE)

# -----------------------------
# CLI options
# -----------------------------
option_list <- list(
  make_option(c("--counts"), type="character", default="data_frame_countTE_cleaned.csv",
              help="TE count matrix (CSV). First column = sample IDs [default: %default]"),
  make_option(c("--meta"), type="character", default="250923_sample_list_TE.csv",
              help="Sample metadata table (CSV/TSV) [default: %default]"),
  make_option(c("--te_structure"), type="character", default="250923_TE_structure.txt",
              help="TE IDs list for STRUCTURE heatmap [default: %default]"),
  make_option(c("--te_residuals"), type="character", default="250923_TE_residuals.txt",
              help="TE IDs list for RESIDUALS heatmap [default: %default]"),
  make_option(c("--outdir"), type="character", default="results_heatmaps_TE",
              help="Output directory [default: %default]"),
  make_option(c("--id_col"), type="character", default="V1",
              help="Sample ID column name in counts matrix [default: %default]"),
  make_option(c("--structure_col"), type="character", default="group_structure_CORE_admix_(K=6)",
              help="Column in meta for structure clusters [default: %default]"),
  make_option(c("--residual_col"), type="character", default="Code_group_geo_residual",
              help="Column in meta for residual geo groups [default: %default]"),
  make_option(c("--min_group_n"), type="integer", default=3,
              help="Minimum samples per residual group to keep [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# -----------------------------
# Helpers
# -----------------------------
stop_if_missing <- function(x) {
  if (!file.exists(x)) stop("Missing input file: ", x)
}

msg <- function(...) cat("[INFO] ", sprintf(...), "\n", sep = "")

# Class/Family/Superfamily parser
# - Class: I for LTR/LINE/SINE, II for TIR/Helitron/MITE, else Unknown
# - Family: LTR/LINE/SINE/TIR/MITE/Helitron/SSR/pararetrovirus/mixture/noCat/Unknown
# - Superfamily: only Gypsy/Copia if present; else Unknown
parse_te_annotation <- function(te_names) {
  te_low <- str_to_lower(te_names)

  class <- case_when(
    str_detect(te_low, "ltr|line|sine") ~ "Class I",
    str_detect(te_low, "tir|helitron|mite") ~ "Class II",
    TRUE ~ "Unknown"
  )

  family <- case_when(
    str_detect(te_low, "mite") ~ "MITE",
    str_detect(te_low, "tir") ~ "TIR",
    str_detect(te_low, "ltr") ~ "LTR",
    str_detect(te_low, "line") ~ "LINE",
    str_detect(te_low, "sine") ~ "SINE",
    str_detect(te_low, "helitron") ~ "Helitron",
    str_detect(te_low, "ssr") ~ "SSR",
    str_detect(te_low, "pararetrovirus") ~ "pararetrovirus",
    str_detect(te_low, "mixture") ~ "mixture",
    str_detect(te_low, "nocat") ~ "noCat",
    TRUE ~ "Unknown"
  )

  # superfamily ONLY Gypsy/Copia, otherwise Unknown
  superfamily <- case_when(
    str_detect(te_low, "gypsy") ~ "Gypsy",
    str_detect(te_low, "copia") ~ "Copia",
    TRUE ~ "Unknown"
  )

  out <- data.frame(
    TE = te_names,
    Class = factor(class, levels = c("Class I", "Class II", "Unknown")),
    Family = factor(family, levels = c("LTR","LINE","SINE","TIR","MITE","Helitron","SSR","pararetrovirus","mixture","noCat","Unknown")),
    Superfamily = factor(superfamily, levels = c("Gypsy","Copia","Unknown")),
    stringsAsFactors = FALSE
  )
  out
}

# Build group mean matrix (log1p), rows = groups, cols = TE
make_group_means <- function(count_mat, meta_df, group_col, te_list) {
  te_found <- intersect(colnames(count_mat), te_list)
  if (length(te_found) == 0) stop("None of the TE IDs from the list were found in the count matrix columns.")

  # align samples
  common_ids <- intersect(rownames(count_mat), meta_df$Sample_ID)
  if (length(common_ids) == 0) stop("No shared Sample_ID between counts and metadata.")

  meta_df <- meta_df %>% filter(Sample_ID %in% common_ids)
  count_sub <- count_mat[meta_df$Sample_ID, te_found, drop = FALSE]

  # log transform
  count_sub <- log1p(count_sub)

  # group means
  df_means <- as.data.frame(count_sub) %>%
    mutate(.group = meta_df[[group_col]]) %>%
    group_by(.group) %>%
    summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")

  mat <- as.data.frame(df_means)
  rownames(mat) <- mat$.group
  mat$.group <- NULL

  mat
}

# Colors for annotations (stable + readable)
annotation_colors <- list(
  Class = c("Class I"="#E64B35FF", "Class II"="#4DBBD5FF", "Unknown"="grey70"),
  Family = c(
    "LTR"="#F39C12", "LINE"="#F7DC6F", "SINE"="#2ECC71", "TIR"="#8E44AD",
    "MITE"="#6C5CE7", "Helitron"="#F78FB3", "SSR"="#8D6E63",
    "pararetrovirus"="#00B8D9", "mixture"="black", "noCat"="grey50", "Unknown"="grey70"
  ),
  Superfamily = c("Gypsy"="#5DADE2", "Copia"="#58D68D", "Unknown"="grey70")
)

plot_heatmap <- function(mat, out_prefix) {
  te_annot <- parse_te_annotation(colnames(mat))
  annotation_col <- te_annot %>% select(Class, Family, Superfamily)
  rownames(annotation_col) <- te_annot$TE

  # pheatmap writes files only if filename is provided
  pheatmap(
    as.matrix(mat),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "correlation",
    clustering_method = "ward.D2",
    scale = "column",                 # scale each TE across groups
    show_colnames = FALSE,
    fontsize_row = 10,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    main = out_prefix,
    filename = paste0(out_prefix, ".pdf"),
    width = 10,
    height = 7
  )

  # Also save as PNG (high-res)
  png(paste0(out_prefix, ".png"), width = 2400, height = 1800, res = 300)
  pheatmap(
    as.matrix(mat),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "correlation",
    clustering_method = "ward.D2",
    scale = "column",
    show_colnames = FALSE,
    fontsize_row = 10,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    main = out_prefix
  )
  dev.off()
}

# -----------------------------
# Main
# -----------------------------
stop_if_missing(opt$counts)
stop_if_missing(opt$meta)
stop_if_missing(opt$te_structure)
stop_if_missing(opt$te_residuals)

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

msg("Reading TE counts: %s", opt$counts)
counts_dt <- fread(opt$counts)

if (!opt$id_col %in% names(counts_dt)) {
  stop("Counts matrix does not contain id_col = '", opt$id_col, "'. Available columns: ",
       paste(names(counts_dt), collapse = ", "))
}

# Fix one known naming artifact if present (optional; harmless otherwise)
names(counts_dt) <- str_replace(names(counts_dt), "LTR_PIF_Harbinger", "LTR")

# to matrix with rownames
counts_df <- as.data.frame(counts_dt)
rownames(counts_df) <- counts_df[[opt$id_col]]
counts_df[[opt$id_col]] <- NULL

# ensure numeric
counts_df[] <- lapply(counts_df, function(x) as.numeric(as.character(x)))
if (anyNA(counts_df)) {
  msg("Warning: NA introduced during numeric conversion (this can be OK if input had missing).")
}

msg("Reading metadata: %s", opt$meta)
meta <- fread(opt$meta) %>% as.data.frame()

if (!"Sample_ID" %in% names(meta)) stop("Metadata must contain a 'Sample_ID' column.")
if (!opt$structure_col %in% names(meta)) stop("Metadata missing structure_col: ", opt$structure_col)
if (!opt$residual_col %in% names(meta)) stop("Metadata missing residual_col: ", opt$residual_col)

# TE lists
te_structure <- readLines(opt$te_structure)
te_residuals <- readLines(opt$te_residuals)

# STRUCTURE groups: remove Admix
meta_structure <- meta %>%
  transmute(Sample_ID = Sample_ID, group = .data[[opt$structure_col]]) %>%
  filter(!grepl("Admix", group))

# Rename Q1..Q6 (if present)
rename_q <- c("Q1"="SYLW","Q2"="DOMD","Q3"="BACC","Q4"="SIOR","Q5"="SYLE","Q6"="DOMC")
meta_structure$group <- ifelse(meta_structure$group %in% names(rename_q),
                              rename_q[meta_structure$group], meta_structure$group)

# RESIDUALS groups: keep groups with >= min_group_n
meta_residuals <- meta %>%
  transmute(Sample_ID = Sample_ID, group = .data[[opt$residual_col]]) %>%
  group_by(group) %>% filter(n() >= opt$min_group_n) %>% ungroup()

# Build matrices
msg("Building STRUCTURE group means matrix...")
mat_structure <- make_group_means(counts_df, meta_structure, "group", te_structure)

msg("Building RESIDUALS group means matrix...")
mat_residuals <- make_group_means(counts_df, meta_residuals, "group", te_residuals)

# Write the group-mean tables too (often useful for supplements)
fwrite(as.data.table(cbind(Group = rownames(mat_structure), mat_structure)),
       file = file.path(opt$outdir, "group_means_STRUCTURE_log1p.tsv"),
       sep = "\t")
fwrite(as.data.table(cbind(Group = rownames(mat_residuals), mat_residuals)),
       file = file.path(opt$outdir, "group_means_RESIDUALS_log1p.tsv"),
       sep = "\t")

# Plot heatmaps
msg("Plotting heatmaps into: %s", opt$outdir)
plot_heatmap(mat_structure, file.path(opt$outdir, "heatmap_STRUCTURE"))
plot_heatmap(mat_residuals, file.path(opt$outdir, "heatmap_RESIDUALS"))

msg("Done. Outputs written in: %s", opt$outdir)