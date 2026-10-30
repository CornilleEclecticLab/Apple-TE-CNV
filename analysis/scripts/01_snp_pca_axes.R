# =============================================================================
# Script:  01_snp_pca_axes.R
# Purpose: Follow-up to 02_candidate_TEs.R (reviewer 2, comment 5),
#          answering three questions raised after the first pass:
#
#          1. Why only three SNP-based principal components? How much extra
#             genetic variance is captured by going to 5 or 10 axes?
#          2. Model A should also be re-run WITH the PC covariates in each
#             M. sylvestris population separately (SYLW = West, SYLE = East)
#             as well as in COMBINED, not only in COMBINED.
#          3. How many candidate TEs does each combination give, and what is
#             shared versus specific between them, at 0.05 / 0.01 / 0.001?
#
#          The published SNP PCA exported only PC1-3, so the same PCA is
#          reproduced here from the same input (dosage_data.rda, prcomp with
#          scale. = TRUE, zero-variance columns dropped) in order to obtain
#          PC1-10 and the eigenvalues. The reproduction is validated against
#          the saved PC1-3 coordinates and the script aborts if they do not
#          match, so this is an extension of the existing PCA, not a new one.
#
#          Model grid: 3 datasets x 4 PC depths (0, 3, 5, 10) x 19 bioclimatic
#          variables x all tested residual TEs, quasi-Poisson, exactly as in
#          Model A. Combinations with too few residual degrees of freedom
#          (notably SYLE, n = 10) are reported as infeasible and skipped
#          rather than fitted unreliably.
#
#          Stringency: two multiple-testing regimes are reported side by side,
#          BH within each bioclimatic variable (what Model A intended) and BH
#          across all TE x variable tests at once (stricter).
#
# Author:  Antho (analysis run with Claude Code)
# Date:    2026-08-24
# Project: TE_Project / review (reviewer response, Genome Biology submission)
#
# Inputs:
#   ../1-PCA/dosage_data.rda            SNP dosage matrix (source of the PCA)
#   data/PCA_SNP_coordinates.csv        published PC1-3, used for validation
#   data/te_geo_SYL_individus_long_with_coords.csv
#   data/file_names_clusters_ok.csv
#   data/worldclim_env_data/climate/wc2.1_10m/wc2.1_10m_bio_*.tif
#
# Outputs (all under outputs/01_snp_pca_axes/):
#   tables/pca_variance_explained.csv       scree / cumulative variance
#   tables/PCA_SNP_coordinates_PC1_10.csv   extended PC coordinates
#   tables/pca_reproduction_check.csv       validation against published PC1-3
#   tables/model_grid_feasibility.csv       which subset x depth were fitted
#   tables/glm_<SUBSET>_pc<DEPTH>.csv       full GLM results per combination
#   tables/candidate_set_sizes_by_grid.csv  n candidate TEs per cell
#   tables/overlap_across_subsets.csv       shared vs specific between datasets
#   tables/overlap_across_depths.csv        effect of adding PC axes
#   tables/jaccard_pairwise.csv             all pairwise Jaccard indices
#   lists/*.txt                             plain TE-ID lists
#   figures/*.pdf + *.svg
#   ../../logs/01_snp_pca_axes_*_sessionInfo.txt
#
# NAMING NOTE: in the manuscript, "Model A" is the climate GLM
# (TE count ~ bioclimatic variable) and "Model B" is the SPATIAL model
# (TE count ~ latitude + longitude, negative binomial), which produces the
# vector-strength (VS) values shown in Fig. 3A-B. To avoid a clash with that
# established terminology, the PC-covariate version of the climate model is
# called "Model A-PC" throughout this folder, never "Model B".
# =============================================================================

# --- Setup -------------------------------------------------------------------
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(terra)     # bioclim raster extraction
  library(svglite)   # SVG device for ggsave
})

set.seed(20260824)

# --- Parameters --------------------------------------------------------------
PC_DEPTHS       <- c(0, 3, 5, 10)   # number of SNP PCs used as covariates
N_PC_EXPORT     <- 10               # PC axes to export
THRESHOLDS      <- c(0.05, 0.01, 0.001)
PADJ_METHOD     <- "BH"
MIN_OBS_PER_TE  <- 10
MIN_RESIDUAL_DF <- 5                # below this a subset x depth cell is skipped

DOSAGE_RDA <- file.path("..", "1-PCA", "dosage_data.rda")

# Dataset definitions, following Model A's fastSTRUCTURE cluster rule.
SUBSET_CLUSTERS <- list(
  SYLW     = c("Q1"),
  SYLE     = c("Q5"),
  COMBINED = c("Q1", "Q1_Admix", "Q5", "Q5_Admix")
)

dir_out  <- file.path(dir_outputs, "01_snp_pca_axes")
dir_tab  <- file.path(dir_out, "tables")
dir_fig  <- file.path(dir_out, "figures")
dir_list <- file.path(dir_out, "lists")
for (d in c(dir_out, dir_tab, dir_fig, dir_list)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# Vector formats only (PDF + SVG) for downstream figure assembly.
save_plot <- function(p, name, width = 7, height = 5) {
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".pdf")), p,
                  width = width, height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".svg")), p,
                  width = width, height = height, device = svglite::svglite)
}

message("=== 01_snp_pca_axes.R ===")

# =============================================================================
# 1. Reproduce the SNP PCA and extend it to PC1-10
# =============================================================================
# The published PCA used prcomp(scale. = TRUE) on the SNP dosage matrix
# restricted to the 189 accessions, after dropping zero-variance columns.
# The manuscript credits FactoMineR; the code that produced the exported
# coordinates is prcomp. Reproducing that exact call is what allows PC4-10 to
# be added without changing the axes already used.
pca_published <- readr::read_csv(file.path(dir_data, "PCA_SNP_coordinates.csv"),
                                 show_col_types = FALSE)
pca_ids <- as.character(pca_published$Sample_ID)

stopifnot(file.exists(DOSAGE_RDA))
dosage_env <- new.env()
load(DOSAGE_RDA, envir = dosage_env)
dosage_data <- get("dosage_data", envir = dosage_env)
dosage_mat  <- as.matrix(dosage_data)

missing_dosage <- setdiff(pca_ids, rownames(dosage_mat))
if (length(missing_dosage) > 0) {
  stop("Accessions missing from the dosage matrix: ",
       paste(missing_dosage, collapse = ", "))
}

# Same order as the published coordinate file, then the same column filters.
pca_input <- dosage_mat[pca_ids, , drop = FALSE]
pca_input <- pca_input[, colSums(is.na(pca_input)) == 0, drop = FALSE]
col_var   <- apply(pca_input, 2, stats::var)
pca_input <- pca_input[, col_var > 0, drop = FALSE]

message(sprintf("PCA input: %d accessions x %d SNP columns",
                nrow(pca_input), ncol(pca_input)))

pca_fit <- stats::prcomp(pca_input, scale. = TRUE)

var_explained <- pca_fit$sdev^2 / sum(pca_fit$sdev^2) * 100
variance_tbl <- tibble::tibble(
  PC                  = seq_along(var_explained),
  variance_explained  = var_explained,
  cumulative_variance = cumsum(var_explained)
)
readr::write_csv(utils::head(variance_tbl, 50),
                 file.path(dir_tab, "pca_variance_explained.csv"))

# --- Validation: PC1-3 must match the published coordinates ------------------
repro_check <- purrr::map_dfr(1:3, function(i) {
  r <- stats::cor(pca_fit$x[, i], pca_published[[paste0("PC", i)]])
  tibble::tibble(PC = i, pearson_r = r,
                 max_abs_diff = max(abs(pca_fit$x[, i] -
                                          pca_published[[paste0("PC", i)]])))
})
readr::write_csv(repro_check, file.path(dir_tab, "pca_reproduction_check.csv"))
message("PCA reproduction check against the published PC1-3:")
print(as.data.frame(repro_check))

if (any(abs(repro_check$pearson_r) < 0.999)) {
  stop("Reproduced PCA does not match the published PC1-3; aborting rather ",
       "than silently using different axes.")
}

pc_coords <- tibble::as_tibble(pca_fit$x[, seq_len(N_PC_EXPORT), drop = FALSE])
pc_coords$sample_id <- toupper(pca_ids)
pc_coords <- pc_coords %>% dplyr::relocate(sample_id)
readr::write_csv(pc_coords, file.path(dir_tab, "PCA_SNP_coordinates_PC1_10.csv"))

message("Variance explained by the SNP PCA:")
for (k in c(2, 3, 5, 10)) {
  message(sprintf("  PC1-%-2d : %5.2f %% cumulative", k,
                  variance_tbl$cumulative_variance[k]))
}
message(sprintf("  going from PC1-3 to PC1-10 adds %.2f percentage points for 7 extra df",
                variance_tbl$cumulative_variance[10] -
                  variance_tbl$cumulative_variance[3]))

# =============================================================================
# 2. TE counts, bioclim and dataset definitions
# =============================================================================
te_long <- readr::read_csv(
  file.path(dir_data, "te_geo_SYL_individus_long_with_coords.csv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    sample_id = toupper(as.character(sample_id)),
    latitude  = as.numeric(x),   # "x" holds latitude in the source file
    longitude = as.numeric(y)    # "y" holds longitude
  )

coords_unique <- te_long %>% dplyr::distinct(sample_id, latitude, longitude)

clim_dir   <- file.path(dir_data, "worldclim_env_data", "climate", "wc2.1_10m")
clim_files <- sprintf("wc2.1_10m_bio_%d.tif", 1:19)
stopifnot(all(file.exists(file.path(clim_dir, clim_files))))

bio_stack <- terra::rast(file.path(clim_dir, clim_files))
names(bio_stack) <- sprintf("bio%d", 1:19)
env_values <- terra::extract(
  bio_stack, as.matrix(coords_unique[, c("longitude", "latitude")])
)
env_full <- dplyr::bind_cols(coords_unique, tibble::as_tibble(env_values))

clusters <- readr::read_csv(file.path(dir_data, "file_names_clusters_ok.csv"),
                            show_col_types = FALSE) %>%
  dplyr::transmute(sample_id = toupper(as.character(name)),
                   cluster   = group_structure_CORE_admix)

design <- env_full %>%
  dplyr::inner_join(pc_coords, by = "sample_id") %>%
  dplyr::left_join(clusters, by = "sample_id") %>%
  dplyr::arrange(sample_id)

count_wide <- te_long %>%
  dplyr::select(sample_id, TE_ID, Count) %>%
  tidyr::pivot_wider(names_from = TE_ID, values_from = Count) %>%
  dplyr::arrange(sample_id)
count_mat <- as.matrix(count_wide[, -1, drop = FALSE])
rownames(count_mat) <- count_wide$sample_id

subset_ids <- lapply(SUBSET_CLUSTERS, function(cl) {
  sort(design$sample_id[design$cluster %in% cl])
})
for (nm in names(subset_ids)) {
  message(sprintf("Dataset %-9s n = %d", nm, length(subset_ids[[nm]])))
}

# TE universe: the same >= 10 observations filter as Model A, applied on the
# COMBINED dataset so all cells of the grid model the same TEs.
counts_combined <- count_mat[subset_ids$COMBINED, , drop = FALSE]
te_universe <- colnames(counts_combined)[colSums(!is.na(counts_combined)) >= MIN_OBS_PER_TE]
message(sprintf("TEs modelled: %d", length(te_universe)))

# =============================================================================
# 3. Model grid feasibility
# =============================================================================
# Each model uses 1 intercept + 1 bioclimatic term + depth PC terms. Cells
# left with fewer than MIN_RESIDUAL_DF residual degrees of freedom are not
# fitted: with n = 10, SYLE simply cannot support 5 or 10 PC covariates.
grid <- tidyr::expand_grid(subset = names(subset_ids), pc_depth = PC_DEPTHS) %>%
  dplyr::mutate(
    n_samples   = purrr::map_int(subset, ~ length(subset_ids[[.x]])),
    n_params    = 2L + as.integer(pc_depth),
    residual_df = n_samples - n_params,
    feasible    = residual_df >= MIN_RESIDUAL_DF
  )
readr::write_csv(grid, file.path(dir_tab, "model_grid_feasibility.csv"))
message("Model grid feasibility:")
print(as.data.frame(grid))

# =============================================================================
# 4. GLM fitting across the grid
# =============================================================================
#' Fit one quasi-Poisson GLM per TE for a single bioclimatic variable.
fit_bio_var <- function(bio_var, design_df, counts, te_ids, pc_names) {
  predictors <- c(bio_var, pc_names)
  model_df   <- design_df[, predictors, drop = FALSE]
  ok_rows    <- stats::complete.cases(model_df)
  model_df   <- model_df[ok_rows, , drop = FALSE]
  counts     <- counts[ok_rows, , drop = FALSE]

  rhs <- paste(sprintf("`%s`", predictors), collapse = " + ")
  fml <- stats::as.formula(paste("Count ~", rhs))

  n_te  <- length(te_ids)
  coefs <- rep(NA_real_, n_te); pvals <- rep(NA_real_, n_te)

  for (i in seq_len(n_te)) {
    dat <- model_df
    dat$Count <- as.numeric(counts[, te_ids[i]])
    fit <- suppressWarnings(tryCatch(
      stats::glm(fml, data = dat, family = stats::quasipoisson()),
      error = function(e) NULL))
    if (is.null(fit)) next
    cf <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
    # Row 2 is the bioclimatic term; the PC terms follow it.
    if (is.null(cf) || nrow(cf) < 2) next
    coefs[i] <- cf[2, 1]; pvals[i] <- cf[2, 4]
  }
  tibble::tibble(TE_ID = te_ids, variable = bio_var, coef = coefs, pval = pvals)
}

#' Fit all 19 bioclimatic variables for one subset x depth cell.
#' Two correction regimes are attached: BH within each bioclimatic variable
#' (padj_within, the Model A intent) and BH across all tests of the cell
#' (padj_global, stricter).
fit_cell <- function(subset_name, pc_depth) {
  ids       <- subset_ids[[subset_name]]
  pc_names  <- if (pc_depth > 0) sprintf("PC%d", seq_len(pc_depth)) else character(0)
  design_df <- design[match(ids, design$sample_id), , drop = FALSE]
  counts    <- count_mat[ids, te_universe, drop = FALSE]
  stopifnot(identical(rownames(counts), design_df$sample_id))

  out <- purrr::map_dfr(sprintf("bio%d", 1:19), function(bv) {
    fit_bio_var(bv, design_df, counts, te_universe, pc_names) %>%
      dplyr::filter(!is.na(pval)) %>%
      dplyr::mutate(padj_within = stats::p.adjust(pval, method = PADJ_METHOD))
  })
  out %>%
    dplyr::mutate(padj_global = stats::p.adjust(pval, method = PADJ_METHOD),
                  subset = subset_name, pc_depth = pc_depth)
}

grid_fit <- grid %>% dplyr::filter(feasible)
results <- list()

for (i in seq_len(nrow(grid_fit))) {
  sn <- grid_fit$subset[i]; pd <- grid_fit$pc_depth[i]
  key <- sprintf("%s_pc%d", sn, pd)
  message(sprintf("Fitting %s (n = %d, %d PC covariates) ...",
                  sn, grid_fit$n_samples[i], pd))
  t0 <- Sys.time()
  res <- fit_cell(sn, pd)
  message(sprintf("  %d TE x variable fits in %.1f min", nrow(res),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  readr::write_csv(res, file.path(dir_tab, sprintf("glm_%s.csv", key)))
  results[[key]] <- res
}

# =============================================================================
# 5. Candidate TE sets
# =============================================================================
#' Union, over the 19 bioclimatic variables, of the TEs significant in one
#' cell at one threshold under one correction regime.
candidate_set <- function(res, thr, regime = "padj_within") {
  sort(unique(res$TE_ID[!is.na(res[[regime]]) & res[[regime]] < thr]))
}

REGIMES <- c("padj_within", "padj_global")

sets <- list()
for (key in names(results)) {
  for (regime in REGIMES) {
    for (thr in THRESHOLDS) {
      sets[[paste(key, regime, thr, sep = "|")]] <-
        candidate_set(results[[key]], thr, regime)
    }
  }
}

set_sizes <- purrr::map_dfr(names(sets), function(k) {
  parts <- strsplit(k, "|", fixed = TRUE)[[1]]
  cell  <- parts[1]
  tibble::tibble(
    subset    = sub("_pc[0-9]+$", "", cell),
    pc_depth  = as.integer(sub("^.*_pc", "", cell)),
    regime    = parts[2],
    threshold = as.numeric(parts[3]),
    n_TE      = length(sets[[k]])
  )
}) %>%
  dplyr::arrange(regime, threshold, subset, pc_depth)

readr::write_csv(set_sizes, file.path(dir_tab, "candidate_set_sizes_by_grid.csv"))
message("Candidate set sizes (BH within variable):")
print(as.data.frame(set_sizes %>% dplyr::filter(regime == "padj_within")))

# Export the plain lists.
for (k in names(sets)) {
  parts <- strsplit(k, "|", fixed = TRUE)[[1]]
  thr_tag <- sub("\\.", "", format(as.numeric(parts[3]), scientific = FALSE))
  writeLines(sets[[k]],
             file.path(dir_list, sprintf("TE_%s_%s_padj%s.txt",
                                         parts[1], parts[2], thr_tag)))
}

# =============================================================================
# 6. What is shared and what is specific BETWEEN DATASETS, at each PC depth
# =============================================================================
jaccard <- function(a, b) {
  u <- length(union(a, b))
  if (u == 0) NA_real_ else length(intersect(a, b)) / u
}

overlap_subsets <- purrr::map_dfr(REGIMES, function(regime) {
  purrr::map_dfr(THRESHOLDS, function(thr) {
    purrr::map_dfr(PC_DEPTHS, function(pd) {
      avail <- names(subset_ids)[
        purrr::map_lgl(names(subset_ids), function(sn) {
          !is.null(sets[[paste(sprintf("%s_pc%d", sn, pd), regime, thr, sep = "|")]])
        })
      ]
      if (length(avail) == 0) return(NULL)
      s <- lapply(avail, function(sn)
        sets[[paste(sprintf("%s_pc%d", sn, pd), regime, thr, sep = "|")]])
      names(s) <- avail

      all_te <- unique(unlist(s))
      # How many of the datasets each candidate TE appears in.
      membership <- if (length(all_te) == 0) integer(0) else
        purrr::map_int(all_te, function(te) sum(purrr::map_lgl(s, ~ te %in% .x)))

      tibble::tibble(
        regime            = regime,
        threshold         = thr,
        pc_depth          = pd,
        datasets_fitted   = paste(avail, collapse = "+"),
        n_datasets        = length(avail),
        n_SYLW            = length(s[["SYLW"]] %||% character(0)),
        n_SYLE            = length(s[["SYLE"]] %||% character(0)),
        n_COMBINED        = length(s[["COMBINED"]] %||% character(0)),
        n_union           = length(all_te),
        n_in_all_datasets = sum(membership == length(avail)),
        n_in_exactly_two  = sum(membership == 2),
        n_dataset_specific= sum(membership == 1),
        jaccard_SYLW_COMB = jaccard(s[["SYLW"]] %||% character(0),
                                    s[["COMBINED"]] %||% character(0)),
        jaccard_SYLE_COMB = jaccard(s[["SYLE"]] %||% character(0),
                                    s[["COMBINED"]] %||% character(0)),
        jaccard_SYLW_SYLE = jaccard(s[["SYLW"]] %||% character(0),
                                    s[["SYLE"]] %||% character(0))
      )
    })
  })
})
readr::write_csv(overlap_subsets, file.path(dir_tab, "overlap_across_subsets.csv"))
message("Overlap across datasets (BH within variable):")
print(as.data.frame(overlap_subsets %>% dplyr::filter(regime == "padj_within")))

# The Model-A-style candidate list (significant in every fitted dataset) at
# each PC depth, which is the list the manuscript would actually report.
for (regime in REGIMES) {
  for (thr in THRESHOLDS) {
    for (pd in PC_DEPTHS) {
      keys <- purrr::keep(
        purrr::map_chr(names(subset_ids), ~ paste(sprintf("%s_pc%d", .x, pd),
                                                  regime, thr, sep = "|")),
        ~ !is.null(sets[[.x]])
      )
      if (length(keys) == 0) next
      inter <- Reduce(intersect, lapply(keys, function(k) sets[[k]]))
      thr_tag <- sub("\\.", "", format(thr, scientific = FALSE))
      writeLines(inter,
                 file.path(dir_list, sprintf("TE_INTERSECT_alldatasets_pc%d_%s_padj%s.txt",
                                             pd, regime, thr_tag)))
    }
  }
}

# =============================================================================
# 7. Effect of ADDING PC AXES, within each dataset
# =============================================================================
# Reference is depth 3, the number of axes used in 02_candidate_TEs.R.
overlap_depths <- purrr::map_dfr(REGIMES, function(regime) {
  purrr::map_dfr(THRESHOLDS, function(thr) {
    purrr::map_dfr(names(subset_ids), function(sn) {
      ref_key <- paste(sprintf("%s_pc3", sn), regime, thr, sep = "|")
      if (is.null(sets[[ref_key]])) return(NULL)
      ref <- sets[[ref_key]]
      purrr::map_dfr(PC_DEPTHS, function(pd) {
        k <- paste(sprintf("%s_pc%d", sn, pd), regime, thr, sep = "|")
        if (is.null(sets[[k]])) return(NULL)
        cur <- sets[[k]]
        tibble::tibble(
          regime = regime, threshold = thr, subset = sn, pc_depth = pd,
          n_TE = length(cur),
          n_shared_with_pc3 = length(intersect(cur, ref)),
          n_lost_vs_pc3     = length(setdiff(ref, cur)),
          n_gained_vs_pc3   = length(setdiff(cur, ref)),
          jaccard_vs_pc3    = jaccard(cur, ref)
        )
      })
    })
  })
})
readr::write_csv(overlap_depths, file.path(dir_tab, "overlap_across_depths.csv"))
message("Effect of PC depth, relative to PC1-3 (BH within variable):")
print(as.data.frame(overlap_depths %>% dplyr::filter(regime == "padj_within")))

# All pairwise Jaccard indices between every fitted cell, per regime/threshold.
jaccard_pairwise <- purrr::map_dfr(REGIMES, function(regime) {
  purrr::map_dfr(THRESHOLDS, function(thr) {
    cells <- names(results)
    combos <- utils::combn(cells, 2, simplify = FALSE)
    purrr::map_dfr(combos, function(cc) {
      a <- sets[[paste(cc[1], regime, thr, sep = "|")]]
      b <- sets[[paste(cc[2], regime, thr, sep = "|")]]
      if (is.null(a) || is.null(b)) return(NULL)
      tibble::tibble(regime = regime, threshold = thr,
                     cell_a = cc[1], cell_b = cc[2],
                     n_a = length(a), n_b = length(b),
                     n_shared = length(intersect(a, b)),
                     jaccard = jaccard(a, b))
    })
  })
})
readr::write_csv(jaccard_pairwise, file.path(dir_tab, "jaccard_pairwise.csv"))

# =============================================================================
# 8. Figures
# =============================================================================
# 8a. Scree / cumulative variance, with the depths actually tested marked.
p_scree <- variance_tbl %>%
  utils::head(20) %>%
  ggplot2::ggplot(ggplot2::aes(x = PC)) +
  ggplot2::geom_col(ggplot2::aes(y = variance_explained), fill = "#4C72B0",
                    width = 0.7) +
  ggplot2::geom_line(ggplot2::aes(y = cumulative_variance / 5), colour = "#C44E52",
                     linewidth = 0.8) +
  ggplot2::geom_point(ggplot2::aes(y = cumulative_variance / 5),
                      colour = "#C44E52", size = 1.5) +
  ggplot2::geom_vline(xintercept = setdiff(PC_DEPTHS, 0) + 0.5, linetype = 3,
                      colour = "grey40") +
  ggplot2::scale_y_continuous(
    name = "Variance explained per axis (%)",
    sec.axis = ggplot2::sec_axis(~ . * 5, name = "Cumulative variance (%)")
  ) +
  ggplot2::scale_x_continuous(breaks = 1:20) +
  ggplot2::labs(
    x = "SNP-based principal component",
    title = "How much genetic variance is gained by adding PC axes",
    subtitle = sprintf("PC1-3 = %.1f %%, PC1-5 = %.1f %%, PC1-10 = %.1f %%",
                       variance_tbl$cumulative_variance[3],
                       variance_tbl$cumulative_variance[5],
                       variance_tbl$cumulative_variance[10])
  ) +
  ggplot2::theme_classic(base_size = 10)
save_plot(p_scree, "pca_scree_cumulative_variance", width = 8, height = 5)

# 8b. Candidate set sizes across the grid.
p_grid <- set_sizes %>%
  dplyr::mutate(
    pc_label  = factor(sprintf("PC1-%d", pc_depth),
                       levels = sprintf("PC1-%d", PC_DEPTHS)),
    pc_label  = forcats::fct_recode(pc_label, `no PC` = "PC1-0"),
    threshold = factor(threshold, levels = THRESHOLDS,
                       labels = paste("padj <", THRESHOLDS)),
    regime    = factor(regime, levels = REGIMES,
                       labels = c("BH within variable", "BH across all tests"))
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = pc_label, y = n_TE, fill = subset)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
  ggplot2::geom_text(ggplot2::aes(label = n_TE),
                     position = ggplot2::position_dodge(width = 0.8),
                     vjust = -0.3, size = 2.2) +
  ggplot2::facet_grid(regime ~ threshold, scales = "free_y") +
  ggplot2::scale_fill_manual(values = c(SYLW = "#3171ad", SYLE = "#eee461",
                                        COMBINED = "#c66526")) +
  ggplot2::labs(x = "SNP-PC covariates included", y = "Candidate residual TEs",
                fill = NULL,
                title = "Candidate climate-associated TEs by dataset, PC depth and stringency") +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(legend.position = "bottom",
                 axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
save_plot(p_grid, "candidate_sets_across_grid", width = 11, height = 7)

# 8c. Shared vs dataset-specific candidates, per PC depth.
p_share <- overlap_subsets %>%
  dplyr::filter(regime == "padj_within") %>%
  dplyr::select(threshold, pc_depth, `In all datasets` = n_in_all_datasets,
                `In exactly two` = n_in_exactly_two,
                `Dataset-specific` = n_dataset_specific) %>%
  tidyr::pivot_longer(-c(threshold, pc_depth), names_to = "class",
                      values_to = "n") %>%
  dplyr::mutate(
    pc_label = factor(sprintf("PC1-%d", pc_depth),
                      levels = sprintf("PC1-%d", PC_DEPTHS)),
    pc_label = forcats::fct_recode(pc_label, `no PC` = "PC1-0"),
    threshold = factor(threshold, levels = THRESHOLDS,
                       labels = paste("padj <", THRESHOLDS)),
    class = factor(class, levels = c("In all datasets", "In exactly two",
                                     "Dataset-specific"))
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = pc_label, y = n, fill = class)) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::facet_wrap(~ threshold, scales = "free_y") +
  ggplot2::scale_fill_manual(values = c("In all datasets" = "#C44E52",
                                        "In exactly two" = "#DD8452",
                                        "Dataset-specific" = "#4C72B0")) +
  ggplot2::labs(x = "SNP-PC covariates included",
                y = "Candidate residual TEs", fill = NULL,
                title = "Shared versus dataset-specific candidates",
                subtitle = "BH within each bioclimatic variable") +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "bottom",
                 axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
save_plot(p_share, "shared_vs_specific_candidates", width = 9, height = 5)

# =============================================================================
# 9. Session log
# =============================================================================
log_session("01_snp_pca_axes")
message("=== 01_snp_pca_axes.R finished ===")
