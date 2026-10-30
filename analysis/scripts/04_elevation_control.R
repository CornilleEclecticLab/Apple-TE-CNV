# =============================================================================
# Script:  04_elevation_control.R
# Purpose: Ask whether residual TE copy number tracks ALTITUDE in
#          M. sylvestris, using exactly the design that produced the climate
#          candidates, so the two are directly comparable.
#
#          Same three-group replication design as the climate analysis:
#          SYLW (Q1, n = 39), SYLE (Q5, n = 10) and SYL_Admix
#          (Q1_Admix + Q5_Admix, n = 21). Per group, a quasi-Poisson GLM of TE
#          count on elevation; TEs reaching the threshold are unioned within
#          each group and the three group lists are intersected. The
#          PC-covariate version (Model A-PC style) is fitted alongside.
#
#          THE QUESTION THAT MATTERS. Elevation and climate are not
#          independent: temperature falls with altitude, so a TE that looks
#          "altitude-associated" may simply be the temperature association
#          already reported. Two things are therefore done rather than
#          assumed away:
#            - the correlation between elevation and each of the 19
#              bioclimatic variables is computed, overall and within group;
#            - every elevation candidate is refitted as
#                  Count ~ elevation + bio_i
#              for all 19 bioclimatic variables in all three groups, and only
#              those keeping a significant elevation term against EVERY
#              climate variable are reported as independent of climate.
#
#          Elevation is WorldClim 2.1 at 30 arc-seconds (~1 km), a finer grid
#          than the 10 arc-minute bioclim layers used for the climate models,
#          because elevation varies far more sharply in space than climate
#          does. The value at 10 arc-minutes is also extracted so the two
#          resolutions can be compared.
#
# NAMING NOTE: in the manuscript, "Model A" is the climate GLM
# (TE count ~ bioclimatic variable) and "Model B" is the SPATIAL model
# (TE count ~ latitude + longitude, negative binomial), which produces the
# vector-strength (VS) values shown in Fig. 3A-B. To avoid a clash with that
# established terminology, the PC-covariate version of a model is called
# "-PC" throughout this folder, never "Model B".
#
# Author:  Antho (analysis run with Claude Code)
# Date:    2026-08-24
# Project: TE_Project / review (reviewer response, Genome Biology submission)
#
# Inputs:
#   data/elevation_data/wc2.1_30s_elev.tif   WorldClim 2.1 elevation, 30 arc-sec
#   data/worldclim_env_data/climate/wc2.1_10m/wc2.1_10m_bio_*.tif
#   data/coord_syl_TE.csv                    residual TE counts
#   data/250923_sample_list_TE.csv           coordinates and K6 clusters
#   outputs/01_snp_pca_axes/tables/PCA_SNP_coordinates_PC1_10.csv
#   outputs/02_candidate_TEs/tables/
#       FINAL_union_candidates.csv           the 38 climate candidates
#
# Outputs (all under outputs/04_elevation_control/):
#   tables/site_elevation.csv                elevation per accession
#   tables/elevation_vs_bioclim_correlation.csv
#   tables/glm_elevation_by_group.csv        full GLM results
#   tables/candidate_counts.csv              per model / threshold
#   tables/elevation_candidates_annotated.csv
#   tables/elevation_vs_climate_overlap.csv  against the 38 climate candidates
#   tables/conditional_on_climate.csv        elevation term given each bio_i
#   lists/*.txt
#   figures/*.pdf + *.svg
#   ../../logs/04_elevation_control_*_sessionInfo.txt
# =============================================================================

# --- Setup -------------------------------------------------------------------
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(terra)
  library(svglite)
})

set.seed(20260824)

# --- Parameters --------------------------------------------------------------
THRESHOLDS  <- c(0.05, 0.01, 0.001)
PADJ_METHOD <- "BH"
PC_DEPTH    <- 3
GROUPS      <- c("SYLW", "SYLE", "SYL_Admix")
MIN_N       <- 3

GROUP_COLOURS <- c(SYLW = "#2E6DA4", SYLE = "#B8A317", SYL_Admix = "#C4571F")

dir_out  <- file.path(dir_outputs, "04_elevation_control")
dir_tab  <- file.path(dir_out, "tables")
dir_fig  <- file.path(dir_out, "figures")
dir_list <- file.path(dir_out, "lists")
for (d in c(dir_out, dir_tab, dir_fig, dir_list)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

save_plot <- function(p, name, width = 7, height = 5) {
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".pdf")), p, width = width,
                  height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".svg")), p, width = width,
                  height = height, device = svglite::svglite)
}

message("=== 04_elevation_control.R ===")

# =============================================================================
# 1. Data, and elevation extraction
# =============================================================================
counts <- readr::read_csv(file.path(dir_data, "coord_syl_TE.csv"),
                          show_col_types = FALSE) %>%
  dplyr::mutate(Sample_ID = toupper(Sample_ID))
te_vars <- setdiff(names(counts), "Sample_ID")

meta <- readr::read_csv(file.path(dir_data, "250923_sample_list_TE.csv"),
                        show_col_types = FALSE) %>%
  dplyr::filter(Species == "M. sylvestris",
                !is.na(latitude), !is.na(longitude),
                group_structure_CORE_admix_K6 != "Q2_Admix") %>%
  dplyr::transmute(
    Sample_ID = toupper(Sample_ID),
    latitude  = as.numeric(latitude),
    longitude = as.numeric(longitude),
    group = dplyr::case_when(
      group_structure_CORE_admix_K6 %in% c("Q1_Admix", "Q5_Admix") ~ "SYL_Admix",
      group_structure_CORE_admix_K6 == "Q5" ~ "SYLE",
      group_structure_CORE_admix_K6 == "Q1" ~ "SYLW",
      TRUE ~ NA_character_)) %>%
  dplyr::filter(!is.na(group)) %>%
  dplyr::semi_join(counts, by = "Sample_ID")

elev_file <- list.files(file.path(dir_data, "elevation_data"),
                        pattern = "elev.*\\.tif$", full.names = TRUE)[1]
stopifnot(!is.na(elev_file), file.exists(elev_file))
elev_r <- terra::rast(elev_file)

xy <- as.matrix(meta[, c("longitude", "latitude")])
meta$elevation <- as.numeric(terra::extract(elev_r, xy)[, 1])

# Bioclim at the same localities, for the confounding checks below.
clim_dir <- file.path(dir_data, "worldclim_env_data", "climate", "wc2.1_10m")
bio_stack <- terra::rast(file.path(clim_dir, sprintf("wc2.1_10m_bio_%d.tif", 1:19)))
names(bio_stack) <- sprintf("bio%d", 1:19)
meta <- dplyr::bind_cols(meta, tibble::as_tibble(terra::extract(bio_stack, xy)))

pc <- readr::read_csv(
  file.path(dir_outputs, "01_snp_pca_axes", "tables",
            "PCA_SNP_coordinates_PC1_10.csv"), show_col_types = FALSE) %>%
  dplyr::select(Sample_ID = sample_id, dplyr::all_of(sprintf("PC%d", 1:PC_DEPTH))) %>%
  dplyr::mutate(Sample_ID = toupper(Sample_ID))

design <- meta %>% dplyr::inner_join(pc, by = "Sample_ID") %>%
  dplyr::arrange(Sample_ID)
count_mat <- as.matrix(counts[match(design$Sample_ID, counts$Sample_ID),
                              te_vars, drop = FALSE])
rownames(count_mat) <- design$Sample_ID

stopifnot(sum(is.na(design$elevation)) == 0)
readr::write_csv(design %>% dplyr::select(Sample_ID, group, longitude, latitude,
                                          elevation),
                 file.path(dir_tab, "site_elevation.csv"))

message("Elevation by group (m):")
print(as.data.frame(design %>% dplyr::group_by(group) %>%
  dplyr::summarise(n = dplyr::n(), min = min(elevation),
                   median = stats::median(elevation),
                   mean = round(mean(elevation)), max = max(elevation),
                   .groups = "drop")))

# =============================================================================
# 2. How confounded is elevation with climate?
# =============================================================================
bio_names <- sprintf("bio%d", 1:19)
cor_tbl <- purrr::map_dfr(bio_names, function(bv) {
  overall <- stats::cor(design$elevation, design[[bv]], use = "complete.obs")
  per_grp <- purrr::map_dbl(GROUPS, function(g) {
    d <- design[design$group == g, ]
    if (nrow(d) < 4) return(NA_real_)
    stats::cor(d$elevation, d[[bv]], use = "complete.obs")
  })
  tibble::tibble(variable = bv, r_overall = overall,
                 r_SYLW = per_grp[1], r_SYLE = per_grp[2],
                 r_SYL_Admix = per_grp[3])
}) %>%
  dplyr::arrange(dplyr::desc(abs(r_overall)))
readr::write_csv(cor_tbl,
                 file.path(dir_tab, "elevation_vs_bioclim_correlation.csv"))
message("Elevation vs bioclim, strongest correlations overall:")
print(as.data.frame(utils::head(cor_tbl %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 3))), 6)))

# =============================================================================
# 3. GLMs: elevation, with and without the SNP-PC covariates
# =============================================================================
pc_names <- sprintf("PC%d", seq_len(PC_DEPTH))

#' Quasi-Poisson GLM of every TE on a set of predictors, within one group.
#' Returns the coefficient and p-value of the FIRST predictor.
fit_group <- function(grp, predictors) {
  idx <- which(design$group == grp)
  d   <- design[idx, , drop = FALSE]
  cm  <- count_mat[idx, , drop = FALSE]

  md <- d[, predictors, drop = FALSE]
  ok <- stats::complete.cases(md)
  md <- md[ok, , drop = FALSE]; cm <- cm[ok, , drop = FALSE]

  rhs <- paste(sprintf("`%s`", predictors), collapse = " + ")
  fml <- stats::as.formula(paste("te_count ~", rhs))

  n_te <- length(te_vars)
  est <- rep(NA_real_, n_te); pv <- rep(NA_real_, n_te)
  for (i in seq_len(n_te)) {
    y <- as.numeric(cm[, te_vars[i]])
    if (length(unique(y)) < 2 || length(y) < MIN_N) next
    dat <- md; dat$te_count <- y
    fit <- suppressWarnings(tryCatch(
      stats::glm(fml, data = dat, family = stats::quasipoisson()),
      error = function(e) NULL))
    if (is.null(fit)) next
    cf <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
    if (is.null(cf) || nrow(cf) < 2) next
    est[i] <- cf[2, 1]; pv[i] <- cf[2, 4]
  }
  tibble::tibble(group = grp, TE_ID = te_vars, estimate = est, p_value = pv)
}

run_model <- function(predictors, label) {
  message(sprintf("Fitting %s ...", label))
  t0 <- Sys.time()
  out <- purrr::map_dfr(GROUPS, fit_group, predictors = predictors) %>%
    dplyr::filter(!is.na(p_value)) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(padj = stats::p.adjust(p_value, method = PADJ_METHOD)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(model = label)
  message(sprintf("  %d fits in %.1f min", nrow(out),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  out
}

elev_A  <- run_model("elevation", "elevation only")
elev_PC <- run_model(c("elevation", pc_names), "elevation + PC1-3")
glm_all <- dplyr::bind_rows(elev_A, elev_PC)
readr::write_csv(glm_all, file.path(dir_tab, "glm_elevation_by_group.csv"))

# =============================================================================
# 4. Candidates: replicated across all three groups
# =============================================================================
replicated <- function(df, thr, p_col = "p_value") {
  s <- lapply(GROUPS, function(g)
    unique(df$TE_ID[df$group == g & !is.na(df[[p_col]]) & df[[p_col]] < thr]))
  names(s) <- GROUPS
  list(per_group = s, shared = sort(Reduce(intersect, s)))
}

counts_tbl <- purrr::map_dfr(unique(glm_all$model), function(m) {
  df <- glm_all[glm_all$model == m, ]
  purrr::map_dfr(c("p_value", "padj"), function(pc_col) {
    purrr::map_dfr(THRESHOLDS, function(thr) {
      r <- replicated(df, thr, pc_col)
      tibble::tibble(model = m,
                     regime = ifelse(pc_col == "p_value", "raw p", "BH within group"),
                     threshold = thr,
                     n_SYLW = length(r$per_group$SYLW),
                     n_SYLE = length(r$per_group$SYLE),
                     n_SYL_Admix = length(r$per_group$SYL_Admix),
                     n_replicated = length(r$shared))
    })
  })
})
readr::write_csv(counts_tbl, file.path(dir_tab, "candidate_counts.csv"))
message("Elevation candidates (n_replicated = shared by all three groups):")
print(as.data.frame(counts_tbl))

for (m in unique(glm_all$model)) {
  for (pc_col in c("p_value", "padj")) {
    for (thr in THRESHOLDS) {
      s <- replicated(glm_all[glm_all$model == m, ], thr, pc_col)$shared
      tag <- sprintf("%s_%s_p%s",
                     ifelse(grepl("PC", m), "elevPC", "elev"),
                     ifelse(pc_col == "p_value", "rawp", "BH"),
                     sub("\\.", "", format(thr, scientific = FALSE)))
      writeLines(s, file.path(dir_list, sprintf("TE_%s.txt", tag)))
    }
  }
}

# The headline set: same threshold and rule as the published climate list.
MAIN_THR <- 0.01
elev_cand <- replicated(elev_A, MAIN_THR, "p_value")$shared
elev_cand_pc <- replicated(elev_PC, MAIN_THR, "p_value")$shared
message(sprintf("Elevation candidates at raw p < %g: %d (elevation only), %d (+ PC1-3)",
                MAIN_THR, length(elev_cand), length(elev_cand_pc)))

# =============================================================================
# 5. Overlap with the climate candidates
# =============================================================================
climate38 <- readr::read_csv(
  file.path(dir_outputs, "02_candidate_TEs", "tables",
            "FINAL_union_candidates.csv"), show_col_types = FALSE)
clim_ids <- sort(unique(climate38$TE_ID))

overlap <- tibble::tibble(
  set = c("elevation only", "elevation + PC1-3"),
  n = c(length(elev_cand), length(elev_cand_pc)),
  n_shared_with_climate38 = c(length(intersect(elev_cand, clim_ids)),
                              length(intersect(elev_cand_pc, clim_ids))),
  n_specific_to_elevation = c(length(setdiff(elev_cand, clim_ids)),
                              length(setdiff(elev_cand_pc, clim_ids)))
)
readr::write_csv(overlap, file.path(dir_tab, "elevation_vs_climate_overlap.csv"))
message("Overlap with the 38 climate candidates:")
print(as.data.frame(overlap))

superfamily_of <- function(x) sub("^TE_[0-9]+_", "", x)
cand_tbl <- tibble::tibble(
  TE_ID = elev_cand,
  superfamily = superfamily_of(elev_cand),
  in_climate38 = elev_cand %in% clim_ids,
  in_elevation_PC = elev_cand %in% elev_cand_pc
) %>%
  dplyr::left_join(
    elev_A %>% dplyr::filter(TE_ID %in% elev_cand, p_value < MAIN_THR) %>%
      dplyr::group_by(TE_ID) %>%
      dplyr::summarise(min_p = min(p_value),
                       mean_slope = mean(estimate), .groups = "drop"),
    by = "TE_ID") %>%
  dplyr::arrange(min_p)
readr::write_csv(cand_tbl,
                 file.path(dir_tab, "elevation_candidates_annotated.csv"))

# =============================================================================
# 6. Is it altitude, or is it climate wearing altitude's clothes?
# =============================================================================
# Each elevation candidate is refitted with elevation AND one bioclimatic
# variable, for all 19 variables in all three groups. A TE only counts as
# independent of climate if its elevation term stays significant in all three
# groups against every single climate variable.
if (length(elev_cand) > 0) {
  message(sprintf("Testing %d elevation candidates against each of the 19 bioclimatic variables ...",
                  length(elev_cand)))

  cond <- purrr::map_dfr(GROUPS, function(g) {
    idx <- which(design$group == g)
    d <- design[idx, , drop = FALSE]
    cm <- count_mat[idx, elev_cand, drop = FALSE]
    purrr::map_dfr(bio_names, function(bv) {
      purrr::map_dfr(seq_along(elev_cand), function(i) {
        dat <- d[, c("elevation", bv)]
        dat$te_count <- as.numeric(cm[, i])
        fml <- stats::as.formula(sprintf("te_count ~ elevation + `%s`", bv))
        fit <- suppressWarnings(tryCatch(
          stats::glm(fml, data = dat, family = stats::quasipoisson()),
          error = function(e) NULL))
        p <- NA_real_
        if (!is.null(fit)) {
          cf <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
          if (!is.null(cf) && "elevation" %in% rownames(cf)) p <- cf["elevation", 4]
        }
        tibble::tibble(group = g, variable = bv, TE_ID = elev_cand[i],
                       p_elevation = p)
      })
    })
  })
  readr::write_csv(cond, file.path(dir_tab, "conditional_on_climate.csv"))

  robust <- cond %>%
    dplyr::group_by(TE_ID) %>%
    dplyr::summarise(
      n_tests = dplyr::n(),
      n_sig   = sum(!is.na(p_elevation) & p_elevation < 0.05),
      worst_p = suppressWarnings(max(p_elevation, na.rm = TRUE)),
      robust_to_all_climate = n_sig == n_tests,
      .groups = "drop") %>%
    dplyr::arrange(worst_p)
  readr::write_csv(robust,
                   file.path(dir_tab, "elevation_robust_to_climate.csv"))
  message(sprintf("Elevation candidates keeping p < 0.05 against ALL 19 climate variables in ALL three groups: %d / %d",
                  sum(robust$robust_to_all_climate), nrow(robust)))
  print(as.data.frame(utils::head(robust, 10)))
}

# =============================================================================
# 7. Figures
# =============================================================================
p_elev <- design %>%
  ggplot2::ggplot(ggplot2::aes(x = group, y = elevation, fill = group)) +
  ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75,
                        show.legend = FALSE) +
  ggplot2::geom_jitter(width = 0.12, height = 0, size = 1.8, alpha = 0.85,
                       shape = 21, colour = "grey20", show.legend = FALSE) +
  ggplot2::scale_fill_manual(values = GROUP_COLOURS) +
  ggplot2::labs(x = NULL, y = "Elevation (m)",
                title = "Elevation of the M. sylvestris sampling localities",
                subtitle = "WorldClim 2.1, 30 arc-seconds") +
  ggplot2::theme_classic(base_size = 11)
save_plot(p_elev, "elevation_by_group", width = 6, height = 5)

p_cor <- cor_tbl %>%
  dplyr::select(variable, Overall = r_overall, SYLW = r_SYLW, SYLE = r_SYLE,
                SYL_Admix = r_SYL_Admix) %>%
  tidyr::pivot_longer(-variable, names_to = "set", values_to = "r") %>%
  dplyr::mutate(variable = factor(variable, levels = rev(bio_names)),
                set = factor(set, levels = c("Overall", "SYLW", "SYLE",
                                             "SYL_Admix"))) %>%
  ggplot2::ggplot(ggplot2::aes(x = set, y = variable, fill = r)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", r)), size = 2.6) +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                               midpoint = 0, limits = c(-1, 1),
                               name = "Pearson r") +
  ggplot2::labs(x = NULL, y = NULL,
                title = "Elevation is strongly correlated with several bioclimatic variables",
                subtitle = "An 'altitude' association can therefore be a temperature association in disguise") +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(panel.grid = ggplot2::element_blank())
save_plot(p_cor, "elevation_vs_bioclim_correlation", width = 6.5, height = 7.5)

p_counts <- counts_tbl %>%
  dplyr::mutate(threshold = factor(threshold, levels = THRESHOLDS,
                                   labels = paste("p <", THRESHOLDS))) %>%
  ggplot2::ggplot(ggplot2::aes(x = threshold, y = n_replicated, fill = model)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
  ggplot2::geom_text(ggplot2::aes(label = n_replicated),
                     position = ggplot2::position_dodge(width = 0.8),
                     vjust = -0.3, size = 3) +
  ggplot2::facet_wrap(~ regime, scales = "free_y") +
  ggplot2::scale_fill_manual(values = c("elevation only" = "#4C72B0",
                                        "elevation + PC1-3" = "#C44E52")) +
  ggplot2::labs(x = NULL, y = "TEs replicated across all three groups",
                fill = NULL,
                title = "Residual TEs associated with elevation",
                subtitle = "Same three-group replication design as the climate analysis") +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")
save_plot(p_counts, "elevation_candidate_counts", width = 9, height = 5)

# =============================================================================
# 8. Session log
# =============================================================================
log_session("04_elevation_control")
message("=== 04_elevation_control.R finished ===")
