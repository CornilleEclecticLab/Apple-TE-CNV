# =============================================================================
# Script:  10_focal_TE_pooled_map.R
# Purpose: Presentation-quality analysis and map for the strict candidate
#          TE_124385_LTR, estimated on ALL M. sylvestris individuals at once
#          instead of separately within the three population subsets.
#
#          WHY POOL NOW. The three-group split was a DISCOVERY device: a TE had
#          to respond to climate in each subset independently before it was
#          kept. It was never the right way to ESTIMATE the effect - it throws
#          away power by splitting n = 70 into 39 / 10 / 21. Once a candidate
#          is chosen, the pooled fit is the better estimator.
#
#          WHY THE PCs STAY. Pooling reintroduces exactly the population-
#          structure confounding that Reviewer 2 raised: without the subset
#          split, nothing separates a climate response from a north-west /
#          south-east genetic cline. The neutral SNP PCs are therefore kept as
#          covariates. At n = 70 with five parameters there is finally enough
#          power to carry them, which was the problem inside SYLE (n = 10).
#          Both versions are fitted and reported, so the cost of the covariates
#          is visible.
#
#          THE MAP. Rather than projecting a log-linear surface over latitude
#          and longitude, copy number is predicted from the fitted climate
#          model over the actual WorldClim temperature raster. The result
#          follows real climate structure - mountains, continentality, coastal
#          buffering - instead of a diagonal ramp. Cells whose temperature
#          falls outside the range actually sampled are masked, so the map
#          never shows an extrapolation dressed as a prediction.
#
# NAMING NOTE: "Model A" is the climate GLM (TE count ~ bioclimatic variable),
# "Model B" is the SPATIAL model (TE count ~ latitude + longitude, negative
# binomial) that yields the vector strength VS. The PC-covariate version of
# the climate model is "Model A-PC".
#
# Author:  Antho (analysis run with Claude Code)
# Date:    2026-08-24
# Project: TE_Project / review (reviewer response, Genome Biology submission)
#
# Inputs:
#   outputs/02_candidate_TEs/lists/
#       TE_STRICT_same_variable_same_sign_p005.txt
#   data/coord_syl_TE.csv, data/250923_sample_list_TE.csv
#   data/worldclim_env_data/climate/wc2.1_10m/wc2.1_10m_bio_*.tif
#   outputs/01_snp_pca_axes/tables/PCA_SNP_coordinates_PC1_10.csv
#
# Outputs (all under outputs/10_focal_TE_pooled_map/):
#   tables/pooled_vs_grouped_slopes.csv   pooled and per-group estimates
#   tables/pooled_model_summary.csv       slope, CI, % change per degree C
#   figures/<TE>_climate_map_<bio>.pdf|.svg  climate-driven prediction map
#   figures/<TE>_pooled_response.pdf|.svg    counts vs temperature, pooled fit
#   figures/<TE>_presentation.pdf|.svg        map + response assembled
#   ../../logs/10_focal_TE_pooled_map_*_sessionInfo.txt
# =============================================================================

# --- Setup -------------------------------------------------------------------
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(sf)
  library(rnaturalearth)
  library(terra)
  library(patchwork)
  library(svglite)
})

set.seed(20260824)

# --- Parameters --------------------------------------------------------------
CRS_EUROPE <- 3035
CRS_WGS84  <- 4326
MAP_LON <- c(-11, 33)      # presentation window
MAP_LAT <- c(35, 62)
MAP_RES_M <- 12000         # map raster resolution after reprojection (m)
PC_DEPTH <- 3
GROUPS <- c("SYLW", "SYLE", "SYL_Admix")
GROUP_COLOURS <- c(SYLW = "#2E6DA4", SYLE = "#B8A317", SYL_Admix = "#C4571F")
HEAT_COLOURS <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))

BIO_LABEL <- c(bio9  = "Mean temperature of the driest quarter (degC)",
               bio11 = "Mean temperature of the coldest quarter (degC)")

dir_out <- file.path(dir_outputs, "10_focal_TE_pooled_map")
dir_fig <- file.path(dir_out, "figures")
dir_tab <- file.path(dir_out, "tables")
for (d in c(dir_out, dir_fig, dir_tab)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

save_plot <- function(p, name, width = 7, height = 5) {
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".pdf")), p, width = width,
                  height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".svg")), p, width = width,
                  height = height, device = svglite::svglite)
}

message("=== 10_focal_TE_pooled_map.R ===")

# =============================================================================
# 1. Data
# =============================================================================
strict_te <- readLines(file.path(
  dir_outputs, "02_candidate_TEs", "lists",
  "TE_STRICT_same_variable_same_sign_p005.txt"))
TE <- strict_te[nzchar(strict_te)][1]
message("Focal TE: ", TE)

counts <- readr::read_csv(file.path(dir_data, "coord_syl_TE.csv"),
                          show_col_types = FALSE) %>%
  dplyr::mutate(Sample_ID = toupper(Sample_ID))

meta <- readr::read_csv(file.path(dir_data, "250923_sample_list_TE.csv"),
                        show_col_types = FALSE) %>%
  dplyr::filter(Species == "M. sylvestris", !is.na(latitude), !is.na(longitude),
                group_structure_CORE_admix_K6 != "Q2_Admix") %>%
  dplyr::transmute(
    Sample_ID = toupper(Sample_ID),
    latitude = as.numeric(latitude), longitude = as.numeric(longitude),
    group = dplyr::case_when(
      group_structure_CORE_admix_K6 %in% c("Q1_Admix", "Q5_Admix") ~ "SYL_Admix",
      group_structure_CORE_admix_K6 == "Q5" ~ "SYLE",
      group_structure_CORE_admix_K6 == "Q1" ~ "SYLW",
      TRUE ~ NA_character_)) %>%
  dplyr::filter(!is.na(group)) %>%
  dplyr::semi_join(counts, by = "Sample_ID")

clim_dir <- file.path(dir_data, "worldclim_env_data", "climate", "wc2.1_10m")
bio_stack <- terra::rast(file.path(clim_dir, sprintf("wc2.1_10m_bio_%d.tif", 1:19)))
names(bio_stack) <- sprintf("bio%d", 1:19)
meta <- dplyr::bind_cols(meta, tibble::as_tibble(terra::extract(
  bio_stack, as.matrix(meta[, c("longitude", "latitude")]))))

pc <- readr::read_csv(file.path(dir_outputs, "01_snp_pca_axes",
                                "tables", "PCA_SNP_coordinates_PC1_10.csv"),
                      show_col_types = FALSE) %>%
  dplyr::select(Sample_ID = sample_id, dplyr::all_of(sprintf("PC%d", 1:PC_DEPTH))) %>%
  dplyr::mutate(Sample_ID = toupper(Sample_ID))

dat <- counts %>%
  dplyr::select(Sample_ID, Count = dplyr::all_of(TE)) %>%
  dplyr::inner_join(meta, by = "Sample_ID") %>%
  dplyr::inner_join(pc, by = "Sample_ID") %>%
  dplyr::mutate(group = factor(group, levels = GROUPS))

message(sprintf("Individuals pooled: %d (previously split %s)",
                nrow(dat),
                paste(sprintf("%s=%d", levels(dat$group), table(dat$group)),
                      collapse = " / ")))

FOCAL_VARS <- intersect(c("bio9", "bio11"), names(dat))
pc_names <- sprintf("PC%d", seq_len(PC_DEPTH))

# =============================================================================
# 2. Pooled models, with and without the structure covariates
# =============================================================================
#' Extract the focal slope from a quasi-Poisson fit, with a Wald CI and the
#' multiplicative effect per unit (per degree C for a temperature variable).
slope_row <- function(fit, term, label, n) {
  cf <- summary(fit)$coefficients
  est <- cf[term, 1]; se <- cf[term, 2]
  tibble::tibble(
    model = label, n = n, term = term,
    slope = est, se = se,
    ci_low = est - 1.96 * se, ci_high = est + 1.96 * se,
    p_value = cf[term, 4],
    pct_change_per_unit = (exp(est) - 1) * 100,
    dispersion = summary(fit)$dispersion)
}

pooled <- purrr::map_dfr(FOCAL_VARS, function(bv) {
  f0 <- stats::glm(stats::as.formula(sprintf("Count ~ `%s`", bv)),
                   data = dat, family = stats::quasipoisson())
  f1 <- stats::glm(stats::as.formula(
    sprintf("Count ~ `%s` + %s", bv, paste(pc_names, collapse = " + "))),
    data = dat, family = stats::quasipoisson())
  dplyr::bind_rows(
    slope_row(f0, bv, "Pooled (climate only)", nrow(dat)),
    slope_row(f1, bv, "Pooled + PC1-3 (structure-controlled)", nrow(dat))
  ) %>% dplyr::mutate(variable = bv)
})

# Per-group estimates, for the comparison table.
per_group <- purrr::map_dfr(FOCAL_VARS, function(bv) {
  purrr::map_dfr(GROUPS, function(g) {
    d <- dat[dat$group == g, ]
    f <- stats::glm(stats::as.formula(sprintf("Count ~ `%s`", bv)),
                    data = d, family = stats::quasipoisson())
    slope_row(f, bv, paste0("Within ", g), nrow(d)) %>%
      dplyr::mutate(variable = bv)
  })
})

all_slopes <- dplyr::bind_rows(pooled, per_group) %>%
  dplyr::select(variable, model, n, slope, se, ci_low, ci_high, p_value,
                pct_change_per_unit, dispersion)
readr::write_csv(all_slopes,
                 file.path(dir_tab, "pooled_vs_grouped_slopes.csv"))
readr::write_csv(pooled, file.path(dir_tab, "pooled_model_summary.csv"))

message("Slopes (pooled first, then within each group):")
print(as.data.frame(all_slopes %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3)))))

# =============================================================================
# 3. Climate-driven prediction map
# =============================================================================
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  sf::st_transform(CRS_EUROPE) %>% sf::st_make_valid()

#' Predict copy number over the temperature raster, holding the genetic
#' covariates at their mean, and mask cells outside the sampled climate range.
build_climate_map <- function(bv) {
  fit <- stats::glm(stats::as.formula(
    sprintf("Count ~ `%s` + %s", bv, paste(pc_names, collapse = " + "))),
    data = dat, family = stats::quasipoisson())

  # The RASTER is reprojected, not the cell centres: transforming a lon/lat
  # grid point by point leaves an irregular mesh that geom_raster cannot draw.
  r <- terra::crop(bio_stack[[bv]],
                   terra::ext(MAP_LON[1], MAP_LON[2], MAP_LAT[1], MAP_LAT[2]))
  r <- terra::project(r, sprintf("EPSG:%d", CRS_EUROPE), method = "bilinear",
                      res = MAP_RES_M)
  df <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df) <- c("x", "y", "bio")

  obs_range <- range(dat[[bv]], na.rm = TRUE)
  df$in_range <- df$bio >= obs_range[1] & df$bio <= obs_range[2]

  nd <- data.frame(bio = df$bio)
  names(nd) <- bv
  # PCs held at the M. sylvestris mean: the map is the prediction for an
  # average genetic background, with climate the only thing varying.
  for (p in pc_names) nd[[p]] <- mean(dat[[p]])
  df$pred <- as.numeric(stats::predict(fit, newdata = nd, type = "response"))

  # Only in-range cells are drawn in colour; the rest stay neutral grey.
  df_in  <- df[df$in_range, ]
  df_out <- df[!df$in_range, ]

  pts <- sf::st_coordinates(sf::st_transform(
    sf::st_as_sf(dat, coords = c("longitude", "latitude"), crs = CRS_WGS84),
    CRS_EUROPE))
  d2 <- dat; d2$x <- pts[, 1]; d2$y <- pts[, 2]

  corners <- sf::st_coordinates(sf::st_transform(
    sf::st_as_sf(data.frame(lon = MAP_LON, lat = MAP_LAT),
                 coords = c("lon", "lat"), crs = CRS_WGS84), CRS_EUROPE))

  rng <- range(c(df_in$pred, d2$Count))
  sl <- pooled %>% dplyr::filter(variable == bv,
                                 grepl("PC1-3", model))

  p <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = df_out, ggplot2::aes(x = x, y = y),
                       fill = "grey86", width = MAP_RES_M, height = MAP_RES_M) +
    ggplot2::geom_tile(data = df_in, ggplot2::aes(x = x, y = y, fill = pred),
                       width = MAP_RES_M, height = MAP_RES_M) +
    ggplot2::geom_sf(data = world, fill = NA, colour = "grey25",
                     linewidth = 0.22) +
    ggplot2::geom_point(data = d2, ggplot2::aes(x = x, y = y, fill = Count),
                        shape = 21, size = 2.5, stroke = 0.5, colour = "black") +
    ggplot2::scale_fill_gradientn(colours = HEAT_COLOURS, limits = rng,
                                  name = "Copy number") +
    ggplot2::coord_sf(xlim = range(corners[, 1]), ylim = range(corners[, 2]),
                      expand = FALSE, crs = sf::st_crs(CRS_EUROPE)) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "#DCE6EC", colour = NA),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey30",
                                           linewidth = 0.4),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 8.5, colour = "grey25"),
      plot.caption = ggplot2::element_text(size = 7, colour = "grey40",
                                           hjust = 0)) +
    ggplot2::labs(
      title = sprintf("%s predicted copy number across Europe", TE),
      subtitle = sprintf(
        "Driven by %s | %.1f%% change per degC (95%% CI %.1f to %.1f), p = %.2g, n = %d",
        BIO_LABEL[[bv]], sl$pct_change_per_unit,
        (exp(sl$ci_low) - 1) * 100, (exp(sl$ci_high) - 1) * 100,
        sl$p_value, nrow(dat)),
      caption = paste0(
        "Quasi-Poisson GLM on all ", nrow(dat),
        " M. sylvestris individuals; neutral SNP PC1-3 held at their mean.\n",
        "Grey = outside the sampled temperature range (",
        sprintf("%.1f to %.1f degC", obs_range[1], obs_range[2]),
        "), not predicted. Points are observed counts on the same scale.\n",
        "Projection ETRS89 / LAEA Europe (EPSG:3035)."))

  list(plot = p, in_range = sum(df$in_range), total = nrow(df),
       obs_range = obs_range)
}

map_info <- list()
for (bv in FOCAL_VARS) {
  mi <- build_climate_map(bv)
  map_info[[bv]] <- mi
  save_plot(mi$plot, sprintf("%s_climate_map_%s", TE, bv),
            width = 8.5, height = 6.4)
  message(sprintf("  %s: %d of %d cells within the sampled range (%.0f%%)",
                  bv, mi$in_range, mi$total, 100 * mi$in_range / mi$total))
}

# =============================================================================
# 4. Pooled response panel
# =============================================================================
resp <- purrr::map_dfr(FOCAL_VARS, function(bv) {
  dat %>% dplyr::transmute(group, Count, value = .data[[bv]],
                           variable = BIO_LABEL[[bv]])
})
resp$variable <- factor(resp$variable, levels = unname(BIO_LABEL[FOCAL_VARS]))

p_resp <- ggplot2::ggplot(resp, ggplot2::aes(x = value, y = Count)) +
  ggplot2::geom_smooth(method = "glm", se = TRUE,
                       method.args = list(family = stats::quasipoisson()),
                       colour = "black", fill = "grey70", linewidth = 0.9,
                       alpha = 0.3) +
  ggplot2::geom_point(ggplot2::aes(fill = group), shape = 21, size = 2.4,
                      colour = "grey20", alpha = 0.9) +
  ggplot2::facet_wrap(~ variable, scales = "free_x", ncol = 1) +
  ggplot2::scale_fill_manual(values = GROUP_COLOURS, name = NULL) +
  ggplot2::labs(x = NULL, y = sprintf("%s copy number", TE),
                title = "Pooled climate response",
                subtitle = sprintf("Single quasi-Poisson fit on all %d individuals; colours show the discovery subsets only",
                                   nrow(dat))) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "bottom",
                 strip.background = ggplot2::element_rect(fill = "grey95",
                                                          colour = "grey70"),
                 strip.text = ggplot2::element_text(size = 8.5),
                 plot.title = ggplot2::element_text(face = "bold", size = 11),
                 plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30"),
                 plot.background = ggplot2::element_rect(fill = "white",
                                                         colour = NA))
save_plot(p_resp, sprintf("%s_pooled_response", TE), width = 6, height = 7)

# =============================================================================
# 5. Presentation composite
# =============================================================================
main_var <- if ("bio11" %in% FOCAL_VARS) "bio11" else FOCAL_VARS[1]
p_pres <- (map_info[[main_var]]$plot | p_resp) +
  patchwork::plot_layout(widths = c(1.45, 1))
save_plot(p_pres, sprintf("%s_presentation", TE), width = 14.5, height = 6.8)

# =============================================================================
# 6. Session log
# =============================================================================
log_session("10_focal_TE_pooled_map")
message("=== 10_focal_TE_pooled_map.R finished ===")
