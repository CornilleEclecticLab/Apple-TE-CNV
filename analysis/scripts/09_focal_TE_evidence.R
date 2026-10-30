# =============================================================================
# Script:  09_focal_TE_evidence.R
# Purpose: A focused figure for the candidate TE that survives the STRICTEST
#          reading of the replication rule - significant for the SAME
#          bioclimatic variable, with the SAME slope sign, in all three
#          M. sylvestris groups.
#
#          The consistency audit in 02_candidate_TEs.R
#          showed that of the 34 candidates behind Figure 3, only one meets
#          that bar. It is TE_124385_LTR, and it does so for two temperature
#          variables at once:
#              bio9  - mean temperature of the driest quarter
#              bio11 - mean temperature of the coldest quarter
#          with negative slopes in all three groups, i.e. more copies where it
#          is colder.
#
#          The script does not hard-code that TE: it reads the strict list
#          produced by script 05, so if the threshold or the design changes,
#          the figure follows.
#
#          HONEST CAVEAT, printed in the outputs. Under the PC-covariate model
#          this TE stays significant only in SYLW (bio11, p = 0.0095); SYLE and
#          SYL_Admix lose significance. The permutation control in script 05
#          showed that the PC-covariate models are heavily power-limited in
#          these small groups, so this is expected rather than a contradiction,
#          but it must be stated if the TE is presented on its own.
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
#   outputs/02_candidate_TEs/tables/glm_modelA_PC_by_group.csv
#   outputs/03_spatial_vs/tables/VS_table_38_candidates.csv
#   data/modelA_fig3_raw_p.csv, data/coord_syl_TE.csv,
#   data/250923_sample_list_TE.csv, data/worldclim_env_data/...
#
# Outputs (all under outputs/09_focal_TE_evidence/):
#   tables/strict_candidate_evidence.csv   every slope and p-value, both models
#   figures/<TE>_prediction_map.pdf|.svg   modelled copy-number surface
#   figures/<TE>_climate_response.pdf|.svg count vs bio9 and bio11, by group
#   figures/<TE>_composite.pdf|.svg        map + climate response assembled
#   ../../logs/09_focal_TE_evidence_*_sessionInfo.txt
# =============================================================================

# --- Setup -------------------------------------------------------------------
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(sf)
  library(rnaturalearth)
  library(glmmTMB)
  library(terra)
  library(patchwork)
  library(svglite)
})

set.seed(20260824)

# --- Parameters --------------------------------------------------------------
CRS_EUROPE <- 3035
CRS_WGS84  <- 4326
GRID_KM    <- 11
HULL_BUFFER_M <- 200000
ARROW_MAX_DEG <- 5.2
GROUPS <- c("SYLW", "SYLE", "SYL_Admix")
GROUP_COLOURS <- c(SYLW = "#2E6DA4", SYLE = "#B8A317", SYL_Admix = "#C4571F")
HEAT_COLOURS <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))

# WorldClim bioclimatic variable definitions, for axis labels.
BIO_LABEL <- c(
  bio1 = "Annual mean temperature (degC)",
  bio5 = "Max temperature of warmest month (degC)",
  bio6 = "Min temperature of coldest month (degC)",
  bio8 = "Mean temperature of wettest quarter (degC)",
  bio9 = "Mean temperature of driest quarter (degC)",
  bio10 = "Mean temperature of warmest quarter (degC)",
  bio11 = "Mean temperature of coldest quarter (degC)",
  bio12 = "Annual precipitation (mm)",
  bio13 = "Precipitation of wettest month (mm)",
  bio14 = "Precipitation of driest month (mm)",
  bio16 = "Precipitation of wettest quarter (mm)",
  bio17 = "Precipitation of driest quarter (mm)",
  bio18 = "Precipitation of warmest quarter (mm)",
  bio19 = "Precipitation of coldest quarter (mm)")

dir_out <- file.path(dir_outputs, "09_focal_TE_evidence")
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

message("=== 09_focal_TE_evidence.R ===")

# =============================================================================
# 1. Which TE, and for which variables?
# =============================================================================
strict_file <- file.path(dir_outputs, "02_candidate_TEs",
                         "lists", "TE_STRICT_same_variable_same_sign_p005.txt")
strict_te <- readLines(strict_file)
strict_te <- strict_te[nzchar(strict_te)]
if (length(strict_te) == 0) {
  stop("No TE passes the strict rule; nothing to draw.")
}
message(sprintf("Strict candidate(s): %s", paste(strict_te, collapse = ", ")))

modelA <- readr::read_csv(file.path(dir_data, "modelA_fig3_raw_p.csv"),
                          show_col_types = FALSE) %>%
  dplyr::mutate(variable = sub("^wc2[.]1_10m_bio_", "bio", bio_var))
modelA_PC <- readr::read_csv(
  file.path(dir_outputs, "02_candidate_TEs", "tables",
            "glm_modelA_PC_by_group.csv"), show_col_types = FALSE)

#' Bioclimatic variables significant in all three groups with a consistent
#' slope sign, for one TE.
strict_variables <- function(te, thr = 0.05) {
  d <- modelA %>% dplyr::filter(TE_ID == te, !is.na(p_value), p_value < thr)
  shared <- Reduce(intersect, lapply(GROUPS, function(g) d$variable[d$group == g]))
  shared[vapply(shared, function(v) {
    sg <- sign(d$estimate[d$variable == v]); length(unique(sg)) == 1
  }, logical(1))]
}

TE <- strict_te[1]
strict_vars <- strict_variables(TE)
message(sprintf("%s is consistently associated with: %s", TE,
                paste(strict_vars, collapse = ", ")))

evidence <- dplyr::bind_rows(
  modelA %>% dplyr::filter(TE_ID == TE, variable %in% strict_vars) %>%
    dplyr::transmute(TE_ID, model = "Model A (bio only)", group, variable,
                     estimate, p_value),
  modelA_PC %>% dplyr::filter(TE_ID == TE, variable %in% strict_vars) %>%
    dplyr::transmute(TE_ID, model = "Model A-PC (bio + PC1-3)", group, variable,
                     estimate, p_value)
) %>%
  dplyr::arrange(variable, model, group)
readr::write_csv(evidence, file.path(dir_tab, "strict_candidate_evidence.csv"))
message("Evidence for the strict candidate:")
print(as.data.frame(evidence %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3)))))

# =============================================================================
# 2. Data
# =============================================================================
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
meta <- dplyr::bind_cols(
  meta, tibble::as_tibble(terra::extract(
    bio_stack, as.matrix(meta[, c("longitude", "latitude")]))))

dat <- counts %>%
  dplyr::select(Sample_ID, Count = dplyr::all_of(TE)) %>%
  dplyr::inner_join(meta, by = "Sample_ID") %>%
  dplyr::mutate(group = factor(group, levels = GROUPS))

# =============================================================================
# 3. Panel A - modelled copy-number surface (Model B, the spatial model)
# =============================================================================
vs_tbl <- readr::read_csv(
  file.path(dir_outputs, "03_spatial_vs", "tables",
            "VS_table_38_candidates.csv"), show_col_types = FALSE)
row <- vs_tbl[vs_tbl$TE_ID == TE, ]

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  sf::st_transform(CRS_EUROPE) %>% sf::st_make_valid()
sites_sf <- sf::st_as_sf(dat, coords = c("longitude", "latitude"),
                         crs = CRS_WGS84, remove = FALSE) %>%
  sf::st_transform(CRS_EUROPE)
hull <- sf::st_convex_hull(sf::st_union(sites_sf)) %>%
  sf::st_buffer(HULL_BUFFER_M)
world_crop <- suppressWarnings(
  sf::st_crop(world, sf::st_bbox(sf::st_buffer(hull, 250000))))
region <- suppressWarnings(
  sf::st_intersection(hull, sf::st_make_valid(sf::st_union(world_crop))))

bb <- sf::st_bbox(hull)
grid <- expand.grid(x = seq(bb["xmin"], bb["xmax"], by = GRID_KM * 1000),
                    y = seq(bb["ymin"], bb["ymax"], by = GRID_KM * 1000))
grid_sf <- sf::st_as_sf(grid, coords = c("x", "y"), crs = CRS_EUROPE)
keep <- lengths(sf::st_intersects(grid_sf, region)) > 0
grid <- grid[keep, ]; grid_sf <- grid_sf[keep, ]
ll <- sf::st_coordinates(sf::st_transform(grid_sf, CRS_WGS84))
grid$longitude <- ll[, 1]; grid$latitude <- ll[, 2]

fit_sp <- glmmTMB::glmmTMB(Count ~ latitude + longitude, data = dat,
                           family = glmmTMB::nbinom2())
grid$pred <- as.numeric(stats::predict(fit_sp, newdata = grid, type = "response"))

dat$x <- sf::st_coordinates(sites_sf)[, 1]
dat$y <- sf::st_coordinates(sites_sf)[, 2]

project_xy <- function(lon, lat) {
  sf::st_coordinates(sf::st_transform(
    sf::st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon", "lat"),
                 crs = CRS_WGS84), CRS_EUROPE))
}
L <- row$true_vs_corrected * ARROW_MAX_DEG
a0 <- project_xy(mean(dat$longitude), mean(dat$latitude))
a1 <- project_xy(mean(dat$longitude) + cos(row$angle_rad_corrected) * L,
                 mean(dat$latitude) + sin(row$angle_rad_corrected) * L)

rng <- range(c(grid$pred, dat$Count))
world_draw <- sf::st_simplify(world_crop, dTolerance = 3000)

p_map <- ggplot2::ggplot() +
  ggplot2::geom_raster(data = grid, ggplot2::aes(x = x, y = y, fill = pred)) +
  ggplot2::geom_contour(data = grid, ggplot2::aes(x = x, y = y, z = pred),
                        colour = "white", linewidth = 0.22, alpha = 0.75,
                        bins = 10) +
  ggplot2::geom_sf(data = world_draw, fill = NA, colour = "grey30",
                   linewidth = 0.22) +
  ggplot2::geom_point(data = dat, ggplot2::aes(x = x, y = y, fill = Count),
                      shape = 21, size = 2.4, stroke = 0.45, colour = "grey15") +
  ggplot2::scale_fill_gradientn(colours = HEAT_COLOURS, limits = rng,
                                name = "Copy number") +
  ggplot2::geom_segment(
    data = data.frame(x = a0[1], y = a0[2], xend = a1[1], yend = a1[2]),
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
    arrow = grid::arrow(length = grid::unit(0.2, "cm"), type = "closed"),
    colour = "black", linewidth = 1.1, lineend = "round") +
  ggplot2::coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
                    ylim = c(bb["ymin"], bb["ymax"]), expand = FALSE,
                    crs = sf::st_crs(CRS_EUROPE)) +
  ggplot2::theme_void(base_size = 10) +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "#EDF1F4", colour = NA),
    panel.border = ggplot2::element_rect(fill = NA, colour = "grey35",
                                         linewidth = 0.4),
    plot.background = ggplot2::element_rect(fill = "white", colour = NA),
    plot.title = ggplot2::element_text(face = "bold", size = 11),
    plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30")) +
  ggplot2::labs(
    title = TE,
    subtitle = sprintf("%s | VS = %.2f | gradient towards %s",
                       row$superfamily, row$true_vs_corrected, row$direction))

save_plot(p_map, paste0(TE, "_prediction_map"), width = 7, height = 5.6)

# =============================================================================
# 4. Panel B - the climate response itself, per group
# =============================================================================
# This is the actual evidence behind the strict rule: the same variable, the
# same direction, in each of the three groups. Fits are the per-group
# quasi-Poisson models of Model A, drawn over each group's own climate range.
resp <- purrr::map_dfr(strict_vars, function(v) {
  dat %>% dplyr::transmute(group, Count, value = .data[[v]], variable = v)
})
resp$variable <- factor(resp$variable, levels = strict_vars,
                        labels = ifelse(strict_vars %in% names(BIO_LABEL),
                                        BIO_LABEL[strict_vars], strict_vars))

p_resp <- ggplot2::ggplot(resp, ggplot2::aes(x = value, y = Count,
                                             colour = group, fill = group)) +
  ggplot2::geom_point(size = 2, alpha = 0.85, shape = 21, colour = "grey20") +
  ggplot2::geom_smooth(method = "glm", se = TRUE,
                       method.args = list(family = stats::quasipoisson()),
                       linewidth = 0.8, alpha = 0.15) +
  ggplot2::facet_wrap(~ variable, scales = "free_x", ncol = 1) +
  ggplot2::scale_colour_manual(values = GROUP_COLOURS, name = NULL) +
  ggplot2::scale_fill_manual(values = GROUP_COLOURS, name = NULL) +
  ggplot2::labs(
    x = NULL, y = sprintf("%s copy number", TE),
    title = "Climate response, fitted separately in each group",
    subtitle = "Quasi-Poisson fits; the slope is negative in all three groups for both variables") +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "bottom",
                 strip.background = ggplot2::element_rect(fill = "grey95",
                                                          colour = "grey70"),
                 strip.text = ggplot2::element_text(size = 8.5),
                 plot.title = ggplot2::element_text(face = "bold", size = 11),
                 plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30"),
                 plot.background = ggplot2::element_rect(fill = "white",
                                                         colour = NA))

save_plot(p_resp, paste0(TE, "_climate_response"), width = 6, height = 7)

# =============================================================================
# 5. Composite
# =============================================================================
p_all <- (p_map | p_resp) + patchwork::plot_layout(widths = c(1.25, 1))
save_plot(p_all, paste0(TE, "_composite"), width = 13, height = 6.4)

message(sprintf("Copy number range: %d to %d observed; %d to %d predicted",
                min(dat$Count), max(dat$Count),
                round(min(grid$pred)), round(max(grid$pred))))

# =============================================================================
# 6. Session log
# =============================================================================
log_session("09_focal_TE_evidence")
message("=== 09_focal_TE_evidence.R finished ===")
