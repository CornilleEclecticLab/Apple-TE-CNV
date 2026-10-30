# =============================================================================
# Script:  08_per_TE_surfaces.R
# Purpose: One map per candidate TE showing the modelled copy-number surface
#          across Europe as a filled colour field, weather-map style, with the
#          gradient arrow of Fig. 3A drawn on top.
#
#          The spatial model already contains everything needed for this. It
#          fits, per TE,
#              Count ~ latitude + longitude        (negative binomial, log link)
#          so it predicts an expected copy number at ANY coordinate, not only
#          at the sampled localities. Colouring that prediction over a grid is
#          simply showing the same model the arrow summarises: the arrow is the
#          direction of steepest increase of the surface drawn underneath it.
#
#          IMPORTANT INTERPRETATION LIMIT, stated on every figure. The model is
#          log-linear in latitude and longitude, so the predicted surface is
#          monotonic and its contours are straight parallel lines. These maps
#          therefore show the fitted large-scale TREND, not local structure -
#          they cannot reveal pockets of high or low copy number. To let the
#          reader see where the data depart from that trend, the observed
#          counts are drawn as points on exactly the same colour scale, so a
#          locality that disagrees with the surface is visible as a point whose
#          colour differs from its background.
#
#          The surface is clipped to a convex hull of the sampling localities
#          buffered by 200 km, intersected with land, so the log-linear trend
#          is never extrapolated far beyond where it was estimated.
#
# NAMING NOTE: in the manuscript, "Model A" is the climate GLM
# (TE count ~ bioclimatic variable) and "Model B" is the SPATIAL model
# (TE count ~ latitude + longitude, negative binomial), which produces the
# vector-strength (VS) values shown in Fig. 3A-B. To avoid a clash with that
# established terminology, the PC-covariate version of the climate model is
# called "Model A-PC" throughout this folder, never "Model B".
#
# Author:  Antho (analysis run with Claude Code)
# Date:    2026-08-24
# Project: TE_Project / review (reviewer response, Genome Biology submission)
#
# Inputs:
#   data/coord_syl_TE.csv           residual TE counts, M. sylvestris
#   data/250923_sample_list_TE.csv  coordinates and K6 clusters
#   outputs/03_spatial_vs/tables/VS_table_38_candidates.csv
#
# Outputs (all under outputs/08_per_TE_surfaces/):
#   figures/per_TE/<TE_ID>.pdf|.svg     one map per candidate TE
#   figures/all_TE_prediction_maps.pdf  every map in one multi-page PDF
#   figures/overview_top12_VS.pdf|.svg  the 12 strongest, relative scale
#   figures/panel_both_models_TEs.pdf|.svg  figure-ready panel, both-model TEs
#   tables/prediction_summary.csv       fold change across the sampled range
#   ../../logs/08_per_TE_surfaces_*_sessionInfo.txt
# =============================================================================

# --- Setup -------------------------------------------------------------------
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(sf)
  library(rnaturalearth)
  library(glmmTMB)
  library(patchwork)
  library(svglite)
})

set.seed(20260824)

# --- Parameters --------------------------------------------------------------
CRS_EUROPE   <- 3035
CRS_WGS84    <- 4326
GRID_KM      <- 11       # prediction grid resolution
HULL_BUFFER_M<- 200000   # how far beyond the sampling hull to extrapolate
ARROW_MAX_DEG<- 5.2      # same convention as script 07
N_OVERVIEW   <- 12       # TEs shown in the overview panel

# Weather-map palette: blue (low) through yellow to red (high).
HEAT_COLOURS <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))

dir_out    <- file.path(dir_outputs, "08_per_TE_surfaces")
dir_fig    <- file.path(dir_out, "figures")
dir_perte  <- file.path(dir_fig, "per_TE")
dir_tab    <- file.path(dir_out, "tables")
for (d in c(dir_out, dir_fig, dir_perte, dir_tab)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

message("=== 08_per_TE_surfaces.R ===")

# =============================================================================
# 1. Data
# =============================================================================
vs <- readr::read_csv(
  file.path(dir_outputs, "03_spatial_vs", "tables",
            "VS_table_38_candidates.csv"),
  show_col_types = FALSE
)

counts <- readr::read_csv(file.path(dir_data, "coord_syl_TE.csv"),
                          show_col_types = FALSE) %>%
  dplyr::mutate(Sample_ID = toupper(Sample_ID))

meta <- readr::read_csv(file.path(dir_data, "250923_sample_list_TE.csv"),
                        show_col_types = FALSE) %>%
  dplyr::filter(Species == "M. sylvestris",
                !is.na(latitude), !is.na(longitude),
                group_structure_CORE_admix_K6 != "Q2_Admix") %>%
  dplyr::transmute(Sample_ID = toupper(Sample_ID),
                   latitude  = as.numeric(latitude),
                   longitude = as.numeric(longitude))

sites <- dplyr::inner_join(counts %>% dplyr::select(Sample_ID), meta,
                           by = "Sample_ID")
message(sprintf("Sampling localities: %d individuals", nrow(sites)))

origin_lon <- mean(sites$longitude)
origin_lat <- mean(sites$latitude)

# =============================================================================
# 2. Prediction grid, clipped to the sampled region
# =============================================================================
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  sf::st_transform(CRS_EUROPE) %>%
  sf::st_make_valid()

sites_sf <- sf::st_as_sf(sites, coords = c("longitude", "latitude"),
                         crs = CRS_WGS84) %>%
  sf::st_transform(CRS_EUROPE)

# The trend is only credible where it was estimated, so the surface is clipped
# to the sampling hull plus a buffer, and to land.
hull <- sf::st_convex_hull(sf::st_union(sites_sf)) %>%
  sf::st_buffer(HULL_BUFFER_M)

world_crop <- suppressWarnings(
  sf::st_crop(world, sf::st_bbox(sf::st_buffer(hull, 250000))))
land <- sf::st_union(world_crop) %>% sf::st_make_valid()
region <- suppressWarnings(sf::st_intersection(hull, land)) %>% sf::st_make_valid()

bb <- sf::st_bbox(hull)
gx <- seq(bb["xmin"], bb["xmax"], by = GRID_KM * 1000)
gy <- seq(bb["ymin"], bb["ymax"], by = GRID_KM * 1000)
grid <- expand.grid(x = gx, y = gy)
message(sprintf("Prediction grid: %d x %d = %d cells at %d km",
                length(gx), length(gy), nrow(grid), GRID_KM))

grid_sf <- sf::st_as_sf(grid, coords = c("x", "y"), crs = CRS_EUROPE)
inside <- lengths(sf::st_intersects(grid_sf, region)) > 0
grid <- grid[inside, , drop = FALSE]
grid_sf <- grid_sf[inside, ]

# Back-transform to lon/lat: the model is parameterised in degrees.
ll <- sf::st_coordinates(sf::st_transform(grid_sf, CRS_WGS84))
grid$longitude <- ll[, 1]
grid$latitude  <- ll[, 2]
message(sprintf("Cells retained inside the sampled region: %d", nrow(grid)))

# Coastline for drawing, simplified so the per-TE SVG files stay small.
world_draw <- sf::st_simplify(world_crop, dTolerance = 3000)
xlim <- c(bb["xmin"], bb["xmax"])
ylim <- c(bb["ymin"], bb["ymax"])

sites$x <- sf::st_coordinates(sites_sf)[, 1]
sites$y <- sf::st_coordinates(sites_sf)[, 2]

# =============================================================================
# 3. Fit, predict and draw, one TE at a time
# =============================================================================
project_xy <- function(lon, lat) {
  p <- sf::st_as_sf(data.frame(lon = lon, lat = lat),
                    coords = c("lon", "lat"), crs = CRS_WGS84) %>%
    sf::st_transform(CRS_EUROPE)
  sf::st_coordinates(p)
}

theme_predmap <- ggplot2::theme_void(base_size = 10) +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "#EDF1F4", colour = NA),
    panel.border = ggplot2::element_rect(fill = NA, colour = "grey35",
                                         linewidth = 0.4),
    legend.position = "right",
    plot.title = ggplot2::element_text(face = "bold", size = 11),
    plot.subtitle = ggplot2::element_text(size = 8.5, colour = "grey30"),
    plot.caption = ggplot2::element_text(size = 6.8, colour = "grey40",
                                         hjust = 0),
    # ggsave exports a transparent background by default, which renders as
    # black in some viewers and is awkward to place in Illustrator.
    plot.background = ggplot2::element_rect(fill = "white", colour = NA),
    plot.margin = ggplot2::margin(8, 8, 8, 8))

#' Build the prediction map for one TE.
build_te_map <- function(te_id, clean = FALSE) {
  row <- vs[vs$TE_ID == te_id, ]

  d <- counts %>%
    dplyr::select(Sample_ID, Count = dplyr::all_of(te_id)) %>%
    dplyr::inner_join(meta, by = "Sample_ID") %>%
    tidyr::drop_na(Count, latitude, longitude)

  fit <- try(glmmTMB::glmmTMB(Count ~ latitude + longitude, data = d,
                              family = glmmTMB::nbinom2()), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)

  g <- grid
  g$pred <- as.numeric(stats::predict(fit, newdata = g, type = "response"))

  obs <- d %>%
    dplyr::left_join(sites %>% dplyr::select(Sample_ID, x, y), by = "Sample_ID")

  # Observed and predicted share one colour scale, so a locality that departs
  # from the fitted trend stands out against its background.
  rng <- range(c(g$pred, obs$Count), na.rm = TRUE)

  arrow_len <- row$true_vs_corrected * ARROW_MAX_DEG
  a_start <- project_xy(origin_lon, origin_lat)
  a_end   <- project_xy(origin_lon + cos(row$angle_rad_corrected) * arrow_len,
                        origin_lat + sin(row$angle_rad_corrected) * arrow_len)

  ggplot2::ggplot() +
    ggplot2::geom_raster(data = g, ggplot2::aes(x = x, y = y, fill = pred)) +
    ggplot2::geom_contour(data = g, ggplot2::aes(x = x, y = y, z = pred),
                          colour = "white", linewidth = 0.22, alpha = 0.75,
                          bins = 10) +
    ggplot2::geom_sf(data = world_draw, fill = NA, colour = "grey30",
                     linewidth = 0.22) +
    ggplot2::geom_point(data = obs,
                        ggplot2::aes(x = x, y = y, fill = Count),
                        shape = 21, size = 2.2, stroke = 0.45,
                        colour = "grey15") +
    ggplot2::scale_fill_gradientn(
      colours = HEAT_COLOURS, limits = rng,
      name = "Copy number", labels = scales::label_number(big.mark = " ")) +
    ggplot2::geom_segment(
      data = data.frame(x = a_start[1], y = a_start[2],
                        xend = a_end[1], yend = a_end[2]),
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      arrow = grid::arrow(length = grid::unit(0.2, "cm"), type = "closed"),
      colour = "black", linewidth = 1.1, lineend = "round") +
    ggplot2::geom_point(
      data = data.frame(x = a_start[1], y = a_start[2]),
      ggplot2::aes(x = x, y = y), shape = 21, size = 2, fill = "white",
      colour = "black", stroke = 0.5) +
    ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE,
                      crs = sf::st_crs(CRS_EUROPE)) +
    theme_predmap +
    ggplot2::labs(
      title = te_id,
      subtitle = sprintf(
        "%s | VS = %.2f | gradient towards %s | predicted %s to %s copies",
        row$superfamily, row$true_vs_corrected, row$direction,
        format(round(min(g$pred)), big.mark = " "),
        format(round(max(g$pred)), big.mark = " ")),
      # `clean = TRUE` drops the methods caption: the figure-ready panel below
      # is meant to sit inside Figure 3, whose legend carries that text.
      caption = if (clean) NULL else paste0(
        "Fill and contours: expected copy number from the spatial model ",
        "(negative binomial, log link, latitude + longitude).\n",
        "The model is log-linear, so contours are straight: this is the fitted ",
        "large-scale trend, not local structure.\n",
        "Points are observed counts on the same scale. Clipped to the sampling ",
        "hull + 200 km. Projection EPSG:3035."))
}

vs <- vs %>% dplyr::arrange(dplyr::desc(true_vs_corrected))
maps <- list()
summ <- list()

for (i in seq_len(nrow(vs))) {
  te <- vs$TE_ID[i]
  p <- build_te_map(te)
  if (is.null(p)) { message("  model failed for ", te); next }
  maps[[te]] <- p

  ggplot2::ggsave(file.path(dir_perte, paste0(te, ".pdf")), p,
                  width = 7.2, height = 6, device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(dir_perte, paste0(te, ".svg")), p,
                  width = 7.2, height = 6, device = svglite::svglite)

  g <- p$layers[[1]]$data
  summ[[te]] <- tibble::tibble(
    TE_ID = te,
    superfamily = vs$superfamily[i],
    true_vs_corrected = vs$true_vs_corrected[i],
    direction = vs$direction[i],
    pred_min = min(g$pred), pred_max = max(g$pred),
    fold_change_across_range = max(g$pred) / min(g$pred)
  )
  if (i %% 10 == 0) message(sprintf("  %d / %d maps done", i, nrow(vs)))
}

message(sprintf("Per-TE maps written: %d", length(maps)))

summary_tbl <- dplyr::bind_rows(summ) %>%
  dplyr::arrange(dplyr::desc(fold_change_across_range))
readr::write_csv(summary_tbl, file.path(dir_tab, "prediction_summary.csv"))
message("Strongest predicted gradients across the sampled range:")
print(as.data.frame(utils::head(summary_tbl %>%
  dplyr::mutate(fold_change_across_range = round(fold_change_across_range, 2),
                pred_min = round(pred_min), pred_max = round(pred_max)), 8)))

# =============================================================================
# 4. All maps in one multi-page PDF, for browsing
# =============================================================================
multipage <- file.path(dir_fig, "all_TE_prediction_maps.pdf")
grDevices::cairo_pdf(multipage, width = 7.2, height = 6, onefile = TRUE)
for (te in names(maps)) print(maps[[te]])
grDevices::dev.off()
message("Multi-page PDF written: ", multipage)

# =============================================================================
# 5. Overview panel - the strongest gradients, on a relative scale
# =============================================================================
# Absolute copy number differs by orders of magnitude between TEs, so a shared
# absolute scale would be unreadable. Each panel is rescaled to its own range,
# which is what the eye compares here: the shape of the gradient, not the level.
top_te <- utils::head(vs$TE_ID, N_OVERVIEW)

ov <- purrr::map_dfr(top_te, function(te) {
  d <- counts %>%
    dplyr::select(Sample_ID, Count = dplyr::all_of(te)) %>%
    dplyr::inner_join(meta, by = "Sample_ID") %>%
    tidyr::drop_na(Count)
  fit <- try(glmmTMB::glmmTMB(Count ~ latitude + longitude, data = d,
                              family = glmmTMB::nbinom2()), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  g <- grid
  g$pred <- as.numeric(stats::predict(fit, newdata = g, type = "response"))
  g$rel  <- (g$pred - min(g$pred)) / (max(g$pred) - min(g$pred))
  r <- vs[vs$TE_ID == te, ]
  g$panel <- sprintf("%s\n%s | VS = %.2f", te, r$superfamily,
                     r$true_vs_corrected)
  g
})
ov$panel <- factor(ov$panel, levels = unique(ov$panel))

arrows_ov <- purrr::map_dfr(top_te, function(te) {
  r <- vs[vs$TE_ID == te, ]
  L <- r$true_vs_corrected * ARROW_MAX_DEG
  s <- project_xy(origin_lon, origin_lat)
  e <- project_xy(origin_lon + cos(r$angle_rad_corrected) * L,
                  origin_lat + sin(r$angle_rad_corrected) * L)
  tibble::tibble(panel = sprintf("%s\n%s | VS = %.2f", te, r$superfamily,
                                 r$true_vs_corrected),
                 x = s[1], y = s[2], xend = e[1], yend = e[2])
})
arrows_ov$panel <- factor(arrows_ov$panel, levels = levels(ov$panel))

p_ov <- ggplot2::ggplot() +
  ggplot2::geom_raster(data = ov, ggplot2::aes(x = x, y = y, fill = rel)) +
  ggplot2::geom_sf(data = world_draw, fill = NA, colour = "grey35",
                   linewidth = 0.18) +
  ggplot2::geom_segment(
    data = arrows_ov,
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
    arrow = grid::arrow(length = grid::unit(0.13, "cm"), type = "closed"),
    colour = "black", linewidth = 0.7, lineend = "round") +
  ggplot2::scale_fill_gradientn(
    colours = HEAT_COLOURS,
    name = "Predicted copy number\n(rescaled within TE)",
    breaks = c(0, 0.5, 1), labels = c("lowest", "mid", "highest")) +
  ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE,
                    crs = sf::st_crs(CRS_EUROPE)) +
  ggplot2::facet_wrap(~ panel, ncol = 4) +
  theme_predmap +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 7.5, lineheight = 1.1),
                 legend.position = "right") +
  ggplot2::labs(
    title = sprintf("Modelled copy-number surfaces, the %d strongest spatial gradients",
                    N_OVERVIEW),
    subtitle = "Each panel rescaled to its own range; arrows as in Fig. 3A",
    caption = paste0(
      "Log-linear spatial model, so each surface is a monotonic trend. ",
      "Clipped to the sampling hull + 200 km."))

ggplot2::ggsave(file.path(dir_fig, "overview_top12_VS.pdf"), p_ov,
                width = 13, height = 9.5, device = grDevices::cairo_pdf)
ggplot2::ggsave(file.path(dir_fig, "overview_top12_VS.svg"), p_ov,
                width = 13, height = 9.5, device = svglite::svglite)

# =============================================================================
# 6. Figure-ready panel: the TEs supported by BOTH candidate models
# =============================================================================
# These are the highest-confidence candidates - significant for at least one
# bioclimatic variable in all three M. sylvestris groups, both with and without
# the neutral SNP-PC covariates. Only a handful qualify, so each panel keeps
# its own absolute colour bar (copy number differs by an order of magnitude
# between them, and a shared scale would flatten one of the two).
final_union <- readr::read_csv(
  file.path(dir_outputs, "02_candidate_TEs", "tables",
            "FINAL_union_candidates.csv"),
  show_col_types = FALSE
)
both_te <- final_union %>%
  dplyr::filter(in_model_A, in_model_A_PC) %>%
  dplyr::pull(TE_ID)

message(sprintf("TEs supported by both models: %d (%s)",
                length(both_te), paste(both_te, collapse = ", ")))

if (length(both_te) > 0) {
  panels <- lapply(both_te, function(te) {
    build_te_map(te, clean = TRUE) +
      ggplot2::theme(
        legend.key.height = grid::unit(0.85, "cm"),
        legend.key.width  = grid::unit(0.32, "cm"),
        legend.title = ggplot2::element_text(size = 8),
        legend.text  = ggplot2::element_text(size = 7),
        plot.title    = ggplot2::element_text(face = "bold", size = 10),
        plot.subtitle = ggplot2::element_text(size = 7.6, colour = "grey30"))
  })

  p_both <- patchwork::wrap_plots(panels, nrow = 1)

  w <- 6.4 * length(panels)
  ggplot2::ggsave(file.path(dir_fig, "panel_both_models_TEs.pdf"), p_both,
                  width = w, height = 5.4, device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(dir_fig, "panel_both_models_TEs.svg"), p_both,
                  width = w, height = 5.4, device = svglite::svglite)
  message("Figure-ready panel written: panel_both_models_TEs.pdf / .svg")

  # Stacked variant, for a portrait figure column.
  p_both_v <- patchwork::wrap_plots(panels, ncol = 1)
  ggplot2::ggsave(file.path(dir_fig, "panel_both_models_TEs_stacked.pdf"),
                  p_both_v, width = 6.6, height = 5.0 * length(panels),
                  device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(dir_fig, "panel_both_models_TEs_stacked.svg"),
                  p_both_v, width = 6.6, height = 5.0 * length(panels),
                  device = svglite::svglite)

  readr::write_csv(
    summary_tbl %>% dplyr::filter(TE_ID %in% both_te),
    file.path(dir_tab, "both_models_TE_prediction_summary.csv"))
}

# =============================================================================
# 6. Session log
# =============================================================================
log_session("08_per_TE_surfaces")
message("=== 08_per_TE_surfaces.R finished ===")
