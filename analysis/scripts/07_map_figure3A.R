# =============================================================================
# Script:  07_map_figure3A.R
# Purpose: Redraw Figure 3A - the geographical projection of the spatial
#          gradients of the candidate residual TEs in M. sylvestris.
#
#          Two problems with the published version are fixed here.
#
#          1. ARROW LENGTH CARRIED NO INFORMATION. The original code computed
#                 end_lat = center_lat + (beta_lat / norm) * 1.5
#             i.e. every vector was rescaled to the same fixed length, so the
#             arrows encoded direction only. (A `scaled_norm` column was
#             computed but never used.) The figure legend nonetheless states
#             that "vector direction and length indicate the orientation and
#             strength of the spatial gradient". Here arrow length is
#             genuinely proportional to the vector strength VS
#             (`true_vs_corrected`), the same quantity plotted in Fig. 3B, and
#             linewidth and opacity reinforce it.
#
#          2. NO MAP PROJECTION. The original used coord_cartesian() on raw
#             longitude/latitude, which stretches Europe horizontally. Here
#             the map is drawn in ETRS89 / LAEA Europe (EPSG:3035), the
#             standard equal-area projection for the continent, so distances
#             and directions look right.
#
#          Arrow geometry stays faithful to the published definition: the
#          angle is atan2(beta_lat, beta_long) in degree space, sign-corrected
#          exactly as in the spatial model, and the endpoint is built in
#          lon/lat before being projected, so the published VS values and
#          directions are the ones being drawn.
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
#   outputs/03_spatial_vs/tables/VS_table_38_candidates.csv
#   data/250923_sample_list_TE.csv     sampling coordinates
#   rnaturalearth country polygons (bundled in rnaturalearthdata, offline)
#
# Outputs (all under outputs/07_map_figure3A/):
#   figures/Fig3A_map_arrows.pdf|.svg          the main panel
#   figures/Fig3A_map_arrows_labelled.pdf|.svg with TE labels on strong vectors
#   figures/Fig3A_windrose.pdf|.svg            direction rose weighted by VS
#   figures/Fig3A_map_by_superfamily.pdf|.svg  small multiples
#   figures/Fig3A_composite.pdf|.svg           map + rose assembled
#   tables/arrow_geometry.csv                  plotted coordinates, for reuse
#   ../../logs/07_map_figure3A_*_sessionInfo.txt
# =============================================================================

# --- Setup -------------------------------------------------------------------
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(sf)
  library(rnaturalearth)
  library(ggrepel)
  library(patchwork)
  library(svglite)
})

set.seed(20260824)

# --- Parameters --------------------------------------------------------------
CRS_EUROPE   <- 3035      # ETRS89 / LAEA Europe - equal-area, standard for EU
CRS_WGS84    <- 4326
# Arrow length is set in projected kilometres, not in degrees. Building the
# endpoints in lon/lat and projecting them afterwards made length depend on
# direction: at 46 N a degree of longitude is about 77 km against 111 km for a
# degree of latitude, so an east-pointing arrow came out ~30% shorter than a
# north-pointing one of the same VS. Measured on the previous version, length
# per unit VS ranged from 681 to 1035 km (CV 11%). Working in the equal-area
# projection makes length exactly proportional to VS in every direction.
ARROW_MAX_KM <- 900       # length in km of an arrow with VS = 1. Raising it
                          # lengthens every arrow by the same factor, so the
                          # proportionality to VS is preserved exactly.
VS_REF       <- c(0.25, 0.5, 0.75)   # reference arrows drawn in the scale box
INSET_SPAN   <- max(VS_REF) + 0.06   # box width in VS units; labels sit above
                          # the arrows rather than after them, so the box is no
                          # wider than the longest reference arrow.
VS_STRONG    <- 0.5       # "strong" threshold used in the Results
MAP_MARGIN_KM<- 210       # padding around the sampling bounding box

# TE palette. Taken verbatim from the Figure 1 heatmaps
# (1-PCA/250925_PCA_TE.rmd, ann_colors: Family and Superfamily), so a
# superfamily carries one colour across the whole paper. The R colour names
# are kept rather than their hex codes, so this is literally the same
# specification as Figure 1.
#   Copia  lightgreen  #90EE90      (Superfamily)
#   Gypsy  lightblue   #ADD8E6      (Superfamily)
#   LTR    orange      #FFA500      (Family)
#   noCat  grey        #BEBEBE      (Family)
#   TIR    purple      #A020F0      (Family)
SF_COLOURS <- c(
  LTR_Copia = "lightgreen",
  LTR_Gypsy = "lightblue",
  LTR       = "orange",
  noCat     = "grey",
  TIR       = "purple"
)
SF_LABELS <- c(LTR_Copia = "Copia", LTR_Gypsy = "Gypsy", LTR = "LTR",
               noCat = "noCat", TIR = "TIR")

# Sampling-group palette, consistent with the rest of the project.
GROUP_COLOURS <- c(SYLW = "#2E6DA4", SYLE = "#B8A317", SYL_Admix = "#C4571F")

dir_out <- file.path(dir_outputs, "07_map_figure3A")
dir_fig <- file.path(dir_out, "figures")
dir_tab <- file.path(dir_out, "tables")
for (d in c(dir_out, dir_fig, dir_tab)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

save_plot <- function(p, name, width = 7, height = 5) {
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".pdf")), p,
                  width = width, height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".svg")), p,
                  width = width, height = height, device = svglite::svglite)
}

message("=== 07_map_figure3A.R ===")

# =============================================================================
# 1. Data
# =============================================================================
vs <- readr::read_csv(
  file.path(dir_outputs, "03_spatial_vs", "tables",
            "VS_table_38_candidates.csv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    superfamily = factor(superfamily, levels = names(SF_COLOURS)),
    sf_label    = SF_LABELS[as.character(superfamily)]
  )

message(sprintf("Candidate TEs to draw: %d", nrow(vs)))

sites <- readr::read_csv(file.path(dir_data, "250923_sample_list_TE.csv"),
                         show_col_types = FALSE) %>%
  dplyr::filter(Species == "M. sylvestris",
                !is.na(latitude), !is.na(longitude),
                group_structure_CORE_admix_K6 != "Q2_Admix") %>%
  dplyr::transmute(
    Sample_ID = toupper(Sample_ID),
    longitude = as.numeric(longitude),
    latitude  = as.numeric(latitude),
    group = dplyr::case_when(
      group_structure_CORE_admix_K6 %in% c("Q1_Admix", "Q5_Admix") ~ "SYL_Admix",
      group_structure_CORE_admix_K6 == "Q5" ~ "SYLE",
      group_structure_CORE_admix_K6 == "Q1" ~ "SYLW",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(group))

message(sprintf("Sampling sites: %d individuals, %d distinct localities",
                nrow(sites), nrow(dplyr::distinct(sites, longitude, latitude))))

# All vectors share one origin: the centroid of the sampling area. Each TE's
# gradient is a property of the whole sampled range, not of a locality, so a
# common origin is the honest representation.
origin_lon <- mean(sites$longitude)
origin_lat <- mean(sites$latitude)

# =============================================================================
# 2. Arrow geometry
# =============================================================================
# The published angle convention is mathematical (0 = East, counterclockwise),
# from atan2(beta_lat, beta_long) with the sign correction already applied in
# the spatial model. Endpoints are built in lon/lat, then projected, so the
# drawn direction is the geographic one.
arrows <- vs %>%
  dplyr::mutate(
    angle     = angle_rad_corrected,
    length_m  = true_vs_corrected * ARROW_MAX_KM * 1000,
    start_lon = origin_lon,
    start_lat = origin_lat,
    # Compass bearing, for the wind rose (0 = North, clockwise).
    bearing   = (90 - angle * 180 / pi) %% 360
  )

project_xy <- function(lon, lat) {
  p <- sf::st_as_sf(data.frame(lon = lon, lat = lat),
                    coords = c("lon", "lat"), crs = CRS_WGS84) %>%
    sf::st_transform(CRS_EUROPE)
  xy <- sf::st_coordinates(p)
  list(x = xy[, 1], y = xy[, 2])
}

# Endpoints are laid out in the projected plane, so the drawn length is
# exactly VS x ARROW_MAX_KM whatever the direction. EPSG:3035 is centred on
# 10 E and the vectors radiate from about 9 E, so meridian convergence is
# under a degree and the grid direction is the geographic one.
s_xy <- project_xy(arrows$start_lon, arrows$start_lat)
arrows$x <- s_xy$x; arrows$y <- s_xy$y
arrows$xend <- arrows$x + cos(arrows$angle) * arrows$length_m
arrows$yend <- arrows$y + sin(arrows$angle) * arrows$length_m

# Back to lon/lat, so the exported table stays readable.
end_ll <- sf::st_as_sf(data.frame(x = arrows$xend, y = arrows$yend),
                       coords = c("x", "y"), crs = CRS_EUROPE) %>%
  sf::st_transform(CRS_WGS84) %>% sf::st_coordinates()
arrows$end_lon <- end_ll[, 1]; arrows$end_lat <- end_ll[, 2]

sites_xy <- project_xy(sites$longitude, sites$latitude)
sites$x <- sites_xy$x; sites$y <- sites_xy$y

readr::write_csv(
  arrows %>% dplyr::select(TE_ID, superfamily, support, true_vs_corrected,
                           direction, angle_rad_corrected, bearing,
                           start_lon, start_lat, end_lon, end_lat,
                           x, y, xend, yend),
  file.path(dir_tab, "arrow_geometry.csv")
)

# =============================================================================
# 3. Base map
# =============================================================================
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  sf::st_transform(CRS_EUROPE)

# Extent: the sampling bounding box plus a margin, and wide enough that the
# longest arrow stays inside the frame.
# The frame is set by the sampling area, not by the arrow tips: letting the
# tips drive it pushed the map far into Scandinavia and the Balkans and left
# most of the panel empty. Arrow tips are still allowed to extend the frame,
# but only by a limited amount.
pad <- MAP_MARGIN_KM * 1000
tip_slack <- 120 * 1000
xlim <- c(min(min(sites$x) - pad, min(arrows$xend) - tip_slack),
          max(max(sites$x) + pad, max(arrows$xend) + tip_slack))
ylim <- c(min(min(sites$y) - pad, min(arrows$yend) - tip_slack),
          max(max(sites$y) + pad, max(arrows$yend) + tip_slack))

# coord_sf keeps the true aspect ratio, so if the saved figure does not have
# the same shape as the data the panel pads the shorter axis with empty ocean.
# The panel aspect is therefore computed from the data and the figure height
# derived from it, rather than guessed.
LEGEND_IN <- 2.6   # horizontal space taken by the legends
CHROME_IN <- 1.5   # vertical space taken by title, subtitle and caption
map_fig_size <- function(fig_width = 9.5) {
  panel_w <- fig_width - LEGEND_IN
  panel_h <- panel_w / (diff(xlim) / diff(ylim))
  c(width = fig_width, height = panel_h + CHROME_IN)
}

theme_map <- ggplot2::theme_void(base_size = 11) +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "#EAF2F7", colour = NA),
    panel.border     = ggplot2::element_rect(fill = NA, colour = "grey35",
                                             linewidth = 0.4),
    panel.grid.major = ggplot2::element_line(colour = "#FFFFFF", linewidth = 0.25),
    legend.position  = "right",
    legend.key       = ggplot2::element_blank(),
    plot.title       = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle    = ggplot2::element_text(size = 9, colour = "grey30"),
    plot.caption     = ggplot2::element_text(size = 7.5, colour = "grey35",
                                             hjust = 0)
  )

base_layers <- list(
  ggplot2::geom_sf(data = world, fill = "#F4F1EA", colour = "#BFB8AC",
                   linewidth = 0.25),
  ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE,
                    crs = sf::st_crs(CRS_EUROPE))
)

# =============================================================================
# 4. Main panel
# =============================================================================
#' Reference arrows showing what a given VS looks like, drawn as a small
#' standalone plot and inset into the map. ggplot cannot legend a length
#' aesthetic directly, so this is built explicitly.
#' Fraction of the panel width taken by an arrow of VS = 1, measured on the
#' projected map. The inset is sized from this so that a reference arrow is
#' physically the same length as a map arrow of the same VS: the box is a real
#' scale bar, not a relative one.
vs_unit_fraction <- function() ARROW_MAX_KM * 1000 / diff(xlim)

make_scale_inset <- function(values = VS_REF) {
  ref <- tibble::tibble(
    vs = values,
    y  = rev(seq_along(values)),
    x  = 0,
    xend = values * 1.0
  )
  ggplot2::ggplot(ref) +
    ggplot2::geom_segment(
      ggplot2::aes(x = x, y = y, xend = xend, yend = y, linewidth = vs),
      arrow = grid::arrow(length = grid::unit(0.16, "cm"), type = "closed"),
      colour = "grey20", lineend = "round") +
    ggplot2::geom_text(ggplot2::aes(x = xend, y = y + 0.33,
                                    label = sprintf("%.2f", vs)),
                       hjust = 1, size = 2.5, colour = "grey20") +
    ggplot2::scale_linewidth_continuous(range = c(0.35, 1.5), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(-0.03, max(values) + 0.03)) +
    ggplot2::scale_y_continuous(limits = c(0.45, length(values) + 1.05)) +
    ggplot2::annotate("text", x = 0, y = length(values) + 0.85,
                      label = "VS", hjust = 0, size = 2.8,
                      fontface = "bold", colour = "grey20") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.background = ggplot2::element_rect(
      fill = scales::alpha("white", 0.82), colour = "grey60", linewidth = 0.3))
}

build_map <- function(labelled = FALSE) {
  p <- ggplot2::ggplot() +
    base_layers +
    # Sampling localities, for geographic context.
    ggplot2::geom_point(
      data = sites,
      ggplot2::aes(x = x, y = y, fill = group),
      shape = 21, size = 2.6, stroke = 0.35, colour = "white", alpha = 0.95) +
    ggplot2::scale_fill_manual(
      values = GROUP_COLOURS, name = "Sampling group",
      guide = ggplot2::guide_legend(order = 2, override.aes = list(size = 3))) +
    ggnewscale::new_scale_fill() +
    # The vectors: length, linewidth and opacity all scale with VS.
    ggplot2::geom_segment(
      data = arrows,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                   colour = superfamily, linewidth = true_vs_corrected,
                   alpha = true_vs_corrected),
      arrow = grid::arrow(length = grid::unit(0.16, "cm"), type = "closed"),
      lineend = "round") +
    ggplot2::scale_colour_manual(
      values = SF_COLOURS, labels = SF_LABELS, name = "TE superfamily",
      drop = FALSE,
      guide = ggplot2::guide_legend(order = 1,
                                    override.aes = list(linewidth = 1.4,
                                                        alpha = 1))) +
    ggplot2::scale_linewidth_continuous(range = c(0.35, 1.5), guide = "none") +
    # Opacity still rises with VS, but from a higher floor. The Figure 1
    # palette is pale (lightgreen, lightblue, grey) and at alpha 0.45 the
    # weakest vectors disappeared into the basemap. The hues are untouched;
    # only the transparency ramp is compressed, and VS remains encoded by
    # length and linewidth as well.
    ggplot2::scale_alpha_continuous(range = c(0.7, 1), guide = "none") +
    # The common origin of every vector.
    ggplot2::geom_point(
      data = data.frame(x = arrows$x[1], y = arrows$y[1]),
      ggplot2::aes(x = x, y = y), shape = 21, size = 2.4,
      fill = "white", colour = "grey15", stroke = 0.6) +
    theme_map +
    ggplot2::labs(
      title = "A  Spatial gradients of climate-associated residual TEs",
      subtitle = sprintf(
        "%d candidate TEs in %s; arrow length, width and opacity are proportional to vector strength (VS)",
        nrow(arrows), "Malus sylvestris"),
      caption = paste0(
        "Vectors radiate from the centroid of the sampling area; direction is the orientation of the modelled\n",
        "copy-number gradient (negative-binomial GLM on latitude and longitude). Points are sampling localities.\n",
        "Projection: ETRS89 / LAEA Europe (EPSG:3035)."))

  if (labelled) {
    p <- p + ggrepel::geom_text_repel(
      data = arrows %>% dplyr::filter(true_vs_corrected > VS_STRONG),
      ggplot2::aes(x = xend, y = yend, label = sub("^TE_0*", "TE", TE_ID),
                   colour = superfamily),
      size = 2.1, segment.size = 0.2, segment.colour = "grey55",
      min.segment.length = 0, max.overlaps = 40, show.legend = FALSE,
      box.padding = 0.25, seed = 20260824)
  }
  p
}

p_map <- build_map(labelled = FALSE)
p_map_inset <- p_map +
  patchwork::inset_element(make_scale_inset(), left = 0.005, bottom = 0.005,
                           right = 0.005 + INSET_SPAN * vs_unit_fraction(),
                           top = 0.235, align_to = "panel")

# Saved for the Fig. 3 assembly (17_figure3_final.R): the base map, the scale
# box and the box width, so the assembly can place the inset itself and still
# keep it a true scale bar.
saveRDS(list(map = p_map,
             aspect = diff(xlim) / diff(ylim),   # so the Fig. 3 assembly can
                                                 # size row 1 to the map and
                                                 # leave no empty band
             make_inset = make_scale_inset,
             unit_fraction = vs_unit_fraction(),
             inset_right = 0.005 + INSET_SPAN * vs_unit_fraction()),
        file.path(dir_out, "panel_A_map.rds"))

mfs <- map_fig_size(9.5)
message(sprintf("Map panel aspect %.2f -> figure %.1f x %.1f in",
                diff(xlim)/diff(ylim), mfs[["width"]], mfs[["height"]]))
save_plot(p_map_inset, "Fig3A_map_arrows", width = mfs[["width"]], height = mfs[["height"]])
save_plot(build_map(labelled = TRUE) +
            patchwork::inset_element(make_scale_inset(), left = 0.005,
                                     bottom = 0.005, right = 0.005 + INSET_SPAN * vs_unit_fraction(),
                                     top = 0.235, align_to = "panel"),
          "Fig3A_map_arrows_labelled", width = mfs[["width"]],
          height = mfs[["height"]])

# =============================================================================
# 5. Direction rose, weighted by VS
# =============================================================================
# 38 vectors from a single origin overlap heavily; the rose summarises the
# same information without occlusion and makes the southeastern bias explicit.
BIN <- 30
rose <- arrows %>%
  dplyr::mutate(bin = floor(bearing / BIN) * BIN + BIN / 2) %>%
  dplyr::group_by(bin, superfamily) %>%
  dplyr::summarise(vs_sum = sum(true_vs_corrected), n = dplyr::n(),
                   .groups = "drop")

p_rose <- ggplot2::ggplot(rose,
                          ggplot2::aes(x = bin, y = vs_sum, fill = superfamily)) +
  ggplot2::geom_col(width = BIN * 0.92, colour = "white", linewidth = 0.2) +
  ggplot2::coord_polar(start = 0, direction = 1) +
  ggplot2::scale_x_continuous(limits = c(0, 360),
                              breaks = seq(0, 315, by = 45),
                              labels = c("N", "NE", "E", "SE",
                                         "S", "SW", "W", "NW")) +
  ggplot2::scale_fill_manual(values = SF_COLOURS, labels = SF_LABELS,
                             name = "TE superfamily", drop = FALSE) +
  ggplot2::labs(
    title = "B  Direction of the gradients",
    subtitle = sprintf("%d-degree bins, stacked by summed VS", as.integer(BIN)),
    x = NULL, y = "Summed VS") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_line(colour = "grey85", linewidth = 0.3),
    panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
    axis.text.x  = ggplot2::element_text(face = "bold", size = 9),
    axis.text.y  = ggplot2::element_text(size = 7, colour = "grey45"),
    plot.title   = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle= ggplot2::element_text(size = 9, colour = "grey30"),
    legend.position = "right")

save_plot(p_rose, "Fig3A_windrose", width = 6.2, height = 5.6)

# =============================================================================
# 6. Small multiples, one facet per superfamily
# =============================================================================
facet_pad <- 90 * 1000
facet_xlim <- range(c(arrows$x, arrows$xend)) + c(-facet_pad, facet_pad)
facet_ylim <- range(c(arrows$y, arrows$yend)) + c(-facet_pad, facet_pad)

p_facets <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = world, fill = "#F4F1EA", colour = "#BFB8AC",
                   linewidth = 0.2) +
  ggplot2::geom_segment(
    data = arrows,
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                 colour = superfamily, linewidth = true_vs_corrected),
    arrow = grid::arrow(length = grid::unit(0.16, "cm"), type = "closed"),
    alpha = 0.85, lineend = "round") +
  ggplot2::scale_colour_manual(values = SF_COLOURS, guide = "none",
                               drop = FALSE) +
  ggplot2::scale_linewidth_continuous(range = c(0.4, 1.6), guide = "none") +
  # Small multiples are cropped to the vector field itself: at the full
  # sampling extent each panel is mostly empty and the arrows are illegible.
  ggplot2::coord_sf(xlim = facet_xlim, ylim = facet_ylim, expand = FALSE,
                    crs = sf::st_crs(CRS_EUROPE)) +
  ggplot2::facet_wrap(~ sf_label, nrow = 1) +
  theme_map +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 10),
                 panel.grid.major = ggplot2::element_blank()) +
  ggplot2::labs(title = "Spatial gradients by TE superfamily",
                subtitle = "Arrow length proportional to vector strength (VS)")

facet_h <- (14 - 0.6) / 5 / (diff(facet_xlim) / diff(facet_ylim)) + 1.1
save_plot(p_facets, "Fig3A_map_by_superfamily", width = 14, height = facet_h)

# =============================================================================
# 7. Composite panel
# =============================================================================
# Guides are left on each panel: the inset makes a collected legend unreliable
# across patchwork versions, and the two panels share the same colour scale.
p_rose_nolegend <- p_rose + ggplot2::theme(legend.position = "none")
p_composite <- (p_map_inset | p_rose_nolegend) +
  patchwork::plot_layout(widths = c(1.55, 1))

save_plot(p_composite, "Fig3A_composite",
          width = mfs[["width"]] + 6, height = max(mfs[["height"]], 6))

# =============================================================================
# 8. Console summary
# =============================================================================
message("Vector strength distribution:")
print(as.data.frame(
  arrows %>%
    dplyr::group_by(sf_label) %>%
    dplyr::summarise(n = dplyr::n(),
                     median_VS = round(stats::median(true_vs_corrected), 3),
                     max_VS = round(max(true_vs_corrected), 3),
                     n_strong = sum(true_vs_corrected > VS_STRONG),
                     .groups = "drop")))
message("Direction breakdown:")
print(as.data.frame(arrows %>% dplyr::count(direction, name = "n")))

# =============================================================================
# 9. Session log
# =============================================================================
log_session("07_map_figure3A")
message("=== 07_map_figure3A.R finished ===")
