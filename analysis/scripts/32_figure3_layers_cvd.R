# =============================================================================
# Script:  32_figure3_layers_cvd.R
# Purpose: Figure 3 as separate LAYERS, for manual assembly.
#
#          Same data, same colours and same geometry as 31_figure3_cvd.R. Each
#          file holds a single kind of element on an identical canvas, so they
#          stack exactly.
#
#          How the geometry is guaranteed: the panel A layers are drawn with
#          theme_void, NO legend, NO margin and NO tag, using the same coord_sf
#          limits. The panel therefore fills the whole canvas in every file and
#          the alignment is true by construction rather than by luck. The
#          legend and the "A" tag live in files of their own.
#
#          The text stays REAL text (svglite), editable in Illustrator.
#
# NOTE:    nothing written by 31_figure3_cvd.R is overwritten; all names are new.
# Author:  Antho
# =============================================================================
source("scripts/00_setup.R")
suppressMessages({library(dplyr); library(readr); library(ggplot2)
                  library(sf); library(grid); library(gtable); library(svglite)})

# Output directory, overridable with FIG_OUT.
OUT <- file.path(Sys.getenv("FIG_OUT", unset = file.path("..", "revision_finale", "figures_cvd")), "FIG3_calques")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

SF_COL <- c(LTR_Copia = "#A50F15", LTR_Gypsy = "#CCE6FF", LTR = "#E69F00",
            noCat     = "#6E6E6E", TIR       = "#332288")
SF_LAB <- c(LTR_Copia = "Copia", LTR_Gypsy = "Gypsy", LTR = "LTR",
            noCat = "noCat", TIR = "TIR")
GRP_COL <- c(SYLW = "#5599FF", SYLE = "#FFD6A5", SYL_Admix = "#000000")
LAND <- "#F4F1EA"; SEA <- "#CBD5DB"; CASING <- "#3A3A3A"

GROUP_OF <- function(x) dplyr::case_when(
  x == "Q1" ~ "SYLW", x == "Q5" ~ "SYLE",
  x %in% c("Q1_Admix", "Q5_Admix") ~ "SYL_Admix", TRUE ~ NA_character_)

mem  <- read_csv(file.path(dir_outputs, "18_climate_pc_diagnostics", "tables",
                           "candidate_membership.csv"), show_col_types = FALSE)
CAND <- mem$TE_ID[mem$in_modelA]
arr  <- read_csv(file.path(dir_outputs, "07_map_figure3A", "tables",
                           "arrow_geometry.csv"), show_col_types = FALSE) %>%
  filter(TE_ID %in% CAND)
sites <- read_csv(file.path(dir_data, "250923_sample_list_TE.csv"),
                  show_col_types = FALSE) %>%
  filter(Species == "M. sylvestris", !is.na(latitude), !is.na(longitude)) %>%
  mutate(group = GROUP_OF(group_structure_CORE_admix_K6)) %>% filter(!is.na(group))
sxy <- st_coordinates(st_transform(
  st_as_sf(sites, coords = c("longitude", "latitude"), crs = 4326), 3035))
sites$x <- sxy[, 1]; sites$y <- sxy[, 2]
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(3035)
pad  <- 210e3
xlim <- c(min(sites$x) - pad, max(sites$x) + pad)
ylim <- c(min(sites$y) - pad, max(sites$y) + pad)

## The basemap is CROPPED to the extent of the figure. Without that a clipping
## path would be needed to hide the rest of the world, and Illustrator turns it
## into a clipping mask that has to be released before anything can be edited.
## Cropping upstream means no mask at all in the SVG.
world_crop <- suppressWarnings(st_crop(
  world, xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2]))

## Canvas shared by EVERY panel A layer.
W <- 5.2; H <- 4.6
ARROW <- grid::arrow(length = unit(.14, "cm"), type = "closed")

## Skeleton: identical coord_sf, no decoration, transparent background.
skeleton <- function(...) {
  ggplot() + list(...) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = 3035) +
    scale_colour_manual(values = SF_COL, labels = SF_LAB, guide = "none") +
    scale_fill_manual(values = GRP_COL, guide = "none") +
    scale_linewidth_continuous(range = c(.4, 1.6), guide = "none",
                               limits = range(c(arr$true_vs_corrected,
                                                arr$true_vs_corrected * .55))) +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0),
          panel.background = element_blank(), plot.background = element_blank())
}

save_svg <- function(p, name, w = W, h = H) {
  ggsave(file.path(OUT, paste0(name, ".svg")), p, width = w, height = h,
         device = svglite::svglite, bg = "transparent")
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h,
         device = grDevices::cairo_pdf, bg = "transparent")
  message("  ", name, ".svg / .pdf")
}

## --- A1: the sea alone ------------------------------------------------------
save_svg(skeleton(annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf,
                           ymax = Inf, fill = SEA)), "A1_mer")

## --- A2: countries alone (land + borders), transparent background -----------
save_svg(skeleton(geom_sf(data = world_crop, fill = LAND, colour = "#BFB8AC",
                          linewidth = .25)), "A2_pays")

## --- A3: sampling localities ------------------------------------------------
save_svg(skeleton(geom_point(data = sites, aes(x, y, fill = group), shape = 21,
                             size = 1.6, colour = "grey20", stroke = .3,
                             alpha = .95)), "A3_sites")

## --- A4: dark casing under the arrows ---------------------------------------
## A file of its own so that it can be deleted in one go if the final rendering
## no longer has a contrast problem.
save_svg(skeleton(geom_segment(data = arr,
           aes(x, y, xend = xend, yend = yend, linewidth = true_vs_corrected),
           colour = CASING, arrow = ARROW, lineend = "round", alpha = .95)),
         "A4_fleches_lisere")

## --- A5: the arrows in colour -----------------------------------------------
save_svg(skeleton(geom_segment(data = arr,
           aes(x, y, xend = xend, yend = yend, colour = superfamily,
               linewidth = true_vs_corrected * .55),
           arrow = ARROW, lineend = "round")),
         "A5_fleches_couleur")

## --- A6: the common origin --------------------------------------------------
save_svg(skeleton(geom_point(data = arr[1, ], aes(x, y), shape = 21, size = 2.2,
                             fill = "white", colour = "grey15", stroke = .6)),
         "A6_origine")

## --- A7: the legend alone ---------------------------------------------------
p_leg <- ggplot() +
  geom_point(data = sites, aes(x, y, fill = group), shape = 21, size = 1.6,
             colour = "grey20", stroke = .3) +
  geom_segment(data = arr, aes(x, y, xend = xend, yend = yend,
                               colour = superfamily), arrow = ARROW) +
  scale_colour_manual(values = SF_COL, labels = SF_LAB, name = "TE superfamily") +
  scale_fill_manual(values = GRP_COL, name = "Sampling group") +
  guides(colour = guide_legend(order = 1, nrow = 2,
                               override.aes = list(linewidth = 1.2, alpha = 1)),
         fill = guide_legend(order = 2, nrow = 1,
                             override.aes = list(size = 2.2))) +
  theme_void(base_size = 8) +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.box.just = "left", legend.spacing.y = unit(.02, "cm"),
        legend.key.size = unit(.26, "cm"),
        legend.text = element_text(size = 6.2),
        legend.title = element_text(size = 6.6))
g   <- ggplotGrob(p_leg)
idx <- grep("guide-box", g$layout$name)
idx <- idx[vapply(idx, function(i) !inherits(g$grobs[[i]], "zeroGrob"), logical(1))]
leg <- g$grobs[[idx[1]]]
for (dev in c("svg", "pdf")) {
  f <- file.path(OUT, paste0("A7_legende.", dev))
  if (dev == "svg") svglite::svglite(f, width = 3.6, height = 0.85, bg = "transparent")
  else grDevices::cairo_pdf(f, width = 3.6, height = 0.85, bg = "transparent")
  grid.newpage(); grid.draw(leg); dev.off()
}
message("  A7_legende.svg / .pdf")

## --- strip the clipping masks -----------------------------------------------
## Everything is already inside the extent (the basemap was cropped upstream),
## so the svglite clipPaths serve no purpose and only complicate editing.
for (f in list.files(OUT, pattern = "\\.svg$", full.names = TRUE)) {
  x <- paste(readLines(f, warn = FALSE), collapse = "\n")
  n <- lengths(regmatches(x, gregexpr("clip-path=", x)))
  x <- gsub("<defs>\\s*<clipPath[^>]*>.*?</clipPath>\\s*</defs>", "", x, perl = TRUE)
  x <- gsub(" clip-path='url\\(#[^']*\\)'", "", x, perl = TRUE)
  writeLines(x, f)
  message(sprintf("  %s : %d masque(s) retire(s)", basename(f), n))
}

message(sprintf("Calques du panneau A ecrits dans %s", OUT))
message(sprintf("  toile commune %.1f x %.1f in | %d pays, %d sites, %d fleches",
                W, H, nrow(world_crop), nrow(sites), nrow(arr)))
