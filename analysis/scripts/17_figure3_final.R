# =============================================================================
# Script:  17_figure3_final.R
# Purpose: Assemble the published Figure 3 from the panels built by the
#          analysis scripts, at journal width and in a layout that reads
#          across rather than down.
#
#          The panels are NOT rebuilt here. Each analysis script saves its
#          ggplot object (panel_A_map.rds, panel_B_vs.rds, panel_C_cpg.rds)
#          and this script only lays them out, so the published figure can
#          never drift from the analysis that produced it.
#
#          Layout, on a 180 mm page width:
#            row 1   A  map of the spatial gradients   |  B  vector strength
#            row 2   C  CpG distributions              |     dinucleotide trio
#          The two columns are aligned, so the map sits above the CpG
#          distributions and the per-TE bars above the per-dinucleotide rows.
#
#          Per-panel titles, subtitles and captions are stripped: the
#          manuscript legend carries that text, and repeating it inside the
#          artwork is what made the stacked draft hard to read.
#
# Inputs:  outputs/03_spatial_vs/panel_B_vs.rds
#          outputs/07_map_figure3A/panel_A_map.rds
#          outputs/15_sequence_composition/panel_C_cpg.rds
# Outputs: outputs/17_figure3_final/figures/Fig3.pdf + .svg
#
# Author:  Antho
# Project: TE_Project / final_pipeline (Genome Biology revision)
# =============================================================================

source("scripts/00_setup.R")
library(ggplot2)
library(patchwork)

TAG     <- "17_figure3_final"
dir_out <- file.path(dir_outputs, TAG)
dir_fig <- file.path(dir_out, "figures")
dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)

# Genome Biology full page width is 180 mm.
FIG_W_IN <- 180 / 25.4
# The column split and the row heights are derived from the map's own aspect
# ratio further down, so the map fills its cell instead of floating in a white
# band, which is what made it look small next to panel B.
COL_A_SHARE  <- 0.58      # share of the width given to the map column
LEGEND_IN    <- 0.60      # vertical room taken by the two keys under the map
ROW2_IN      <- 2.75      # height of the CpG row

need <- c(A = file.path(dir_outputs, "07_map_figure3A", "panel_A_map.rds"),
          B = file.path(dir_outputs, "03_spatial_vs", "panel_B_vs.rds"),
          C = file.path(dir_outputs, "15_sequence_composition", "panel_C_cpg.rds"))
missing <- need[!file.exists(need)]
if (length(missing)) {
  stop("Missing panel objects. Run scripts 03, 07 and 15 first:\n  ",
       paste(missing, collapse = "\n  "))
}
A <- readRDS(need[["A"]]); B <- readRDS(need[["B"]]); C <- readRDS(need[["C"]])

#' Remove the text a journal legend already carries. Deliberately does NOT
#' touch axis.title or axis.text: the map panel blanks those in its own theme,
#' and re-declaring them here resurrected the graticule labels and a stray
#' "x" / "y" on the map. Type sizes are set per panel instead.
strip_text <- function(p) {
  p + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(size = 10.5, face = "bold"),
      # Left in the default position, above the panel: placed inside, the B
      # tag landed on the first bar label.
      plot.margin = ggplot2::margin(1, 2, 1, 2))
}

# --- A : the map ------------------------------------------------------------
# Panel B lists the VS of every element, so the full three-arrow scale box is
# redundant here; one reference arrow is enough to read the map. Its width is
# recomputed from the map scale so it stays a true scale bar.
REF_VS <- 0.5
inset_right <- 0.005 + (REF_VS + 0.06) * A$unit_fraction

# The two keys are moved under the row and shared with panel B, which uses the
# same superfamily colours. Left inside panel A they took about 40% of its
# width and squeezed the map to a thumbnail.
pA <- (strip_text(A$map) + ggplot2::labs(tag = "A") +
         # Stacked, not side by side: laid out horizontally the two keys were
         # wider than the map and ran into panel B.
         ggplot2::guides(
           colour = ggplot2::guide_legend(nrow = 2, order = 1,
                                          override.aes = list(linewidth = 1.2, alpha = 1)),
           fill = ggplot2::guide_legend(nrow = 1, order = 2)) +
         ggplot2::theme(legend.position = "bottom",
                        legend.box = "vertical",
                        legend.box.just = "left",
                        legend.spacing.y = ggplot2::unit(0.02, "cm"),
                        legend.text = ggplot2::element_text(size = 6.1),
                        legend.title = ggplot2::element_text(size = 6.5),
                        legend.key.size = ggplot2::unit(0.26, "cm"),
                        legend.margin = ggplot2::margin(0, 4, 0, 0),
                        legend.box.margin = ggplot2::margin(-4, 0, 0, 0))) +
  patchwork::inset_element(A$make_inset(REF_VS),
                           left = 0.008, bottom = 0.008,
                           right = inset_right, top = 0.135,
                           align_to = "panel")

# --- B : vector strength ----------------------------------------------------
# 38 rows in half a page: the labels have to be small. The superfamily key is
# on panel A, so the legend goes.
pB <- strip_text(B) + ggplot2::labs(tag = "B") +
  ggplot2::theme(legend.position = "none",
                 axis.title.x = ggplot2::element_text(size = 7.4),
                 axis.text.x = ggplot2::element_text(size = 6.4),
                 axis.text.y = ggplot2::element_text(size = 4.2),
                 axis.ticks.y = ggplot2::element_blank())

# --- C : CpG ----------------------------------------------------------------
# One legend for the pair, on the left half.
csize <- ggplot2::theme(axis.title = ggplot2::element_text(size = 7.4),
                        axis.text = ggplot2::element_text(size = 6.4),
                        strip.text = ggplot2::element_text(size = 7),
                        plot.margin = ggplot2::margin(1, 2, 1, 2))
pC1 <- strip_text(C$density) + ggplot2::labs(tag = "C") + csize +
  ggplot2::theme(legend.position = "top",
                 legend.text = ggplot2::element_text(size = 6.6),
                 legend.key.size = ggplot2::unit(0.30, "cm"),
                 legend.margin = ggplot2::margin(0, 0, 0, 0))
pC2 <- strip_text(C$trio) + csize +
  ggplot2::theme(legend.position = "none")

# guides = "collect" combined with `&` fails on a patchwork that carries an
# inset_element, so the keys stay attached to panel A. They sit under the map
# and read for panel B as well, which uses the same superfamily colours.
# Row 1 is exactly as tall as the map needs at this width, plus the keys.
# Anything more and the map floats; anything less and it is squeezed while the
# panel pads the sides with ocean.
col_a_in <- FIG_W_IN * COL_A_SHARE
row1_in  <- col_a_in / A$aspect + LEGEND_IN
FIG_H_IN <- row1_in + ROW2_IN

fig3 <- (pA | pB) / (pC1 | pC2) +
  patchwork::plot_layout(heights = c(row1_in, ROW2_IN),
                         widths = c(COL_A_SHARE, 1 - COL_A_SHARE))

message(sprintf("Map aspect %.2f -> row 1 = %.2f in, row 2 = %.2f in",
                A$aspect, row1_in, ROW2_IN))
ggplot2::ggsave(file.path(dir_fig, "Fig3.pdf"), fig3,
                width = FIG_W_IN, height = FIG_H_IN, device = grDevices::cairo_pdf)
ggplot2::ggsave(file.path(dir_fig, "Fig3.svg"), fig3,
                width = FIG_W_IN, height = FIG_H_IN, device = svglite::svglite)
message(sprintf("Fig3 written: %.2f x %.2f in (%.0f x %.0f mm)",
                FIG_W_IN, FIG_H_IN, FIG_W_IN * 25.4, FIG_H_IN * 25.4))

log_session(TAG)
message("=== 17_figure3_final.R finished ===")
