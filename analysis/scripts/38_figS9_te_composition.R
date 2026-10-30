# =============================================================================
# Script:  38_figS9_te_composition.R
# Purpose: Supplementary Figure S9 - distribution of the 2,420 TE consensus
#          sequences identified as structural TEs, per TE category.
#
#          Redrawn with the colour-blind-safe palette of the paper, so that a
#          category keeps the same colour here as in Figures 1, 2 and 3.
#
#          NOTHING is recomputed: the categories are read off the TE names in
#          250923_TE_structure.txt, exactly as the published figure does, and
#          the counts are printed to the console so they can be checked against
#          it (noCat 667, LTR 652, Gypsy 555, LINE 209, Copia 168, TIR 157,
#          Helitron 4, SINE 3, SSR 2, MITE 1, mixture 1, pararetrovirus 1).
#
# Author:  Antho
# Outputs: FigS09_TE_composition_structural, SVG + PDF
# =============================================================================
suppressPackageStartupMessages({library(svglite)})

OUT  <- Sys.getenv("FIG_OUT", unset = file.path("..", "revision_finale", "figures_cvd"))
DATA <- Sys.getenv("TE_DATA_DIR", unset = file.path("..", "250923_TE-2"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## Palette of the paper. The first ten are the colours already used in
## Figures 1, 2 and 3; SINE and mixture are the two extra slots this figure
## needs. Twelve categories is past the point where every pair can be told
## apart under red-green colour blindness - the measured ceiling is a minimum
## Delta E of 12 - but here each bar is separate and labelled on the axis, so
## the colour only has to carry the link to the other figures.
TE_COL <- c(
  noCat          = "#6E6E6E",
  LTR            = "#E69F00",
  Gypsy          = "#CCE6FF",
  LINE           = "#56B4E9",
  Copia          = "#A50F15",
  TIR            = "#332288",
  Helitron       = "#0072B2",
  SINE           = "#EE4444",
  SSR            = "#FFD6A5",
  MITE           = "#1A1A1A",
  mixture        = "#4D4D4D",
  pararetrovirus = "#AAAAAA")

te <- readLines(file.path(DATA, "250923_TE_structure.txt"))

## Order matters: Gypsy and Copia are LTR subfamilies, so they must be matched
## before the generic LTR pattern, exactly as in the published figure.
classify <- function(x) {
  ifelse(grepl("Gypsy", x),          "Gypsy",
  ifelse(grepl("Copia", x),          "Copia",
  ifelse(grepl("LTR", x),            "LTR",
  ifelse(grepl("LINE", x),           "LINE",
  ifelse(grepl("SINE", x),           "SINE",
  ifelse(grepl("TIR", x),            "TIR",
  ifelse(grepl("Helitron", x),       "Helitron",
  ifelse(grepl("SSR", x),            "SSR",
  ifelse(grepl("MITE", x),           "MITE",
  ifelse(grepl("mixture", x),        "mixture",
  ifelse(grepl("pararetrovirus", x), "pararetrovirus",
  ifelse(grepl("noCat", x),          "noCat", NA_character_))))))))))))
}

cat_te <- classify(te)
stopifnot(!any(is.na(cat_te)))
cnt <- sort(table(cat_te), decreasing = TRUE)
stopifnot(all(names(cnt) %in% names(TE_COL)))

draw <- function() {
  par(mar = c(7.5, 4.6, 3.2, 1.2), cex = 0.9)
  b <- barplot(cnt, col = TE_COL[names(cnt)], border = "grey25",
               las = 2, ylim = c(0, max(cnt) * 1.12),
               ylab = "Number of sequences", xlab = "",
               main = sprintf("TE composition - structural TEs\nTotal elements: %d",
                              length(te)))
  text(b, as.numeric(cnt), labels = as.numeric(cnt), pos = 3, cex = 0.75)
  mtext("TE category", side = 1, line = 6.2, cex = 0.9)
}

W <- 170 / 25.4; H <- 5.4
f <- file.path(OUT, "FigS09_TE_composition_structural")
svglite::svglite(paste0(f, ".svg"), width = W, height = H); draw(); dev.off()
grDevices::cairo_pdf(paste0(f, ".pdf"), width = W, height = H); draw(); dev.off()

cat(sprintf("Figure S9: %d structural TEs in %d categories\n",
            length(te), length(cnt)))
print(cnt)
cat("  written:", f, ".pdf / .svg\n")
