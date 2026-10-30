# =============================================================================
# Script:  35_figS5_snp_pca.R
# Purpose: Supplementary Figure S5 - PCA on the 1,297 unlinked synonymous SNPs
#          genotyped in the 189 Malus accessions, first two components.
#
#          The published version was embedded in the supplementary Word file as
#          a 909 x 910 px bitmap, i.e. 136 dpi at 170 mm wide, which the editor
#          asked us to improve. This script redraws it as true vector output,
#          so the question of resolution does not arise at all.
#
#          NOTHING is recomputed: the coordinates come from the same table that
#          produced the published figure (data/PCA_SNP_coordinates.csv, also
#          distributed as Table S4), and the axis percentages come from
#          outputs/01_snp_pca_axes/tables/pca_variance_explained.csv. Colours
#          and symbols are those measured in the published PDF.
#
# Author:  Antho
# Inputs:  data/PCA_SNP_coordinates.csv
#          outputs/01_snp_pca_axes/tables/pca_variance_explained.csv
# Outputs: FigS05_PCA_SNP_PC1_PC2, SVG + PDF
# =============================================================================
suppressPackageStartupMessages({library(svglite)})

OUT <- Sys.getenv("FIG_OUT", unset = file.path("..", "revision_finale", "figures_cvd"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## Q cluster -> group name, the same mapping as everywhere else in the project.
GROUP <- c(Q1 = "SYLW", Q2 = "DOMD", Q3 = "BACC",
           Q4 = "SIOR", Q5 = "SYLE", Q6 = "DOMC")

## Colours measured in the published PDF (sources_supplementaires/
## PCA_SNP_PC1_vs_PC2.pdf), not guessed from the R colour names.
COL <- c(BACC = "#6EB2E4", DOMC = "#C17DA5", DOMD = "#DCA327",
         SIOR = "#C66526", SYLE = "#EEE461", SYLW = "#3171AD")
ORDER <- c("BACC", "DOMC", "DOMD", "SIOR", "SYLE", "SYLW")

pca <- read.csv("data/PCA_SNP_coordinates.csv", stringsAsFactors = FALSE)
var <- read.csv(file.path("outputs", "01_snp_pca_axes", "tables",
                          "pca_variance_explained.csv"))

q            <- sub("_Admix$", "", pca$group_structure_CORE_admix_K6)
pca$group    <- GROUP[q]
pca$admixed  <- grepl("_Admix$", pca$group_structure_CORE_admix_K6)
stopifnot(!any(is.na(pca$group)))

## Admixed individuals are drawn as hollow circles, as in the published figure.
pt_col <- COL[pca$group]
pca$bg  <- ifelse(pca$admixed, NA, pt_col)
pca$fg  <- pt_col

lab <- function(i) sprintf("PC%d (%s%%)", i,
                           format(round(var$variance_explained[i], 1), trim = TRUE))

draw <- function() {
  layout(matrix(1:2, ncol = 1), heights = c(4.2, 1))
  par(mar = c(4.2, 4.4, 2.2, 1.2), cex = 0.9)
  plot(pca$PC1, pca$PC2, type = "n", xlab = lab(1), ylab = lab(2),
       main = "PCA SNP - PC1 vs PC2")
  points(pca$PC1, pca$PC2, pch = 21, cex = 1.15, lwd = 1.1,
         col = pca$fg, bg = pca$bg)

  ## Each group and its admixed counterpart sit one above the other, three
  ## columns of four, exactly as in the published legend. legend() fills
  ## column by column, so the entries are interleaved here.
  lab_leg <- as.vector(rbind(ORDER, paste0(ORDER, "_Admix")))
  col_leg <- as.vector(rbind(COL[ORDER], COL[ORDER]))
  bg_leg  <- as.vector(rbind(COL[ORDER], rep(NA, length(ORDER))))
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", ncol = 3, bty = "n", cex = 0.8, pt.cex = 1.1,
         legend = lab_leg, pch = 21, pt.lwd = 1.1,
         col = col_leg, pt.bg = bg_leg)
}

W <- 170 / 25.4; H <- 6.4
f <- file.path(OUT, "FigS05_PCA_SNP_PC1_PC2")
svglite::svglite(paste0(f, ".svg"), width = W, height = H); draw(); dev.off()
grDevices::cairo_pdf(paste0(f, ".pdf"), width = W, height = H); draw(); dev.off()

cat(sprintf("Figure S5 redrawn: %d accessions, %d admixed\n",
            nrow(pca), sum(pca$admixed)))
print(table(pca$group, ifelse(pca$admixed, "admixed", "core")))
cat(sprintf("axes: %s / %s\n", lab(1), lab(2)))
cat("  written:", f, ".pdf / .svg\n")
