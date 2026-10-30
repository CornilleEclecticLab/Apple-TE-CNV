# =============================================================================
# Script:  36_figS8_te_pca.R
# Purpose: Supplementary Figure S8 - PCA on the copy number of the 157,816
#          non-null TE variants genotyped across the 189 Malus accessions,
#          shown in THREE DIMENSIONS, coloured by SNP-based population
#          assignment, admixed individuals as hollow circles.
#
#          S7 and S8 are the same PCA. S7 is the pair of 2D projections
#          coloured by geographical group; S8 is the 3D view coloured by SNP
#          assignment, mirroring the S5 / S6 pair on the SNP side.
#
#          NOTHING is recomputed. The coordinates are the published ones,
#          distributed as Table S5 of Additional file 2, and the axis
#          percentages are those printed on the published figure
#          (22.6 / 9 / 6 %, read from PCA_TE_PC1_vs_PC2.pdf). The colours are
#          the ones used for the SNP PCA in Figure S5, so that a reader can
#          move between the two figures.
#
# Author:  Antho
# Inputs:  data/PCA_TE_coordinates_published.csv   (= Table S5)
#          data/250923_sample_list_TE.csv
# Outputs: FigS08_PCA_TE_SNPgroups, SVG + PDF
# =============================================================================
suppressPackageStartupMessages({library(scatterplot3d); library(svglite)})

OUT <- Sys.getenv("FIG_OUT", unset = file.path("..", "revision_finale", "figures_cvd"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

GROUP <- c(Q1 = "SYLW", Q2 = "DOMD", Q3 = "BACC",
           Q4 = "SIOR", Q5 = "SYLE", Q6 = "DOMC")
COL <- c(BACC = "#6EB2E4", DOMC = "#C17DA5", DOMD = "#DCA327",
         SIOR = "#C66526", SYLE = "#EEE461", SYLW = "#3171AD")
ORDER <- c("BACC", "DOMC", "DOMD", "SIOR", "SYLE", "SYLW")

## Variance explained, as printed on the published figure.
VAR <- c(22.6, 9, 6)

pca  <- read.csv("data/PCA_TE_coordinates_published.csv", stringsAsFactors = FALSE)
meta <- read.csv("data/250923_sample_list_TE.csv", check.names = FALSE,
                 stringsAsFactors = FALSE)
names(meta)[1] <- "Sample_ID"              # the file carries a UTF-8 BOM
pca  <- merge(pca, meta[, c("Sample_ID", "group_structure_CORE_admix_K6")],
              by = "Sample_ID", all.x = TRUE)
stopifnot(nrow(pca) == 189, !any(is.na(pca$group_structure_CORE_admix_K6)))

q           <- sub("_Admix$", "", pca$group_structure_CORE_admix_K6)
pca$group   <- GROUP[q]
pca$admixed <- grepl("_Admix$", pca$group_structure_CORE_admix_K6)
pca$fg      <- COL[pca$group]
pca$bg      <- ifelse(pca$admixed, NA, pca$fg)

lab <- function(i) sprintf("PC%d (%s%%)", i, format(VAR[i], trim = TRUE))

draw <- function() {
  layout(matrix(1:2, ncol = 1), heights = c(3.5, 1))
  par(mar = c(2.2, 2.4, 2.0, 1.2), cex = 0.9)
  s3 <- scatterplot3d(pca$PC1, pca$PC2, pca$PC3,
                      xlab = lab(1), ylab = lab(2), zlab = lab(3),
                      main = "PCA TE - PC1, PC2 and PC3",
                      pch = "", box = TRUE, grid = TRUE,
                      angle = 55, scale.y = 0.85,
                      col.axis = "grey65", col.grid = "grey88",
                      cex.axis = 0.7, cex.lab = 0.85)
  p <- s3$xyz.convert(pca$PC1, pca$PC2, pca$PC3)
  points(p$x, p$y, pch = 21, cex = 1.05, lwd = 1.1,
         col = pca$fg, bg = pca$bg)

  lab_leg <- as.vector(rbind(ORDER, paste0(ORDER, "_Admix")))
  col_leg <- as.vector(rbind(COL[ORDER], COL[ORDER]))
  bg_leg  <- as.vector(rbind(COL[ORDER], rep(NA, length(ORDER))))
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", ncol = 4, bty = "n", cex = 0.78, pt.cex = 1.05,
         legend = lab_leg, pch = 21, pt.lwd = 1.1, col = col_leg, pt.bg = bg_leg)
}

W <- 170 / 25.4; H <- 5.9
f <- file.path(OUT, "FigS08_PCA_TE_3D")
svglite::svglite(paste0(f, ".svg"), width = W, height = H); draw(); dev.off()
grDevices::cairo_pdf(paste0(f, ".pdf"), width = W, height = H); draw(); dev.off()

cat(sprintf("Figure S8: %d accessions, %d admixed\n", nrow(pca), sum(pca$admixed)))
cat(sprintf("PC1 %.0f..%.0f | PC2 %.0f..%.0f | PC3 %.0f..%.0f\n",
            min(pca$PC1), max(pca$PC1), min(pca$PC2), max(pca$PC2),
            min(pca$PC3), max(pca$PC3)))
cat("  written:", f, ".pdf / .svg\n")
