# =============================================================================
# Script:  36_figS8_te_pca.R
# Purpose: Supplementary Figure S8 - PCA on the copy number of the 157,816
#          non-null TE variants genotyped across the 189 Malus accessions,
#          coloured by SNP-based population assignment, admixed individuals as
#          hollow circles. Two projections, PC1 vs PC2 and PC1 vs PC3.
#
#          S7 and S8 are the SAME PCA; only the colouring differs. S7 colours
#          by geographical group, S8 by SNP-based assignment.
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
suppressPackageStartupMessages({library(svglite)})

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

panel <- function(yi) {
  plot(pca$PC1, pca[[paste0("PC", yi)]], type = "n",
       xlab = lab(1), ylab = lab(yi),
       main = sprintf("PCA TE - PC1 vs PC%d", yi))
  points(pca$PC1, pca[[paste0("PC", yi)]], pch = 21, cex = 1.15, lwd = 1.1,
         col = pca$fg, bg = pca$bg)
}

draw <- function() {
  layout(matrix(1:3, ncol = 1), heights = c(4, 4, 1.1))
  par(mar = c(4.2, 4.4, 2.2, 1.2), cex = 0.9)
  panel(2); panel(3)
  lab_leg <- as.vector(rbind(ORDER, paste0(ORDER, "_Admix")))
  col_leg <- as.vector(rbind(COL[ORDER], COL[ORDER]))
  bg_leg  <- as.vector(rbind(COL[ORDER], rep(NA, length(ORDER))))
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", ncol = 3, bty = "n", cex = 0.8, pt.cex = 1.1,
         legend = lab_leg, pch = 21, pt.lwd = 1.1, col = col_leg, pt.bg = bg_leg)
}

W <- 170 / 25.4; H <- 9.4
f <- file.path(OUT, "FigS08_PCA_TE_SNPgroups")
svglite::svglite(paste0(f, ".svg"), width = W, height = H); draw(); dev.off()
grDevices::cairo_pdf(paste0(f, ".pdf"), width = W, height = H); draw(); dev.off()

cat(sprintf("Figure S8: %d accessions, %d admixed\n", nrow(pca), sum(pca$admixed)))
cat(sprintf("PC1 %.0f..%.0f | PC2 %.0f..%.0f | PC3 %.0f..%.0f\n",
            min(pca$PC1), max(pca$PC1), min(pca$PC2), max(pca$PC2),
            min(pca$PC3), max(pca$PC3)))
cat("  written:", f, ".pdf / .svg\n")
