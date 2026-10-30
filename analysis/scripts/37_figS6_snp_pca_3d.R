# =============================================================================
# Script:  37_figS6_snp_pca_3d.R
# Purpose: Supplementary Figure S6 - PCA on the 1,297 unlinked synonymous SNPs
#          genotyped in the 189 Malus accessions, shown in three dimensions,
#          individuals coloured by taxonomic and geographic group.
#
#          The published version was a 650 x 972 px bitmap, 97 dpi at 170 mm
#          wide, the lowest resolution of the whole supplement. This redraws it
#          as vector output.
#
#          NOTHING is recomputed: the coordinates come from
#          data/PCA_SNP_coordinates.csv (= Table S4) and the axis percentages
#          from outputs/01_snp_pca_axes/tables/pca_variance_explained.csv.
#
#          Colours and symbols are MEASURED in the published Figure S7
#          (sources_supplementaires/PCA_TE_PC1_vs_PC2.pdf), so that S6 and S7
#          share one scheme: filled triangles for M. baccata, hollow triangles
#          for the cultivated groups, filled circles for M. orientalis, filled
#          squares for M. sieversii, filled diamonds for M. sylvestris.
#
# Author:  Antho
# Outputs: FigS06_PCA_SNP_3D, SVG + PDF
# =============================================================================
suppressPackageStartupMessages({library(scatterplot3d); library(svglite)})

OUT <- Sys.getenv("FIG_OUT", unset = file.path("..", "revision_finale", "figures_cvd"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## Measured in PCA_TE_PC1_vs_PC2.pdf, legend "Geographical Groups".
## NOTE: BAC_FSU and SIE_CN carry the same colour in the published figure.
## That is reproduced here rather than silently corrected; BAC_FSU has a
## single accession, so the two are not confusable in practice.
GEO <- rbind(
  c("BAC_CN",      "#000000", "fill",   "24"),
  c("BAC_FSU",     "#A9746E", "fill",   "24"),
  c("DOMC",        "#FF0000", "hollow", "24"),
  c("DOMD",        "#8B0000", "hollow", "24"),
  c("DOMR",        "#F08080", "hollow", "24"),
  c("OR_AM",       "#7E57C2", "fill",   "21"),
  c("OR_TK",       "#512DA8", "fill",   "21"),
  c("SIE_CN",      "#A9746E", "fill",   "22"),
  c("SIE_KZ",      "#5D4037", "fill",   "22"),
  c("SYL_AT",      "#28B463", "fill",   "23"),
  c("SYL_DE",      "#2ECC71", "fill",   "23"),
  c("SYL_DK",      "#1D8348", "fill",   "23"),
  c("SYL_FR_EAST", "#ABEBC6", "fill",   "23"),
  c("SYL_FR_MID",  "#D5F5E3", "fill",   "23"),
  c("SYL_FR_WEST", "#EAFBF2", "fill",   "23"),
  c("SYL_IT",      "#82E0AA", "fill",   "23"),
  c("SYL_MKD",     "#58D68D", "fill",   "23"),
  c("SYL_NW",      "#3B6440", "fill",   "23"),
  c("SYL_RO",      "#239B56", "fill",   "23"))
colnames(GEO) <- c("group", "colour", "style", "pch")
GEO <- as.data.frame(GEO, stringsAsFactors = FALSE)
GEO$pch <- as.integer(GEO$pch)

pca  <- read.csv("data/PCA_SNP_coordinates.csv", stringsAsFactors = FALSE)
var  <- read.csv(file.path("outputs", "01_snp_pca_axes", "tables",
                           "pca_variance_explained.csv"))
meta <- read.csv("data/250923_sample_list_TE.csv", check.names = FALSE,
                 stringsAsFactors = FALSE)
names(meta)[1] <- "Sample_ID"               # the file carries a UTF-8 BOM
pca  <- merge(pca, meta[, c("Sample_ID", "Code_group_geo_residual")],
              by = "Sample_ID", all.x = TRUE)
pca  <- pca[!is.na(pca$Code_group_geo_residual), ]
stopifnot(all(pca$Code_group_geo_residual %in% GEO$group))

## Filled symbols get a thin dark edge. Three of the sylvestris greens are
## nearly white (#EAFBF2, #D5F5E3, #ABEBC6) and vanish on the white panel
## without one. The fill colours themselves are untouched.
i        <- match(pca$Code_group_geo_residual, GEO$group)
pca$pch  <- GEO$pch[i]
pca$fg   <- ifelse(GEO$style[i] == "fill", "grey30", GEO$colour[i])
pca$bg   <- ifelse(GEO$style[i] == "fill", GEO$colour[i], NA)

lab <- function(k) sprintf("PC%d (%s%%)", k,
                           format(round(var$variance_explained[k], 1), trim = TRUE))

draw <- function() {
  layout(matrix(1:2, ncol = 1), heights = c(3.5, 1))
  par(mar = c(2.2, 2.4, 2.0, 1.2), cex = 0.9)
  s3 <- scatterplot3d(pca$PC1, pca$PC2, pca$PC3,
                      xlab = lab(1), ylab = lab(2), zlab = lab(3),
                      main = "PCA SNP - PC1, PC2 and PC3",
                      pch = "", box = TRUE, grid = TRUE,
                      angle = 55, scale.y = 0.85,
                      col.axis = "grey65", col.grid = "grey88",
                      cex.axis = 0.7, cex.lab = 0.85)
  p <- s3$xyz.convert(pca$PC1, pca$PC2, pca$PC3)
  points(p$x, p$y, pch = pca$pch, cex = 1.05, lwd = 1.1,
         col = pca$fg, bg = pca$bg)

  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", ncol = 5, bty = "n", cex = 0.72, pt.cex = 1,
         title = "Geographical groups",
         legend = GEO$group, pch = GEO$pch, pt.lwd = 1.1,
         col = ifelse(GEO$style == "fill", "grey30", GEO$colour),
         pt.bg = ifelse(GEO$style == "fill", GEO$colour, NA))
}

W <- 170 / 25.4; H <- 5.9
f <- file.path(OUT, "FigS06_PCA_SNP_3D")
svglite::svglite(paste0(f, ".svg"), width = W, height = H); draw(); dev.off()
grDevices::cairo_pdf(paste0(f, ".pdf"), width = W, height = H); draw(); dev.off()

cat(sprintf("Figure S6: %d accessions, %d geographical groups\n",
            nrow(pca), length(unique(pca$Code_group_geo_residual))))
cat(sprintf("axes: %s / %s / %s\n", lab(1), lab(2), lab(3)))
cat("  written:", f, ".pdf / .svg\n")
