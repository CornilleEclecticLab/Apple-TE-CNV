# =============================================================================
# Script:  39_figS7_te_pca_geo.R
# Purpose: Supplementary Figure S7 - PCA on the copy number of the 157,816
#          non-null TE variants genotyped across the 189 Malus accessions,
#          two projections (PC1 vs PC2 and PC1 vs PC3), individuals coloured
#          by taxonomic and geographic group.
#
#          S7 and S8 are the same PCA: S7 is this pair of 2D projections
#          coloured by geographical group, S8 the 3D view coloured by SNP
#          assignment. S6 and S7 share one colour and symbol scheme, as do
#          S5 and S8.
#
#          NOTHING is recomputed. The coordinates are the published ones
#          (Table S5 of Additional file 2) and the axis percentages those
#          printed on the published figure.
#
# Author:  Antho
# Outputs: FigS07_PCA_TE_geo, SVG + PDF
# =============================================================================
suppressPackageStartupMessages({library(svglite)})

OUT <- Sys.getenv("FIG_OUT", unset = file.path("..", "revision_finale", "figures_cvd"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## Colours and symbols measured in the published Figure S7, shared with S6.
## BAC_FSU and SIE_CN carry the same colour in the published figure; that is
## reproduced rather than silently corrected, and the symbols separate them.
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

## Variance explained, as printed on the published figure.
VAR <- c(22.6, 9, 6)

pca  <- read.csv("data/PCA_TE_coordinates_published.csv", stringsAsFactors = FALSE)
meta <- read.csv("data/250923_sample_list_TE.csv", check.names = FALSE,
                 stringsAsFactors = FALSE)
names(meta)[1] <- "Sample_ID"               # the file carries a UTF-8 BOM
pca  <- merge(pca, meta[, c("Sample_ID", "Code_group_geo_residual")],
              by = "Sample_ID", all.x = TRUE)
pca  <- pca[!is.na(pca$Code_group_geo_residual), ]
stopifnot(nrow(pca) == 189, all(pca$Code_group_geo_residual %in% GEO$group))

## Filled symbols get a thin dark edge: three of the sylvestris greens are
## nearly white and vanish on the white panel without one. The fill colours
## themselves are untouched.
i       <- match(pca$Code_group_geo_residual, GEO$group)
pca$pch <- GEO$pch[i]
pca$fg  <- ifelse(GEO$style[i] == "fill", "grey30", GEO$colour[i])
pca$bg  <- ifelse(GEO$style[i] == "fill", GEO$colour[i], NA)

lab <- function(k) sprintf("PC%d (%s%%)", k, format(VAR[k], trim = TRUE))

panel <- function(yk) {
  plot(pca$PC1, pca[[paste0("PC", yk)]], type = "n",
       xlab = lab(1), ylab = lab(yk),
       main = sprintf("PCA TE - PC1 vs PC%d", yk))
  points(pca$PC1, pca[[paste0("PC", yk)]], pch = pca$pch, cex = 1.05,
         lwd = 1.1, col = pca$fg, bg = pca$bg)
}

draw <- function() {
  layout(matrix(1:3, ncol = 1), heights = c(4, 4, 1.5))
  par(mar = c(4.2, 4.4, 2.2, 1.2), cex = 0.9)
  panel(2); panel(3)
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", ncol = 5, bty = "n", cex = 0.72, pt.cex = 1,
         title = "Geographical groups",
         legend = GEO$group, pch = GEO$pch, pt.lwd = 1.1,
         col = ifelse(GEO$style == "fill", "grey30", GEO$colour),
         pt.bg = ifelse(GEO$style == "fill", GEO$colour, NA))
}

W <- 170 / 25.4; H <- 9.8
f <- file.path(OUT, "FigS07_PCA_TE_geo")
svglite::svglite(paste0(f, ".svg"), width = W, height = H); draw(); dev.off()
grDevices::cairo_pdf(paste0(f, ".pdf"), width = W, height = H); draw(); dev.off()

cat(sprintf("Figure S7: %d accessions, %d geographical groups\n",
            nrow(pca), length(unique(pca$Code_group_geo_residual))))
cat(sprintf("axes: %s / %s / %s\n", lab(1), lab(2), lab(3)))
cat("  written:", f, ".pdf / .svg\n")
