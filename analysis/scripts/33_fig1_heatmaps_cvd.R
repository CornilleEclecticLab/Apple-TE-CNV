# =============================================================================
# Script:  33_fig1_heatmaps_cvd.R
# Purpose: Figure 1E and 1F, the structural and residual TE heatmaps.
#          The annotation bands were rebuilt for red-green colour blindness at
#          the request of Genome Biology. The blue-white-red gradient of the
#          published figure is unchanged; only the three Class / Family /
#          Superfamily bands differ.
# Author:  Antho
# Outputs: Fig1_heatmap_STRUCTURE_CVD and Fig1_heatmap_RESIDUALS_CVD, PDF + SVG
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(stringr); library(pheatmap); library(grid)
})
# Directory holding the TE count matrix and the sample list; override with
# the TE_DATA_DIR environment variable.
setwd(Sys.getenv("TE_DATA_DIR", unset = file.path("..", "250923_TE-2")))
OUT <- Sys.getenv("FIG_OUT", unset = file.path("..", "..", "revision_finale", "figures_cvd"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---------- palette of the annotation bands ----------
## One colour per category, the same one in Figures 1, 2 and 3. Constraints:
##   - no green at all: under deuteranopia and protanopia green folds onto
##     yellow and merges with orange, whereas red keeps its own lightness;
##   - every "n.a." and "noCat" shares the same grey;
##   - Class I and Class II keep the red and blue of the published figure;
##   - each legend validated with CIEDE2000 on SIMULATED colours (normal
##     vision + deuteranopia + protanopia + tritanopia);
##   - lightness (L*) separation on top of hue: two colours of equal L* stay
##     hard to tell apart at a glance even at a high Delta E.
## Gradient: blue -> white -> red, as published. It diverges properly (poles
## 53 Delta E apart under all three kinds of colour blindness) and passes
## through no brown step.
GRADIENT <- colorRampPalette(c("blue", "white", "red"))(50)

## thickening factor for the three annotation bands (1 = stock pheatmap)
ANNOT_SCALE <- 2

TE_COLS <- c(
  "LTR"      = "#E69F00",  # Okabe-Ito orange    L* 71
  "LINE"     = "#56B4E9",  # sky blue            L* 70
  "TIR"      = "#332288",  # Paul Tol indigo     L* 22
  "grey"     = "#6E6E6E",  # mid grey            L* 46   noCat and every n.a.
  "other"    = "#000000",  # black               L*  0
  "Gypsy"    = "#CCE6FF",  # very pale blue      L* 90
  "Copia"    = "#A50F15",  # dark red            L* 35
  "Class I"  = "#FF0000",  # red                 L* 53   as published
  "Class II" = "#0000FF"   # blue                L* 32   as published
)

ann_colors <- list(
  Class = c("Class I"  = TE_COLS[["Class I"]],
            "Class II" = TE_COLS[["Class II"]],
            "n.a."     = TE_COLS[["grey"]]),
  Family = c("LTR" = TE_COLS[["LTR"]], "LINE" = TE_COLS[["LINE"]], "TIR" = TE_COLS[["TIR"]],
             "noCat" = TE_COLS[["grey"]], "other" = TE_COLS[["other"]]),
  Superfamily = c("Gypsy" = TE_COLS[["Gypsy"]],
                  "Copia" = TE_COLS[["Copia"]],
                  "n.a."  = TE_COLS[["grey"]])
)

## ---------- data ----------
DATA <- as.data.frame(fread("data_frame_countTE_cleaned.csv"))
names(DATA) <- str_replace(names(DATA), "LTR_PIF_Harbinger", "LTR")
rownames(DATA) <- DATA$V1; DATA$V1 <- NULL

GROUPS_structure <- fread("250923_sample_list_TE.csv") %>%
  dplyr::select(Sample_ID, cluster = group_structure_CORE_admix_K6) %>%
  filter(!grepl("Admix", cluster))
GROUPS_residuals <- fread("250923_sample_list_TE.csv") %>%
  group_by(Code_group_geo_residual) %>% filter(n() >= 3) %>% ungroup() %>%
  dplyr::select(Sample_ID, cluster = Code_group_geo_residual)
structure <- readLines("250923_TE_structure.txt")
residuals <- readLines("250923_TE_residuals.txt")

parse_TE_info <- function(te_names) {
  ## The six families with <=5 TEs (SINE, SSR, Helitron, pararetrovirus, MITE
  ## and mixture) amount to 12 TEs out of 2420 (0.5%) in the structural set and
  ## 2 out of 1360 in the residual set: one to four one-pixel columns, invisible
  ## whatever colour they are given. They used six of the thirteen colour slots
  ## and forced indistinguishable hues, so they are pooled into "other".
  data.frame(TE = te_names, stringsAsFactors = FALSE) %>%
    mutate(
      Class = case_when(str_detect(TE,"LTR|LINE|SINE") ~ "Class I",
                        str_detect(TE,"TIR|Helitron")  ~ "Class II",
                        TRUE ~ "n.a."),
      Family = case_when(str_detect(TE,"LTR")   ~ "LTR",
                         str_detect(TE,"LINE")  ~ "LINE",
                         str_detect(TE,"TIR")   ~ "TIR",
                         str_detect(TE,"noCat") ~ "noCat",
                         TRUE ~ "other"),
      Superfamily = case_when(str_detect(TE,"Gypsy") ~ "Gypsy",
                              str_detect(TE,"Copia") ~ "Copia",
                              TRUE ~ "n.a."))
}
make_group_means <- function(DATA, GROUPS, TE_list) {
  df <- DATA[, intersect(colnames(DATA), TE_list), drop = FALSE]
  ids <- intersect(rownames(df), GROUPS$Sample_ID)
  df <- df[ids, , drop = FALSE]; GROUPS <- GROUPS %>% filter(Sample_ID %in% ids)
  df <- log1p(df[match(GROUPS$Sample_ID, rownames(df)), , drop = FALSE])
  m <- df %>% as.data.frame() %>% mutate(Group = GROUPS$cluster) %>%
    group_by(Group) %>% summarise(across(everything(), mean), .groups="drop") %>% as.data.frame()
  rownames(m) <- m$Group; m$Group <- NULL; m
}

make_heatmap <- function(mat, output_name, w = 10, h = 5.2) {
  ai <- parse_TE_info(colnames(mat))
  ac <- ai %>% dplyr::select(Class, Family, Superfamily) %>% as.data.frame()
  rownames(ac) <- ai$TE
  ph <- pheatmap(as.matrix(mat), cluster_rows=TRUE, cluster_cols=TRUE,
    clustering_distance_rows="euclidean", clustering_distance_cols="correlation",
    clustering_method="ward.D2", scale="column", show_colnames=FALSE,
    color=GRADIENT, fontsize_row=10, annotation_col=ac, annotation_colors=ann_colors,
    border_color=NA, silent=TRUE)

  ## pheatmap pins the height of an annotation track at 8 bigpts. With 2420
  ## one-pixel columns that is too thin, so the "col_annotation" row is scaled
  ## up in the gtable. CAREFUL: pheatmap computes the "matrix" row height as
  ## 1npc minus the other rows, using the ORIGINAL value. Without subtracting
  ## the extra height from the matrix, the total exceeds 1npc and the whole
  ## figure shifts down, clipping the last row.
  ia <- which(ph$gtable$layout$name == "col_annotation")
  if (length(ia)) {
    r     <- unique(ph$gtable$layout$t[ia])
    mr    <- unique(ph$gtable$layout$t[ph$gtable$layout$name == "matrix"])
    extra <- ph$gtable$heights[r] * (ANNOT_SCALE - 1)
    ph$gtable$heights[r]  <- ph$gtable$heights[r] * ANNOT_SCALE
    ph$gtable$heights[mr] <- ph$gtable$heights[mr] - extra
  }

  for (dev in c("pdf","svg")) {
    f <- file.path(OUT, sprintf("Fig1_heatmap_%s_CVD.%s", output_name, dev))
    if (dev=="pdf") cairo_pdf(f, width=w, height=h) else svg(f, width=w, height=h)
    grid.newpage(); grid.draw(ph$gtable); dev.off()
    cat("  written:", f, "\n")
  }
  invisible(ph)
}

new_labels <- c(Q1="SYLW", Q2="DOMD", Q3="BACC", Q4="SIOR", Q5="SYLE", Q6="DOMC")
ms <- make_group_means(DATA, GROUPS_structure, structure)
rownames(ms) <- new_labels[rownames(ms)]
cat("STRUCTURE:", nrow(ms), "groups x", ncol(ms), "TEs\n"); make_heatmap(ms, "STRUCTURE")

mr <- make_group_means(DATA, GROUPS_residuals, residuals)
cat("RESIDUALS:", nrow(mr), "groups x", ncol(mr), "TEs\n"); make_heatmap(mr, "RESIDUALS")
