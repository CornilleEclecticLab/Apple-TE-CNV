# =============================================================================
# Script:  30_pca_combat_sensitivity.R
# Purpose: The Methods now state that batch and species are strongly confounded
#          and that ComBat-seq was run without a biological covariate. A reviewer
#          will ask what that correction does to the structure the whole paper
#          rests on. This recomputes the TE PCA with and without the correction
#          and compares them to each other and to the published coordinates.
# Author:  Antho
# =============================================================================
source("scripts/00_setup.R")
suppressMessages({library(readr); library(dplyr); library(tibble)})
TAG <- "30_pca_combat_sensitivity"
dir_tab <- file.path(dir_outputs, TAG, "tables"); dir_fig <- file.path(dir_outputs, TAG, "figures")
for (d in c(dir_tab, dir_fig)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

M <- as.matrix(read.csv(file.path(dir_data, "TE_counts_clean_wide.csv"),
                        row.names = 1, check.names = FALSE))
message(sprintf("matrice brute : %d accessions x %d TEs", nrow(M), ncol(M)))

meta <- read_csv(file.path(dir_data, "250923_sample_list_TE.csv"), show_col_types = FALSE)
meta <- meta[match(rownames(M), meta$Sample_ID), ]
batch <- as.factor(meta$Dataset_source)
message("batch : "); print(table(batch, meta$Species))

# --- 1. uncorrected: library-size normalisation, then PCA --------------------
cpm <- sweep(M, 1, rowSums(M), "/") * 1e6
keep <- apply(cpm, 2, function(x) sd(x) > 0)
message(sprintf("TEs with non-zero variance: %d", sum(keep)))
pca_raw <- prcomp(cpm[, keep], center = TRUE, scale. = TRUE)

# --- 2. with ComBat-seq, exactly as in the original pipeline -----------------
suppressMessages(library(sva))
adj <- sva::ComBat_seq(t(round(M)), batch = as.numeric(batch), group = NULL)
cpm_adj <- sweep(t(adj), 1, rowSums(t(adj)), "/") * 1e6
keep2 <- apply(cpm_adj, 2, function(x) sd(x) > 0)
pca_adj <- prcomp(cpm_adj[, keep2], center = TRUE, scale. = TRUE)

ve <- function(p, k = 5) round(100 * (p$sdev^2 / sum(p$sdev^2))[1:k], 2)
message("\nvariance explained (%):")
message("  sans ComBat : ", paste(ve(pca_raw), collapse = "  "))
message("  avec ComBat : ", paste(ve(pca_adj), collapse = "  "))

# --- 3. do the two PCAs tell the same story? ---------------------------------
cmp <- tibble(axis = paste0("PC", 1:5),
              pearson_abs = sapply(1:5, function(k)
                abs(cor(pca_raw$x[, k], pca_adj$x[, k]))),
              spearman_abs = sapply(1:5, function(k)
                abs(cor(pca_raw$x[, k], pca_adj$x[, k], method = "spearman"))))
print(as.data.frame(cmp)); write_csv(cmp, file.path(dir_tab, "pc_correlation.csv"))

# --- 4. which one matches the published coordinates (Table S5)? --------------
pub <- read_csv(file.path(dir_tab, "published_S5.csv"), show_col_types = FALSE)
pub <- pub[match(rownames(M), pub$Sample_ID), ]
for (nm in c("sans ComBat", "avec ComBat")) {
  X <- if (nm == "sans ComBat") pca_raw$x else pca_adj$x
  r <- sapply(1:3, function(k) abs(cor(X[, k], pub[[paste0("PC", k)]], use = "complete.obs")))
  message(sprintf("  %s vs Table S5 : PC1 %.3f  PC2 %.3f  PC3 %.3f", nm, r[1], r[2], r[3]))
}

# --- 5. does the separation between species survive? -------------------------
sep <- function(p, lbl) {
  d <- data.frame(PC1 = p$x[,1], PC2 = p$x[,2], sp = meta$Species)
  a1 <- summary(aov(PC1 ~ sp, d))[[1]]; a2 <- summary(aov(PC2 ~ sp, d))[[1]]
  tibble(model = lbl,
         R2_PC1 = a1[1,2] / sum(a1[,2]), p_PC1 = a1[1,5],
         R2_PC2 = a2[1,2] / sum(a2[,2]), p_PC2 = a2[1,5])
}
S <- bind_rows(sep(pca_raw, "sans ComBat"), sep(pca_adj, "avec ComBat"))
print(as.data.frame(S)); write_csv(S, file.path(dir_tab, "species_separation.csv"))

write_csv(as_tibble(pca_raw$x[,1:5], rownames = "Sample_ID"),
          file.path(dir_tab, "pca_scores_raw.csv"))
write_csv(as_tibble(pca_adj$x[,1:5], rownames = "Sample_ID"),
          file.path(dir_tab, "pca_scores_combat.csv"))
log_session(TAG)

# --- 6. do PC1 and PC2 rotate, or does the structure disappear? --------------
# When two axes have similar variances their order and orientation are
# unstable, so comparing PC1 with PC1 is meaningless. We compare the sub-
# espaces.
cc12 <- cancor(pca_raw$x[,1:2], pca_adj$x[,1:2])$cor
cc35 <- cancor(pca_raw$x[,1:3], pca_adj$x[,1:3])$cor
message(sprintf("\ncorrelations canoniques, plan PC1-PC2 : %s",
                paste(round(cc12,4), collapse = "  ")))
message(sprintf("correlations canoniques, PC1-PC3      : %s",
                paste(round(cc35,4), collapse = "  ")))

suppressMessages(library(vegan))
pr <- vegan::protest(pca_raw$x[,1:5], pca_adj$x[,1:5], permutations = 999)
message(sprintf("Procruste sur PC1-PC5 : correlation = %.4f, p = %.3g",
                pr$scale * 0 + pr$t0, pr$signif))

# the same question for the separation between species, in the plane rather
# than axis by axis
d1 <- data.frame(pca_raw$x[,1:2], sp = meta$Species)
d2 <- data.frame(pca_adj$x[,1:2], sp = meta$Species)
m1 <- vegan::adonis2(dist(d1[,1:2]) ~ sp, data = d1, permutations = 999)
m2 <- vegan::adonis2(dist(d2[,1:2]) ~ sp, data = d2, permutations = 999)
message(sprintf("PERMANOVA espece dans le plan PC1-PC2 : sans ComBat R2 = %.3f (p = %.3f) | avec ComBat R2 = %.3f (p = %.3f)",
                m1$R2[1], m1$`Pr(>F)`[1], m2$R2[1], m2$`Pr(>F)`[1]))

# total variance of the scores explained by species in the PC1-PC2 plane
tot <- function(p) {
  X <- p$x[,1:2]
  ss <- sum(sapply(1:2, function(k) { a <- summary(aov(X[,k] ~ meta$Species))[[1]]; a[1,2] }))
  sst <- sum(sapply(1:2, function(k) { a <- summary(aov(X[,k] ~ meta$Species))[[1]]; sum(a[,2]) }))
  ss/sst }
message(sprintf("part de variance PC1-PC2 expliquee par l espece : sans ComBat %.3f | avec ComBat %.3f",
                tot(pca_raw), tot(pca_adj)))

# --- 7. figure: the two PCAs side by side ------------------------------------
suppressMessages({library(ggplot2); library(patchwork)})
SP_COL <- c("M. baccata"    = "#4C72B0", "M. domestica" = "#C44E52",
            "M. orientalis" = "#8172B2", "M. sieversii" = "#CCB974",
            "M. sylvestris" = "#55A868")
mk <- function(p, lbl, tag) {
  v <- 100 * p$sdev^2 / sum(p$sdev^2)
  d <- data.frame(PC1 = p$x[,1], PC2 = p$x[,2],
                  Species = meta$Species, batch = meta$Dataset_source)
  ggplot(d, aes(PC1, PC2, fill = Species, shape = batch)) +
    geom_point(size = 2.1, colour = "grey25", stroke = .25, alpha = .9) +
    scale_fill_manual(values = SP_COL) + scale_shape_manual(values = c(21, 24)) +
    labs(title = lbl, tag = tag,
         x = sprintf("PC1 (%.1f%%)", v[1]), y = sprintf("PC2 (%.1f%%)", v[2])) +
    guides(fill = guide_legend(override.aes = list(shape = 21), order = 1)) +
    theme_bw(base_size = 9) +
    theme(legend.position = "right", legend.key.size = unit(.35, "cm"),
          legend.text = element_text(size = 7, face = "italic"),
          legend.title = element_text(size = 7.5),
          plot.title = element_text(size = 9, face = "bold"))
}
pA <- mk(pca_raw, "Without ComBat-seq", "A") + theme(legend.position = "none")
pB <- mk(pca_adj, "With ComBat-seq", "B") +
      theme(legend.position = "bottom", legend.box = "horizontal")
fig <- pA / pB + plot_layout(heights = c(1, 1.28))
ggsave(file.path(dir_fig, "FigS24_PCA_combat_sensitivity.pdf"), fig,
       width = 4.6, height = 7.2, device = cairo_pdf)
ggsave(file.path(dir_fig, "FigS24_PCA_combat_sensitivity.svg"), fig,
       width = 4.6, height = 7.2)
message("figure ecrite (PDF + SVG)")
