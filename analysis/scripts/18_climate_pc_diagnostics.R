# =============================================================================
# Script:  18_climate_pc_diagnostics.R
# Purpose: Answer the questions raised on the climate analysis after Reviewer 2:
#
#          1. candidates from the original subset/replication approach
#          2. candidates from the climate + PC1-3 model
#          3. their overlap
#          4. those specific to each approach
#          5. how many survive an FDR correction, for information
#
#          and, most importantly, WHY the count drops from 34 to 6:
#          is it a loss of power (climate effects unchanged, standard errors
#          inflated by collinearity with the PCs), or genuine confounding
#          (the climate effects themselves shrink once structure is modelled)?
#
#          Both models are fitted on EXACTLY the same rows, so the two sets of
#          coefficients are comparable. The previous script fitted them
#          separately and stored no standard errors, which is what made the
#          question unanswerable.
#
# Author:  Antho
# Project: TE_Project / final_pipeline (Genome Biology revision)
# =============================================================================

source("scripts/00_setup.R")
library(dplyr); library(tidyr); library(purrr); library(readr); library(tibble)

TAG      <- "18_climate_pc_diagnostics"
dir_out  <- file.path(dir_outputs, TAG)
dir_tab  <- file.path(dir_out, "tables")
dir_fig  <- file.path(dir_out, "figures")
for (d in c(dir_out, dir_tab, dir_fig)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

PC_DEPTH <- 3
GROUPS   <- c("SYLW", "SYLE", "SYL_Admix")
THR      <- 0.01
MIN_N    <- 8
GROUP_OF <- function(x) dplyr::case_when(
  x == "Q1" ~ "SYLW", x == "Q5" ~ "SYLE",
  x %in% c("Q1_Admix", "Q5_Admix") ~ "SYL_Admix", TRUE ~ NA_character_)

# --- Data (same inputs as 02_candidate_TEs.R) --------------------------------
counts <- read_csv(file.path(dir_data, "coord_syl_TE.csv"), show_col_types = FALSE) %>%
  mutate(Sample_ID = toupper(Sample_ID))
te_vars <- setdiff(names(counts), "Sample_ID")

sample_info <- read_csv(file.path(dir_data, "250923_sample_list_TE.csv"),
                        show_col_types = FALSE) %>%
  filter(Species == "M. sylvestris", !is.na(latitude), !is.na(longitude)) %>%
  transmute(Sample_ID = toupper(Sample_ID), latitude, longitude,
            group_label = GROUP_OF(group_structure_CORE_admix_K6)) %>%
  filter(!is.na(group_label))

pc_tbl <- read_csv(file.path(dir_outputs, "01_snp_pca_axes", "tables",
                             "PCA_SNP_coordinates_PC1_10.csv"), show_col_types = FALSE) %>%
  select(Sample_ID = sample_id, all_of(sprintf("PC%d", 1:PC_DEPTH))) %>%
  mutate(Sample_ID = toupper(Sample_ID))

clim_dir  <- file.path(dir_data, "worldclim_env_data", "climate", "wc2.1_10m")
bio_stack <- terra::rast(file.path(clim_dir, sprintf("wc2.1_10m_bio_%d.tif", 1:19)))
names(bio_stack) <- sprintf("bio%d", 1:19)

meta <- inner_join(sample_info, pc_tbl, by = "Sample_ID")
meta <- bind_cols(meta, as_tibble(terra::extract(bio_stack,
          as.matrix(meta[, c("longitude", "latitude")]))))
design <- meta %>% filter(Sample_ID %in% counts$Sample_ID) %>% arrange(Sample_ID)
count_mat <- as.matrix(counts[match(design$Sample_ID, counts$Sample_ID), te_vars, drop = FALSE])
message(sprintf("Design: %d individuals, %d residual TEs", nrow(design), length(te_vars)))
print(table(design$group_label))

# =============================================================================
# 1. Collinearity between the bioclimatic variables and the neutral PCs
# =============================================================================
# If climate and structure are strongly collinear, adding the PCs must inflate
# the standard error of the climate coefficient even when the effect is real.
collin <- map_dfr(GROUPS, function(g) {
  d <- design[design$group_label == g, , drop = FALSE]
  map_dfr(sprintf("bio%d", 1:19), function(bv) {
    ok <- complete.cases(d[, c(bv, sprintf("PC%d", 1:PC_DEPTH))])
    if (sum(ok) < MIN_N) return(tibble())
    fit <- lm(as.formula(sprintf("`%s` ~ PC1 + PC2 + PC3", bv)), data = d[ok, ])
    r2  <- summary(fit)$r.squared
    tibble(group = g, variable = bv, n = sum(ok), r2_bio_on_PC = r2,
           vif = 1 / (1 - r2))
  })
})
write_csv(collin, file.path(dir_tab, "collinearity_bio_vs_PC.csv"))
message("\nCollinearity of each bioclimatic variable with PC1-3 (R^2):")
print(as.data.frame(collin %>% group_by(group) %>%
  summarise(median_r2 = median(r2_bio_on_PC), max_r2 = max(r2_bio_on_PC),
            median_vif = median(vif), .groups = "drop")))

# =============================================================================
# 2. Both models on identical rows, keeping standard errors
# =============================================================================
fit_both <- function(grp) {
  d  <- design[design$group_label == grp, , drop = FALSE]
  cm <- count_mat[design$group_label == grp, , drop = FALSE]
  map_dfr(sprintf("bio%d", 1:19), function(bv) {
    preds <- c(bv, sprintf("PC%d", 1:PC_DEPTH))
    ok <- complete.cases(d[, preds, drop = FALSE])
    md <- d[ok, , drop = FALSE]; cmv <- cm[ok, , drop = FALSE]
    if (nrow(md) < MIN_N) return(tibble())
    f0 <- as.formula(sprintf("te_count ~ `%s`", bv))
    f1 <- as.formula(paste("te_count ~",
            paste(sprintf("`%s`", preds), collapse = " + ")))
    n <- length(te_vars)
    b0 <- s0 <- p0 <- b1 <- s1 <- p1 <- rep(NA_real_, n)
    for (i in seq_len(n)) {
      y <- as.numeric(cmv[, te_vars[i]])
      if (length(unique(y)) < 2) next
      dat <- md; dat$te_count <- y
      m0 <- suppressWarnings(tryCatch(glm(f0, data = dat, family = quasipoisson()),
                                      error = function(e) NULL))
      m1 <- suppressWarnings(tryCatch(glm(f1, data = dat, family = quasipoisson()),
                                      error = function(e) NULL))
      if (!is.null(m0)) { c0 <- summary(m0)$coefficients
        if (nrow(c0) >= 2) { b0[i] <- c0[2,1]; s0[i] <- c0[2,2]; p0[i] <- c0[2,4] } }
      if (!is.null(m1)) { c1 <- summary(m1)$coefficients
        if (nrow(c1) >= 2) { b1[i] <- c1[2,1]; s1[i] <- c1[2,2]; p1[i] <- c1[2,4] } }
    }
    tibble(group = grp, variable = bv, TE_ID = te_vars, n = nrow(md),
           beta_A = b0, se_A = s0, p_A = p0,
           beta_PC = b1, se_PC = s1, p_PC = p1)
  })
}
message("\nFitting both models on identical rows ...")
t0 <- Sys.time()
fits <- map_dfr(GROUPS, function(g) { message("  ", g, " ..."); fit_both(g) }) %>%
  filter(!is.na(p_A) | !is.na(p_PC))
message(sprintf("  %d TE x variable x group fits in %.1f min", nrow(fits),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
write_csv(fits, file.path(dir_tab, "glm_both_models_with_SE.csv"))

# FDR within each group x variable, for both models
fits <- fits %>% group_by(group, variable) %>%
  mutate(fdr_A = p.adjust(p_A, "BH"), fdr_PC = p.adjust(p_PC, "BH")) %>% ungroup()

# =============================================================================
# 3. The five sets the supervisor asked for
# =============================================================================
replicated <- function(pcol, thr = THR) {
  per <- lapply(GROUPS, function(g) {
    s <- fits[fits$group == g & !is.na(fits[[pcol]]) & fits[[pcol]] < thr, ]
    unique(s$TE_ID)
  })
  names(per) <- GROUPS
  list(per_group = per, shared = sort(Reduce(intersect, per)))
}
setA   <- replicated("p_A")
setPC  <- replicated("p_PC")
setA_f <- replicated("fdr_A", 0.05)
setPC_f<- replicated("fdr_PC", 0.05)

both     <- sort(intersect(setA$shared, setPC$shared))
onlyA    <- sort(setdiff(setA$shared, setPC$shared))
onlyPC   <- sort(setdiff(setPC$shared, setA$shared))

sets <- tibble(
  set = c("Model A (subset replication, p<0.01)",
          "Model A-PC (climate + PC1-3, p<0.01)",
          "overlap",
          "Model A only",
          "Model A-PC only",
          "union",
          "Model A, FDR<0.05 within group x variable",
          "Model A-PC, FDR<0.05 within group x variable"),
  n = c(length(setA$shared), length(setPC$shared), length(both),
        length(onlyA), length(onlyPC),
        length(union(setA$shared, setPC$shared)),
        length(setA_f$shared), length(setPC_f$shared)))
write_csv(sets, file.path(dir_tab, "candidate_sets.csv"))
message("\nCandidate sets:"); print(as.data.frame(sets))

write_csv(tibble(TE_ID = union(setA$shared, setPC$shared)) %>%
  mutate(in_modelA = TE_ID %in% setA$shared,
         in_modelA_PC = TE_ID %in% setPC$shared,
         support = case_when(in_modelA & in_modelA_PC ~ "both",
                             in_modelA_PC ~ "Model A-PC only",
                             TRUE ~ "Model A only"),
         passes_FDR_modelA = TE_ID %in% setA_f$shared,
         passes_FDR_modelA_PC = TE_ID %in% setPC_f$shared),
  file.path(dir_tab, "candidate_membership.csv"))

# =============================================================================
# 4. Power loss or confounding?
# =============================================================================
# Restricted to the TE x variable x group cells that drove the original
# candidate list: those significant under Model A at p < 0.01.
drivers <- fits %>%
  filter(TE_ID %in% setA$shared, !is.na(p_A), p_A < THR,
         !is.na(beta_PC), !is.na(se_A), se_A > 0)

diag <- drivers %>%
  mutate(abs_ratio = abs(beta_PC) / abs(beta_A),
         se_ratio  = se_PC / se_A,
         sign_kept = sign(beta_PC) == sign(beta_A))

summ <- tibble(
  n_cells                = nrow(diag),
  median_abs_beta_ratio  = median(diag$abs_ratio, na.rm = TRUE),
  q25_abs_beta_ratio     = quantile(diag$abs_ratio, .25, na.rm = TRUE),
  q75_abs_beta_ratio     = quantile(diag$abs_ratio, .75, na.rm = TRUE),
  pct_sign_preserved     = 100 * mean(diag$sign_kept, na.rm = TRUE),
  median_se_ratio        = median(diag$se_ratio, na.rm = TRUE),
  q25_se_ratio           = quantile(diag$se_ratio, .25, na.rm = TRUE),
  q75_se_ratio           = quantile(diag$se_ratio, .75, na.rm = TRUE))
write_csv(summ, file.path(dir_tab, "power_vs_confounding.csv"))
write_csv(diag, file.path(dir_tab, "power_vs_confounding_cells.csv"))
message("\nEffect size and standard error, Model A-PC relative to Model A,")
message("on the cells that produced the original candidate list:")
print(as.data.frame(summ))

log_session(TAG)
message("=== 18_climate_pc_diagnostics.R finished ===")

# =============================================================================
# 5. The model the reviewer actually proposed: PCs INSTEAD OF subsetting
# =============================================================================
# R2.5 says the subset design is an unclear control and that "including the
# neutral PCs as covariates may be a more appropriate approach". We had run
# PCs *within* each subset and still required replication across all three,
# which is the most conservative combination and not what was asked. Here the
# model is fitted once on all individuals, with the PCs carrying the structure.
pooled <- local({
  preds <- c(sprintf("PC%d", 1:PC_DEPTH))
  map_dfr(sprintf("bio%d", 1:19), function(bv) {
    ok <- complete.cases(design[, c(bv, preds), drop = FALSE])
    md <- design[ok, , drop = FALSE]; cmv <- count_mat[ok, , drop = FALSE]
    fml <- as.formula(paste("te_count ~",
             paste(sprintf("`%s`", c(bv, preds)), collapse = " + ")))
    n <- length(te_vars); b <- s <- p <- rep(NA_real_, n)
    for (i in seq_len(n)) {
      y <- as.numeric(cmv[, te_vars[i]])
      if (length(unique(y)) < 2) next
      dat <- md; dat$te_count <- y
      m <- suppressWarnings(tryCatch(glm(fml, data = dat, family = quasipoisson()),
                                     error = function(e) NULL))
      if (is.null(m)) next
      cf <- summary(m)$coefficients
      if (nrow(cf) >= 2) { b[i] <- cf[2,1]; s[i] <- cf[2,2]; p[i] <- cf[2,4] }
    }
    tibble(variable = bv, TE_ID = te_vars, n = nrow(md),
           beta = b, se = s, p = p)
  })
}) %>% filter(!is.na(p)) %>%
  group_by(variable) %>% mutate(fdr_within = p.adjust(p, "BH")) %>% ungroup() %>%
  mutate(fdr_global = p.adjust(p, "BH"))
write_csv(pooled, file.path(dir_tab, "glm_pooled_PC.csv"))

pooled_sets <- tibble(
  criterion = c("p < 0.01, any bioclimatic variable",
                "FDR < 0.05 within variable",
                "FDR < 0.05 across all TE x variable tests"),
  n_TE = c(n_distinct(pooled$TE_ID[pooled$p < 0.01]),
           n_distinct(pooled$TE_ID[pooled$fdr_within < 0.05]),
           n_distinct(pooled$TE_ID[pooled$fdr_global < 0.05])))
write_csv(pooled_sets, file.path(dir_tab, "pooled_PC_sets.csv"))
message("\nPooled model on all individuals, climate + PC1-3:")
print(as.data.frame(pooled_sets))

pooled_ids <- unique(pooled$TE_ID[pooled$fdr_within < 0.05])
message(sprintf("\nOverlap of the FDR-significant pooled set with the 34: %d",
                length(intersect(pooled_ids, setA$shared))))
write_csv(tibble(TE_ID = pooled_ids,
                 also_in_modelA_34 = pooled_ids %in% setA$shared),
          file.path(dir_tab, "pooled_PC_FDR_candidates.csv"))
