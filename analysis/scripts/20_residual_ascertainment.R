# =============================================================================
# Script:  20_residual_ascertainment.R
# Purpose: The residual TE set was not built as "TEs without structure". Methods
#          (Defining structural and residual TE consensus) state that structural
#          TEs were first removed, then the remaining TEs were screened with
#          SPECIES AND GEOGRAPHICAL GROUPS as supervised categories, retaining
#          those "that exhibited significant structure associated with species
#          identity and/or geographical distribution".
#
#          So residual TEs are ascertained ON geography. Two consequences:
#
#          (a) They should be depleted of the neutral SYLW/SYLE split, because
#              Q1 and Q5 were separate clusters in the K = 6 run used to define
#              structural TEs. Test 1 checks this.
#          (b) Their p-value distribution against climate cannot be used as a
#              calibration check, because climate covaries with the geography
#              they were selected on. Test 2 separates climate from geography by
#              putting latitude and longitude in the model.
#
# Author:  Antho
# Project: TE_Project / final_pipeline (Genome Biology revision)
# =============================================================================

source("scripts/00_setup.R")
library(dplyr); library(readr); library(tibble); library(purrr)

TAG <- "20_residual_ascertainment"
dir_out <- file.path(dir_outputs, TAG); dir_tab <- file.path(dir_out, "tables")
for (d in c(dir_out, dir_tab)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

PC_DEPTH <- 3
GROUP_OF <- function(x) dplyr::case_when(
  x == "Q1" ~ "SYLW", x == "Q5" ~ "SYLE",
  x %in% c("Q1_Admix", "Q5_Admix") ~ "SYL_Admix", TRUE ~ NA_character_)

counts <- read_csv(file.path(dir_data, "coord_syl_TE.csv"), show_col_types = FALSE) %>%
  mutate(Sample_ID = toupper(Sample_ID))
te_vars <- setdiff(names(counts), "Sample_ID")
sample_info <- read_csv(file.path(dir_data, "250923_sample_list_TE.csv"), show_col_types = FALSE) %>%
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
cm <- as.matrix(counts[match(design$Sample_ID, counts$Sample_ID), te_vars, drop = FALSE])
message(sprintf("%d individuals, %d residual TEs", nrow(design), length(te_vars)))

# =============================================================================
# Test 1 — are residual TEs purged of the neutral SYLW / SYLE split?
# =============================================================================
ew <- design$group_label %in% c("SYLW", "SYLE")
p_ew <- vapply(seq_along(te_vars), function(i) {
  y <- as.numeric(cm[ew, i]); g <- factor(design$group_label[ew])
  if (length(unique(y)) < 2) return(NA_real_)
  m <- suppressWarnings(tryCatch(glm(y ~ g, family = quasipoisson()), error = function(e) NULL))
  if (is.null(m)) return(NA_real_)
  cf <- summary(m)$coefficients
  if (nrow(cf) < 2) NA_real_ else cf[2, 4]
}, numeric(1))
p_ew <- p_ew[!is.na(p_ew)]
message(sprintf("\nTest 1 — residual TE CNV against the SYLW/SYLE split (n = %d TEs)", length(p_ew)))
message(sprintf("  p < 0.05 : %.3f   (0.05 expected if purged)", mean(p_ew < 0.05)))
message(sprintf("  median p : %.3f   (0.50 expected if purged)", median(p_ew)))
write_csv(tibble(p = p_ew), file.path(dir_tab, "residual_vs_EW_split.csv"))

# =============================================================================
# Test 2 — climate, once latitude and longitude are in the model
# =============================================================================
# Residual TEs were ascertained on geography, so a bare climate model cannot
# tell climate from the geography it was selected on. Adding lat/long asks
# whether climate explains anything beyond position.
fit_set <- function(extra, label) {
  res <- map_dfr(sprintf("bio%d", 1:19), function(bv) {
    preds <- c(bv, extra)
    ok <- complete.cases(design[, preds, drop = FALSE])
    md <- design[ok, , drop = FALSE]; cmv <- cm[ok, , drop = FALSE]
    fml <- as.formula(paste("te_count ~", paste(sprintf("`%s`", preds), collapse = " + ")))
    p <- rep(NA_real_, length(te_vars))
    for (i in seq_along(te_vars)) {
      y <- as.numeric(cmv[, i]); if (length(unique(y)) < 2) next
      dat <- md; dat$te_count <- y
      m <- suppressWarnings(tryCatch(glm(fml, data = dat, family = quasipoisson()),
                                     error = function(e) NULL))
      if (is.null(m)) next
      cf <- summary(m)$coefficients
      if (nrow(cf) >= 2) p[i] <- cf[2, 4]
    }
    tibble(variable = bv, TE_ID = te_vars, p = p)
  }) %>% filter(!is.na(p))
  tibble(model = label,
         prop_p05 = mean(res$p < 0.05), median_p = median(res$p),
         n_TE_p01 = n_distinct(res$TE_ID[res$p < 0.01]),
         n_TE_fdr = n_distinct(res$TE_ID[p.adjust(res$p, "BH") < 0.05]))
}
message("\nTest 2 — pooled models on all individuals:")
grid <- bind_rows(
  fit_set(character(0),                              "climate only"),
  fit_set(sprintf("PC%d", 1:PC_DEPTH),               "climate + PC1-3"),
  fit_set(c("latitude", "longitude"),                "climate + lat/long"),
  fit_set(c(sprintf("PC%d", 1:PC_DEPTH), "latitude", "longitude"),
                                                     "climate + PC1-3 + lat/long"))
write_csv(grid, file.path(dir_tab, "pooled_model_grid.csv"))
print(as.data.frame(grid))

log_session(TAG)
message("=== 20_residual_ascertainment.R finished ===")

# --- Export the candidate list for the climate + lat/long model --------------
res_ll <- purrr::map_dfr(sprintf("bio%d", 1:19), function(bv) {
  preds <- c(bv, "latitude", "longitude")
  ok <- complete.cases(design[, preds, drop = FALSE])
  md <- design[ok, , drop = FALSE]; cmv <- cm[ok, , drop = FALSE]
  fml <- as.formula(paste("te_count ~", paste(sprintf("`%s`", preds), collapse = " + ")))
  b <- p <- rep(NA_real_, length(te_vars))
  for (i in seq_along(te_vars)) {
    y <- as.numeric(cmv[, i]); if (length(unique(y)) < 2) next
    dat <- md; dat$te_count <- y
    m <- suppressWarnings(tryCatch(glm(fml, data = dat, family = quasipoisson()),
                                   error = function(e) NULL))
    if (is.null(m)) next
    cf <- summary(m)$coefficients
    if (nrow(cf) >= 2) { b[i] <- cf[2,1]; p[i] <- cf[2,4] }
  }
  tibble(variable = bv, TE_ID = te_vars, beta = b, p = p)
}) %>% filter(!is.na(p)) %>% mutate(fdr = p.adjust(p, "BH"))
write_csv(res_ll, file.path(dir_tab, "glm_pooled_latlong_all.csv"))
cand_ll <- res_ll %>% filter(fdr < 0.05) %>%
  group_by(TE_ID) %>%
  summarise(n_variables = n(), best_fdr = min(fdr),
            best_variable = variable[which.min(fdr)], .groups = "drop") %>%
  arrange(best_fdr)
write_csv(cand_ll, file.path(dir_tab, "candidates_latlong_FDR.csv"))
message(sprintf("\nclimate + lat/long, FDR < 0.05 : %d TEs", nrow(cand_ll)))
