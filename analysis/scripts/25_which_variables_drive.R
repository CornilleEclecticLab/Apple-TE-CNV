# =============================================================================
# Script:  25_which_variables_drive.R
# Purpose: R2.5 asks whether analysing three subsets controls for population
#          structure. It does not control for confounding shared by all three
#          strata, and climate, geography and structure are collinear here. So
#          the useful question is narrower: are the candidates driven by the
#          bioclimatic variables that are MOST confounded with structure, or by
#          the ones that are LEAST confounded?
#
#          Each of the 19 variables has a different R^2 on PC1-3 within each
#          subset. If the candidate-driving cells sit preferentially at low
#          R^2, the association is unlikely to be a by-product of structure.
#          If they sit at high R^2, it probably is.
#
#          The same logic applies to the 34 -> 6 attenuation: if the loss of
#          effect is concentrated in the collinear variables, that is
#          confounding being removed; if it is uniform, it is the cost of
#          spending three degrees of freedom.
# Author:  Antho
# =============================================================================
source("scripts/00_setup.R")
suppressMessages({library(dplyr); library(readr); library(tibble)})
TAG <- "25_which_variables_drive"
dir_tab <- file.path(dir_outputs, TAG, "tables")
dir.create(dir_tab, showWarnings = FALSE, recursive = TRUE)

fits   <- read_csv(file.path(dir_outputs,"18_climate_pc_diagnostics","tables",
                             "glm_both_models_with_SE.csv"), show_col_types=FALSE)
collin <- read_csv(file.path(dir_outputs,"18_climate_pc_diagnostics","tables",
                             "collinearity_bio_vs_PC.csv"), show_col_types=FALSE)
mem    <- read_csv(file.path(dir_outputs,"18_climate_pc_diagnostics","tables",
                             "candidate_membership.csv"), show_col_types=FALSE)
s34 <- mem$TE_ID[mem$in_modelA]

d <- fits %>% left_join(collin %>% select(group, variable, r2_bio_on_PC),
                        by = c("group","variable"))

# --- 1. collinearity of the cells that drive the candidates ------------------
drivers <- d %>% filter(TE_ID %in% s34, !is.na(p_A), p_A < 0.01)
allcell <- d %>% filter(!is.na(p_A))
cmp <- tibble(
  set = c("cells driving the 34 candidates", "all TE x variable x group cells"),
  n = c(nrow(drivers), nrow(allcell)),
  median_r2 = c(median(drivers$r2_bio_on_PC, na.rm=TRUE),
                median(allcell$r2_bio_on_PC, na.rm=TRUE)),
  mean_r2 = c(mean(drivers$r2_bio_on_PC, na.rm=TRUE),
              mean(allcell$r2_bio_on_PC, na.rm=TRUE)))
w <- wilcox.test(drivers$r2_bio_on_PC, allcell$r2_bio_on_PC, exact=FALSE)
write_csv(cmp, file.path(dir_tab, "driver_collinearity.csv"))
message("Collinearity with PC1-3 of the cells driving the candidates:")
print(as.data.frame(cmp))
message(sprintf("  Wilcoxon drivers vs all cells: p = %.3g", w$p.value))

# per group, since SYLE is the collinear one
per_g <- bind_rows(
  drivers %>% group_by(group) %>% summarise(set="drivers", n=n(),
              median_r2=median(r2_bio_on_PC), .groups="drop"),
  allcell %>% group_by(group) %>% summarise(set="all", n=n(),
              median_r2=median(r2_bio_on_PC), .groups="drop")) %>% arrange(group, set)
write_csv(per_g, file.path(dir_tab, "driver_collinearity_by_group.csv"))
message("\nBy group:"); print(as.data.frame(per_g))

# --- 2. is the attenuation concentrated in the collinear variables? ----------
att <- drivers %>%
  filter(!is.na(beta_PC), !is.na(beta_A), beta_A != 0, se_A > 0) %>%
  mutate(beta_ratio = abs(beta_PC)/abs(beta_A), se_ratio = se_PC/se_A)
ct_b <- cor.test(att$r2_bio_on_PC, att$beta_ratio, method="spearman", exact=FALSE)
ct_s <- cor.test(att$r2_bio_on_PC, att$se_ratio,   method="spearman", exact=FALSE)
res <- tibble(
  relation = c("R2(bio~PC) vs |beta_PC|/|beta_A|", "R2(bio~PC) vs SE_PC/SE_A"),
  rho = c(ct_b$estimate, ct_s$estimate),
  p   = c(ct_b$p.value,  ct_s$p.value))
write_csv(res, file.path(dir_tab, "attenuation_vs_collinearity.csv"))
message("\nIs attenuation concentrated where climate and structure are collinear?")
print(as.data.frame(res))

# effect by tercile of collinearity, easier to read
terc <- att %>% mutate(band = cut(r2_bio_on_PC, quantile(r2_bio_on_PC, c(0,1/3,2/3,1)),
                                  include.lowest=TRUE,
                                  labels=c("low collinearity","medium","high collinearity"))) %>%
  group_by(band) %>% summarise(n=n(), median_beta_ratio=median(beta_ratio),
                               median_se_ratio=median(se_ratio), .groups="drop")
write_csv(terc, file.path(dir_tab, "attenuation_by_tercile.csv"))
message("\nBy tercile of collinearity:"); print(as.data.frame(terc))
log_session(TAG)
