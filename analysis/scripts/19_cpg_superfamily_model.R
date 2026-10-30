# =============================================================================
# Script:  19_cpg_superfamily_model.R
# Purpose: Test whether the CpG difference between climate-associated and
#          background residual TEs survives explicit control for TE superfamily.
#
#          The manuscript reports a global Wilcoxon test plus separate
#          within-superfamily tests, but the Abstract calls the comparison
#          "superfamily-matched", which is a stronger claim than a global test
#          followed by stratified ones. This script fits the model the claim
#          would need:
#              CpG observed/expected ~ climate status + superfamily
#              CpG observed/expected ~ climate status * superfamily
#          and compares them, so the wording can follow the evidence.
#
# Author:  Antho
# Project: TE_Project / final_pipeline (Genome Biology revision)
# =============================================================================

source("scripts/00_setup.R")
library(dplyr); library(readr); library(tibble); library(ggplot2)

TAG     <- "19_cpg_superfamily_model"
dir_out <- file.path(dir_outputs, TAG)
dir_tab <- file.path(dir_out, "tables"); dir_fig <- file.path(dir_out, "figures")
for (d in c(dir_out, dir_tab, dir_fig)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

comp <- read_csv(file.path(dir_outputs, "15_sequence_composition", "tables",
                           "composition_per_sequence.csv"), show_col_types = FALSE) %>%
  filter(is.finite(cpg_obs_exp)) %>%
  mutate(status = factor(grp, levels = c("background", "climate_associated"),
                         labels = c("background", "climate")),
         superfamily = factor(superfamily))
message(sprintf("n = %d sequences", nrow(comp)))
print(table(comp$status, comp$superfamily))

# --- 1. the global test currently reported -----------------------------------
w <- wilcox.test(cpg_obs_exp ~ status, data = comp, exact = FALSE)
message(sprintf("\nGlobal Wilcoxon (as reported): W = %.0f, p = %.3g", w$statistic, w$p.value))

# --- 2. additive model: does status survive control for superfamily? ---------
m_null <- lm(cpg_obs_exp ~ superfamily, data = comp)
m_add  <- lm(cpg_obs_exp ~ status + superfamily, data = comp)
m_int  <- lm(cpg_obs_exp ~ status * superfamily, data = comp)

a_add <- anova(m_null, m_add)      # status, after superfamily
a_int <- anova(m_add,  m_int)      # interaction, after both main effects

co <- summary(m_add)$coefficients["statusclimate", ]
ci <- confint(m_add)["statusclimate", ]
res <- tibble(
  test = c("status | superfamily (additive)", "status x superfamily (interaction)"),
  df   = c(a_add$Df[2], a_int$Df[2]),
  F    = c(a_add$F[2], a_int$F[2]),
  p    = c(a_add$`Pr(>F)`[2], a_int$`Pr(>F)`[2]))
write_csv(res, file.path(dir_tab, "cpg_superfamily_anova.csv"))
message("\nSequential F tests:"); print(as.data.frame(res))
message(sprintf("\nAdjusted effect of climate status (additive model): %+.4f (95%% CI %.4f to %.4f), p = %.3g",
                co[1], ci[1], ci[2], co[4]))
write_csv(tibble(estimate = co[1], se = co[2], t = co[3], p = co[4],
                 ci_low = ci[1], ci_high = ci[2]),
          file.path(dir_tab, "cpg_status_adjusted_effect.csv"))

# --- 3. a rank-based equivalent, since CpG o/e is skewed ---------------------
# Same question without the normality assumption: rank the values, then test.
comp$rk <- rank(comp$cpg_obs_exp)
m_rk_null <- lm(rk ~ superfamily, data = comp)
m_rk_add  <- lm(rk ~ status + superfamily, data = comp)
a_rk <- anova(m_rk_null, m_rk_add)
message(sprintf("Rank-based check, status after superfamily: F = %.2f, p = %.3g",
                a_rk$F[2], a_rk$`Pr(>F)`[2]))
write_csv(tibble(F = a_rk$F[2], p = a_rk$`Pr(>F)`[2]),
          file.path(dir_tab, "cpg_status_rank_based.csv"))

# --- 4. per-superfamily effect, for the interaction picture ------------------
per_fam <- comp %>% group_by(superfamily) %>%
  group_modify(~{
    x <- .x$cpg_obs_exp[.x$status == "climate"]
    y <- .x$cpg_obs_exp[.x$status == "background"]
    if (length(x) < 3) return(tibble(n_climate = length(x), n_background = length(y),
                                     median_climate = median(x), median_background = median(y),
                                     difference = median(x) - median(y), p = NA_real_))
    tt <- wilcox.test(x, y, exact = FALSE)
    tibble(n_climate = length(x), n_background = length(y),
           median_climate = median(x), median_background = median(y),
           difference = median(x) - median(y), p = tt$p.value)
  }) %>% ungroup() %>% mutate(padj = p.adjust(p, "BH"))
write_csv(per_fam, file.path(dir_tab, "cpg_per_superfamily.csv"))
message("\nPer superfamily:"); print(as.data.frame(per_fam))

log_session(TAG)
message("=== 19_cpg_superfamily_model.R finished ===")
