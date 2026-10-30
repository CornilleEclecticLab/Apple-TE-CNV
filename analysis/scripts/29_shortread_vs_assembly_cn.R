# =============================================================================
# Script:  29_shortread_vs_assembly_cn.R
# Purpose: Validate the central variable of the paper. TE "copy number" is
#          estimated from short reads mapped to a consensus library; a reviewer
#          can object that a copy diverged from its consensus maps less well and
#          therefore looks less abundant, which would bias between-species
#          comparisons. We test this directly against the five annotated
#          assemblies, where copies are counted, not inferred.
#
#          Two tests:
#            A. across families, within each assembly: does a family with more
#               annotated copies carry more short reads?
#            B. across species, within each family: does the between-species
#               ranking of short-read abundance match the between-species
#               ranking of annotated copy number? This is the test that speaks
#               to the mapping-bias objection.
# Author:  Antho
# =============================================================================
source("scripts/00_setup.R")
suppressMessages({library(dplyr); library(readr); library(tidyr); library(tibble)})
TAG <- "29_shortread_vs_assembly_cn"
dir_tab <- file.path(dir_outputs, TAG, "tables")
dir_fig <- file.path(dir_outputs, TAG, "figures")
for (d in c(dir_tab, dir_fig)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

asm <- read_csv(file.path(dir_tab, "assembly_copy_counts.csv"), show_col_types = FALSE)
sr  <- read_csv(file.path(dir_tab, "shortread_species_means.csv"), show_col_types = FALSE)

ASM_SPECIES <- c(Mdo = "M. domestica", MdoR = "M. domestica (rootstock)",
                 Mori = "M. orientalis", Msie = "M. sieversii", Msyl = "M. sylvestris")

d <- asm %>% mutate(species = ASM_SPECIES[assembly]) %>%
  inner_join(sr, by = c("consensus", "species"))
message(sprintf("%d family x assembly pairs, %d families",
                nrow(d), n_distinct(d$consensus)))

# --- A. across families, within assembly -------------------------------------
A <- d %>% group_by(assembly) %>%
  summarise(n_families = n(),
            rho = cor(n_copies, mean_reads, method = "spearman"),
            p   = cor.test(n_copies, mean_reads, method = "spearman",
                           exact = FALSE)$p.value, .groups = "drop")
write_csv(A, file.path(dir_tab, "A_across_families.csv")); print(as.data.frame(A))

# --- B. across species, within family ----------------------------------------
B <- d %>% group_by(consensus, TE_type) %>%
  filter(n() >= 4, sd(n_copies) > 0, sd(mean_reads) > 0) %>%
  summarise(n_assemblies = n(),
            rho = cor(n_copies, mean_reads, method = "spearman"), .groups = "drop")
write_csv(B, file.path(dir_tab, "B_across_species.csv"))
message(sprintf("\nB: %d families testable; median rho = %.3f; %.1f%% positive",
                nrow(B), median(B$rho), 100*mean(B$rho > 0)))
w <- wilcox.test(B$rho, mu = 0, exact = FALSE)
message(sprintf("   Wilcoxon signed-rank against rho = 0: V = %.0f, p = %.3g",
                w$statistic, w$p.value))

# the objection is specifically about structural vs residual, which differ in
# identity to consensus, so the two sets are compared explicitly
bt <- B %>% group_by(TE_type) %>%
  summarise(n = n(), median_rho = median(rho), pct_positive = 100*mean(rho > 0), .groups = "drop")
write_csv(bt, file.path(dir_tab, "B_by_TE_type.csv")); print(as.data.frame(bt))
if (n_distinct(B$TE_type) == 2) {
  wt <- wilcox.test(rho ~ TE_type, B, exact = FALSE)
  message(sprintf("   structural vs residual rho: W = %.0f, p = %.3g", wt$statistic, wt$p.value))
}
log_session(TAG)

# --- C. the direct test of the objection -------------------------------------
# If a copy that has diverged from the consensus recruits fewer reads, then the
# families whose copies are least identical to the consensus must carry fewer
# reads per annotated copy. That is exactly what this test measures.
idt <- read_csv(file.path(dir_tab, "mean_identity.csv"), show_col_types = FALSE)
C <- d %>% inner_join(idt, by = c("assembly", "consensus")) %>%
  filter(n_copies > 0, mean_reads > 0) %>%
  mutate(reads_per_copy = mean_reads / n_copies)
res <- C %>% group_by(assembly) %>%
  summarise(n = n(),
            rho = cor(mean_identity, reads_per_copy, method = "spearman"),
            p   = cor.test(mean_identity, reads_per_copy, method = "spearman",
                           exact = FALSE)$p.value, .groups = "drop")
write_csv(res, file.path(dir_tab, "C_identity_vs_reads_per_copy.csv"))
message("\nC. identity to the consensus vs reads per annotated copy:")
print(as.data.frame(res))

# The identity contrast between the two categories is reported at the insertion
# level in Table S11 (residual 86.2% vs structural 81.0%, Welch t-test). It is
# not recomputed here at the family x assembly scale: that unit weights the
# families differently and flips the sign, which would contradict Table S11
# without measuring anything more.
rpc <- C %>% group_by(TE_type) %>%
  summarise(n = n(), median_reads_per_copy = median(reads_per_copy), .groups = "drop")
write_csv(rpc, file.path(dir_tab, "C_reads_per_copy_by_type.csv"))
print(as.data.frame(rpc))
message(sprintf("reads/copie, structural vs residual : p = %.3g",
        wilcox.test(reads_per_copy ~ TE_type, C, exact = FALSE)$p.value))
log_session(TAG)
