# =============================================================================
# Script:  28_composition_final34.R
# Purpose: Re-run the sequence-composition analysis on the 34 replicated
#          candidates. Script 15 was written when the candidate set was the
#          union of 38; Table S15 and the Results numbers still come from it,
#          while Figure 3C is already drawn on the 34. This closes that gap.
# Author:  Antho
# =============================================================================
source("scripts/00_setup.R")
suppressMessages({library(dplyr); library(readr); library(tibble); library(Biostrings)})
TAG <- "28_composition_final34"
dir_tab <- file.path(dir_outputs, TAG, "tables")
dir.create(dir_tab, showWarnings = FALSE, recursive = TRUE)

mem  <- read_csv(file.path(dir_outputs,"18_climate_pc_diagnostics","tables",
                           "candidate_membership.csv"), show_col_types = FALSE)
CAND <- mem$TE_ID[mem$in_modelA]
message(sprintf("candidate set: %d TEs", length(CAND)))

base <- read_csv(file.path(dir_outputs,"15_sequence_composition","tables",
                           "composition_per_sequence.csv"), show_col_types = FALSE)
sf_of <- setNames(base$superfamily, base$TE_ID)
BG    <- base$TE_ID[base$grp == "background"]

db <- readDNAStringSet(file.path(dir_data,"TE_consensus_db.fasta"))
names(db) <- sub(" .*$","",names(db))
ids <- intersect(c(BG, CAND), names(db))
s   <- db[ids]
mono <- letterFrequency(s, c("A","C","G","T")); di <- dinucleotideFrequency(s); L <- width(s)
oe <- function(d) di[,d]*L/(mono[,substr(d,1,1)]*mono[,substr(d,2,2)])

DIN <- c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT")
comp <- tibble(TE_ID = names(s),
               length_bp = L,
               gc_percent = 100*(mono[,"G"]+mono[,"C"])/L,
               superfamily = sf_of[names(s)],
               grp = factor(ifelse(names(s) %in% CAND, "climate_associated", "background"),
                            levels = c("background","climate_associated")))
for (d in DIN) comp[[d]] <- oe(d)
comp <- comp %>% filter(if_all(all_of(DIN), is.finite))
message(sprintf("n background = %d, n climate = %d",
                sum(comp$grp=="background"), sum(comp$grp=="climate_associated")))
write_csv(comp, file.path(dir_tab,"composition_per_sequence_34.csv"))

# --- dinucleotide obs/exp tests ----------------------------------------------
din <- bind_rows(lapply(DIN, function(d) {
  x <- comp[[d]]; g <- comp$grp
  w <- wilcox.test(x ~ g, exact = FALSE)
  tibble(dinucleotide = d,
         oe_candidates = median(x[g=="climate_associated"]),
         oe_background = median(x[g=="background"]),
         difference    = median(x[g=="climate_associated"]) - median(x[g=="background"]),
         wilcoxon_p    = w$p.value)
})) %>% mutate(padj = p.adjust(wilcoxon_p, "BH"),
               role = ifelse(dinucleotide=="CG","CpG",
                      ifelse(dinucleotide %in% c("TG","CA"),"deamination product","other"))) %>%
  arrange(wilcoxon_p)
write_csv(din, file.path(dir_tab,"dinucleotide_oe_tests_34.csv"))
print(as.data.frame(din))

# --- length / GC -------------------------------------------------------------
ct <- bind_rows(lapply(c("length_bp","gc_percent"), function(v) {
  x <- comp[[v]]; g <- comp$grp
  tibble(variable = v,
         median_candidates = median(x[g=="climate_associated"]),
         median_background = median(x[g=="background"]),
         wilcoxon_p = wilcox.test(x ~ g, exact = FALSE)$p.value)
}))
write_csv(ct, file.path(dir_tab,"composition_tests_34.csv")); print(as.data.frame(ct))

# --- CpG within superfamily --------------------------------------------------
cs <- comp %>% group_by(superfamily) %>%
  filter(sum(grp=="climate_associated") >= 3, sum(grp=="background") >= 3) %>%
  summarise(n_climate = sum(grp=="climate_associated"),
            n_background = sum(grp=="background"),
            median_climate = median(CG[grp=="climate_associated"]),
            median_background = median(CG[grp=="background"]),
            p = wilcox.test(CG ~ grp, exact = FALSE)$p.value, .groups="drop") %>%
  mutate(padj = p.adjust(p,"BH"))
write_csv(cs, file.path(dir_tab,"cpg_per_superfamily_34.csv")); print(as.data.frame(cs))

# --- superfamily-adjusted model ----------------------------------------------
d2 <- comp %>% mutate(status = factor(grp, labels=c("background","climate")),
                      superfamily = factor(superfamily))
m_null <- lm(CG ~ superfamily, d2); m_add <- lm(CG ~ status + superfamily, d2)
m_int  <- lm(CG ~ status * superfamily, d2)
co <- summary(m_add)$coefficients["statusclimate", ]; ci <- confint(m_add)["statusclimate", ]
a_add <- anova(m_null, m_add); a_int <- anova(m_add, m_int)
res <- tibble(estimate = co[1], se = co[2], t = co[3], p = co[4],
              ci_low = ci[1], ci_high = ci[2],
              p_status_after_superfamily = a_add$`Pr(>F)`[2],
              p_interaction = a_int$`Pr(>F)`[2])
write_csv(res, file.path(dir_tab,"cpg_status_adjusted_effect_34.csv")); print(as.data.frame(res))

message(sprintf("\nCpG: median %.1f%% (climate) vs %.1f%% (background), p = %.3g; BH padj = %.3g",
        100*din$oe_candidates[din$dinucleotide=="CG"],
        100*din$oe_background[din$dinucleotide=="CG"],
        din$wilcoxon_p[din$dinucleotide=="CG"], din$padj[din$dinucleotide=="CG"]))
log_session(TAG)
