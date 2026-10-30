# =============================================================================
# Script:  22_optionC_check.R
# Purpose: Option C in full. The candidate set becomes the residual TEs whose
#          climate association survives latitude and longitude, at FDR < 0.05
#          (n = 82). Everything the figure needs is recomputed for that set:
#          vector strength, CpG composition, and the superfamily-adjusted test.
#          Nothing is carried over from the 38, because only 9 of the 82 were
#          ever in it.
# Author:  Antho
# =============================================================================
source("scripts/00_setup.R")
suppressMessages({library(dplyr); library(readr); library(tibble); library(purrr)
                  library(Biostrings); library(glmmTMB)})

TAG <- "22_optionC_check"
dir_out <- file.path(dir_outputs, TAG); dir_tab <- file.path(dir_out, "tables")
dir.create(dir_tab, showWarnings = FALSE, recursive = TRUE)

cand <- read_csv(file.path(dir_outputs, "20_residual_ascertainment", "tables",
                           "candidates_latlong_FDR.csv"), show_col_types = FALSE)$TE_ID
counts <- read_csv(file.path(dir_data, "coord_syl_TE.csv"), show_col_types = FALSE) %>%
  mutate(Sample_ID = toupper(Sample_ID))
te_all <- setdiff(names(counts), "Sample_ID")
message(sprintf("candidates: %d | residual TEs: %d", length(cand), length(te_all)))

meta_spatial <- read_csv(file.path(dir_data, "250923_sample_list_TE.csv"), show_col_types = FALSE) %>%
  filter(Species == "M. sylvestris", !is.na(latitude), !is.na(longitude)) %>%
  transmute(Sample_ID = toupper(Sample_ID), latitude, longitude)

# --- 1. Vector strength for the candidates -----------------------------------
fit_vs <- function(te_id) {
  d <- counts %>% select(Sample_ID, Count = all_of(te_id)) %>%
    inner_join(meta_spatial, by = "Sample_ID") %>% tidyr::drop_na()
  fit <- try(glmmTMB(Count ~ latitude + longitude, data = d, family = nbinom2()), silent = TRUE)
  if (inherits(fit, "try-error")) return(tibble(TE_ID = te_id, true_vs = NA_real_))
  fx <- try(fixef(fit)$cond, silent = TRUE)
  if (inherits(fx, "try-error") || !all(c("latitude","longitude") %in% names(fx)))
    return(tibble(TE_ID = te_id, true_vs = NA_real_))
  ang <- atan2(fx[["latitude"]], fx[["longitude"]])
  proj <- d$longitude * cos(ang) + d$latitude * sin(ang)
  tibble(TE_ID = te_id, beta_lat = fx[["latitude"]], beta_long = fx[["longitude"]],
         angle_rad = ang, true_vs = abs(suppressWarnings(cor(d$Count, proj, use = "complete.obs"))))
}
message("Fitting spatial models ...")
vs <- map_dfr(cand, fit_vs) %>% filter(!is.na(true_vs)) %>%
  mutate(superfamily = case_when(grepl("Copia$", TE_ID) ~ "LTR_Copia",
                                 grepl("Gypsy$", TE_ID) ~ "LTR_Gypsy",
                                 grepl("TIR$", TE_ID) ~ "TIR",
                                 grepl("noCat$", TE_ID) ~ "noCat", TRUE ~ "LTR"))
write_csv(vs, file.path(dir_tab, "VS_optionC.csv"))
message(sprintf("  VS computed for %d; %d have VS > 0.5 (%.0f%%)",
                nrow(vs), sum(vs$true_vs > .5), 100*mean(vs$true_vs > .5)))

# --- 2. CpG composition for every residual TE --------------------------------
db <- readDNAStringSet(file.path(dir_data, "TE_consensus_db.fasta"))
names(db) <- sub(" .*$", "", names(db))
have <- intersect(te_all, names(db))
message(sprintf("Sequences found for %d of %d residual TEs", length(have), length(te_all)))
s  <- db[have]
nC <- as.numeric(letterFrequency(s, "C")); nG <- as.numeric(letterFrequency(s, "G"))
L  <- width(s); nCG <- as.numeric(vcountPattern("CG", s))
comp <- tibble(TE_ID = have, length_bp = L,
               gc = 100*(nC+nG)/L,
               cpg_oe = ifelse(nC*nG > 0, nCG*L/(nC*nG), NA_real_)) %>%
  filter(is.finite(cpg_oe)) %>%
  mutate(status = factor(ifelse(TE_ID %in% cand, "climate", "background"),
                         levels = c("background","climate")),
         superfamily = factor(case_when(grepl("Copia$", TE_ID) ~ "Copia",
                                        grepl("Gypsy$", TE_ID) ~ "Gypsy",
                                        grepl("TIR$", TE_ID) ~ "TIR",
                                        grepl("noCat$", TE_ID) ~ "noCat", TRUE ~ "LTR")))
write_csv(comp, file.path(dir_tab, "composition_optionC.csv"))
print(table(comp$status, comp$superfamily))

w  <- wilcox.test(cpg_oe ~ status, data = comp, exact = FALSE)
a  <- anova(lm(cpg_oe ~ superfamily, comp), lm(cpg_oe ~ status + superfamily, comp))
ai <- anova(lm(cpg_oe ~ status + superfamily, comp), lm(cpg_oe ~ status * superfamily, comp))
co <- summary(lm(cpg_oe ~ status + superfamily, comp))$coefficients["statusclimate", ]
res <- tibble(n_climate = sum(comp$status=="climate"),
              n_background = sum(comp$status=="background"),
              median_climate = median(comp$cpg_oe[comp$status=="climate"]),
              median_background = median(comp$cpg_oe[comp$status=="background"]),
              wilcoxon_p = w$p.value,
              adjusted_effect = co[1], adjusted_p = a$`Pr(>F)`[2],
              interaction_p = ai$`Pr(>F)`[2])
write_csv(res, file.path(dir_tab, "cpg_optionC.csv"))
message("\nCpG under option C:"); print(as.data.frame(res))
log_session(TAG)
