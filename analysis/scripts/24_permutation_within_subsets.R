# =============================================================================
# Script:  24_permutation_within_subsets.R
# Purpose: Answer R2.5 empirically rather than by argument.
#
#          The reviewer asks how analysing three subsets controls for
#          population structure. The direct test is to measure the false
#          positive rate of the actual rule: permute the climatic values
#          AMONG INDIVIDUALS WITHIN EACH SUBSET, which destroys any real
#          climate association while leaving intact
#            - the correlation structure among TEs,
#            - the residual genetic structure inside each subset,
#            - the group sizes and the replication requirement,
#          then run the published rule and count the candidates it returns.
#
#          If stratification controls confounding, the permuted rule should
#          return almost nothing, whatever structure remains within subsets.
#
# Author:  Antho
# =============================================================================
source("scripts/00_setup.R")
suppressMessages({library(dplyr); library(readr); library(tibble); library(purrr)})
set.seed(20260824)

TAG <- "24_permutation_within_subsets"
dir_tab <- file.path(dir_outputs, TAG, "tables")
dir.create(dir_tab, showWarnings = FALSE, recursive = TRUE)

N_PERM <- 20
THR    <- 0.01
GROUPS <- c("SYLW", "SYLE", "SYL_Admix")
GROUP_OF <- function(x) case_when(x == "Q1" ~ "SYLW", x == "Q5" ~ "SYLE",
  x %in% c("Q1_Admix","Q5_Admix") ~ "SYL_Admix", TRUE ~ NA_character_)

counts <- read_csv(file.path(dir_data,"coord_syl_TE.csv"), show_col_types=FALSE) %>%
  mutate(Sample_ID = toupper(Sample_ID))
te_vars <- setdiff(names(counts), "Sample_ID")
si <- read_csv(file.path(dir_data,"250923_sample_list_TE.csv"), show_col_types=FALSE) %>%
  filter(Species=="M. sylvestris", !is.na(latitude), !is.na(longitude)) %>%
  transmute(Sample_ID=toupper(Sample_ID), latitude, longitude,
            group_label=GROUP_OF(group_structure_CORE_admix_K6)) %>% filter(!is.na(group_label))
clim <- file.path(dir_data,"worldclim_env_data","climate","wc2.1_10m")
bio <- terra::rast(file.path(clim, sprintf("wc2.1_10m_bio_%d.tif", 1:19)))
names(bio) <- sprintf("bio%d", 1:19)
si <- bind_cols(si, as_tibble(terra::extract(bio, as.matrix(si[,c("longitude","latitude")]))))
design <- si %>% filter(Sample_ID %in% counts$Sample_ID) %>% arrange(Sample_ID)
cm <- as.matrix(counts[match(design$Sample_ID, counts$Sample_ID), te_vars, drop=FALSE])
message(sprintf("%d individuals, %d TEs", nrow(design), length(te_vars)))

#' The published rule: significant at p < THR for at least one bioclimatic
#' variable in every one of the three subsets, not necessarily the same one.
run_rule <- function(dsg) {
  per <- lapply(GROUPS, function(g) {
    idx <- which(dsg$group_label == g)
    d <- dsg[idx, , drop=FALSE]; cmv <- cm[idx, , drop=FALSE]
    hit <- rep(FALSE, length(te_vars))
    for (bv in sprintf("bio%d", 1:19)) {
      x <- d[[bv]]; ok <- !is.na(x)
      if (sum(ok) < 8) next
      for (i in seq_along(te_vars)) {
        if (hit[i]) next
        y <- as.numeric(cmv[ok, i]); if (length(unique(y)) < 2) next
        m <- suppressWarnings(tryCatch(glm(y ~ x[ok], family=quasipoisson()),
                                       error=function(e) NULL))
        if (is.null(m)) next
        cf <- summary(m)$coefficients
        if (nrow(cf) >= 2 && !is.na(cf[2,4]) && cf[2,4] < THR) hit[i] <- TRUE
      }
    }
    te_vars[hit]
  })
  Reduce(intersect, per)
}

obs <- run_rule(design)
message(sprintf("\nObserved candidates under the published rule: %d", length(obs)))

bio_cols <- sprintf("bio%d", 1:19)
null_counts <- integer(N_PERM)
for (k in seq_len(N_PERM)) {
  d <- design
  for (g in GROUPS) {                       # permute climate WITHIN each subset
    idx <- which(d$group_label == g)
    d[idx, bio_cols] <- d[sample(idx), bio_cols]
  }
  null_counts[k] <- length(run_rule(d))
  message(sprintf("  permutation %2d/%d : %d candidates", k, N_PERM, null_counts[k]))
}
res <- tibble(permutation = seq_len(N_PERM), n_candidates = null_counts)
write_csv(res, file.path(dir_tab, "null_counts.csv"))
summ <- tibble(observed = length(obs), n_perm = N_PERM,
               null_mean = mean(null_counts), null_median = median(null_counts),
               null_max = max(null_counts),
               empirical_p = (sum(null_counts >= length(obs)) + 1) / (N_PERM + 1))
write_csv(summ, file.path(dir_tab, "summary.csv"))
message("\nFalse positive rate of the three-subset replication rule:")
print(as.data.frame(summ))
log_session(TAG)
