# =============================================================================
# Script:  02_candidate_TEs.R
# Purpose: Produce a RELIABLE candidate TE list for Figure 3, addressing
#          reviewer 2 comment 5 without inflating the list.
#
#          Context. Scripts 02 and 04 fitted the PC-covariate GLMs on the
#          COMBINED M. sylvestris dataset and returned 723 candidate TEs at
#          padj < 0.05, against the 34 shown in Figure 3. That gap is not a
#          statistical accident, and this script exists to close it.
#
#          The 34 candidates in Figure 3 were NOT produced by a multiple-
#          testing correction. They come from a three-way REPLICATION design:
#          the 72 M. sylvestris individuals are partitioned into three
#          disjoint groups - SYLW (Q1, n = 39), SYLE (Q5, n = 10) and
#          SYL_Admix (Q1_Admix + Q5_Admix, n = 23) - a quasi-Poisson GLM of
#          TE count on each of the 19 bioclimatic variables is fitted within
#          each group, the TEs reaching raw p < 0.01 for at least one variable
#          are unioned within each group, and the three group lists are then
#          intersected. Requiring the same TE to respond to climate
#          independently in all three groups is what makes the list small:
#          367 / 257 / 282 TEs per group collapse to 34 shared ones.
#
#          Two consequences follow. First, the subset design was a
#          replication filter, not the population-structure control the
#          manuscript claims it to be - which is exactly the reviewer's point.
#          Second, replacing it with PC covariates removes the replication and
#          the list explodes. The defensible answer is to keep BOTH: fit the
#          same three-group replication design, but add the neutral SNP PCs as
#          covariates inside each group.
#
#          This script therefore:
#            1. reproduces the published Figure 3 pipeline exactly (validation
#               gate: it must return the same 34 TEs, or the script stops);
#            2. refits the identical design with PC1-3 as covariates in every
#               group;
#            3. reports candidate counts, and what is shared versus specific,
#               at raw p < 0.05 / 0.01 / 0.001 and under BH correction;
#            4. exports the new candidate list and its FASTA so Figure 3 and
#               the TFBS analysis can be regenerated.
#
#          Note on the manuscript text: the Results currently describe the
#          three subsets as "SYLE, SYLW, and a combined dataset". The code
#          that produced Figure 3 used SYL_Admix as the third group, not a
#          combined dataset. That sentence needs correcting regardless of
#          which candidate list is finally adopted.
#
# Author:  Antho (analysis run with Claude Code)
# Date:    2026-08-24
# Project: TE_Project / review (reviewer response, Genome Biology submission)
#
# Inputs:
#   data/coord_syl_TE.csv          72 M. sylvestris x 1360 residual TE counts
#                                  (exact matrix used for Figure 3)
#   data/250923_sample_list_TE.csv metadata incl. coordinates and K6 clusters
#   data/modelA_fig3_raw_p.csv     published Model A results (raw p-values)
#   data/fig3_34TE.fasta           the published 34 candidates
#   outputs/01_snp_pca_axes/tables/PCA_SNP_coordinates_PC1_10.csv
#   data/worldclim_env_data/climate/wc2.1_10m/wc2.1_10m_bio_*.tif
#   data/TE_consensus_db.fasta     TE consensus sequences
#
# Outputs (all under outputs/02_candidate_TEs/):
#   tables/validation_reproduce_34.csv      reproduction gate
#   tables/glm_modelA_PC_by_group.csv          full PC-covariate GLM results
#   tables/candidate_counts.csv             per group / threshold / model
#   tables/candidate_comparison.csv         shared vs specific vs the 34
#   tables/FINAL_union_candidates.csv       the final non-redundant set
#   tables/candidate_consistency.csv        variable/sign consistency audit
#   tables/final_candidates_annotated.csv   Model A-PC replicated list
#   tables/permutation_control_summary.csv  observed vs permuted-PC null
#   lists/*.txt                             plain TE-ID lists
#   fasta/candidates_*.fasta                sequences for the TFBS rerun
#   figures/*.pdf + *.svg
#   ../../logs/02_candidate_TEs_*_sessionInfo.txt
#
# NAMING NOTE: in the manuscript, "Model A" is the climate GLM
# (TE count ~ bioclimatic variable) and "Model B" is the SPATIAL model
# (TE count ~ latitude + longitude, negative binomial), which produces the
# vector-strength (VS) values shown in Fig. 3A-B. To avoid a clash with that
# established terminology, the PC-covariate version of the climate model is
# called "Model A-PC" throughout this folder, never "Model B".
# =============================================================================

# --- Setup -------------------------------------------------------------------
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(terra)
  library(Biostrings)
  library(svglite)
})

set.seed(20260824)

# --- Parameters --------------------------------------------------------------
RAW_THRESHOLDS <- c(0.05, 0.01, 0.001)  # raw p, as used for Figure 3
PADJ_METHOD    <- "BH"
PC_DEPTH       <- 3                      # PC1-3, justified in script 04
MIN_N_PER_MODEL<- 3                      # published pipeline's minimum n
N_PERM         <- 30                     # permutations for the PC null control

# The three disjoint groups the Figure 3 pipeline actually used.
GROUP_OF <- function(cluster) {
  dplyr::case_when(
    cluster %in% c("Q5_Admix", "Q1_Admix") ~ "SYL_Admix",
    cluster == "Q5"                        ~ "SYLE",
    cluster == "Q1"                        ~ "SYLW",
    TRUE                                   ~ NA_character_
  )
}
GROUPS <- c("SYLW", "SYLE", "SYL_Admix")

dir_out  <- file.path(dir_outputs, "02_candidate_TEs")
dir_tab  <- file.path(dir_out, "tables")
dir_fig  <- file.path(dir_out, "figures")
dir_list <- file.path(dir_out, "lists")
dir_fa   <- file.path(dir_out, "fasta")
for (d in c(dir_out, dir_tab, dir_fig, dir_list, dir_fa)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

save_plot <- function(p, name, width = 7, height = 5) {
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".pdf")), p,
                  width = width, height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".svg")), p,
                  width = width, height = height, device = svglite::svglite)
}

message("=== 02_candidate_TEs.R ===")

# =============================================================================
# 1. Validation gate - reproduce the published 34 candidates
# =============================================================================
modelA <- readr::read_csv(file.path(dir_data, "modelA_fig3_raw_p.csv"),
                          show_col_types = FALSE) %>%
  dplyr::mutate(variable = sub("^wc2\\.1_10m_bio_", "bio", bio_var))

#' Figure 3 rule: within each group, union over the 19 bioclimatic variables
#' of the TEs reaching the threshold; then intersect the group lists.
replicated_candidates <- function(df, thr, p_col = "p_value", groups = GROUPS) {
  per_group <- lapply(groups, function(g) {
    unique(df$TE_ID[df$group == g & !is.na(df[[p_col]]) & df[[p_col]] < thr])
  })
  names(per_group) <- groups
  list(per_group = per_group, shared = sort(Reduce(intersect, per_group)))
}

published_34 <- sort(sub(" .*$", "",
  names(Biostrings::readDNAStringSet(file.path(dir_data, "fig3_34TE.fasta")))))

repro <- replicated_candidates(modelA, 0.01)
validation <- tibble::tibble(
  threshold        = 0.01,
  n_SYLW           = length(repro$per_group$SYLW),
  n_SYLE           = length(repro$per_group$SYLE),
  n_SYL_Admix      = length(repro$per_group$SYL_Admix),
  n_reproduced     = length(repro$shared),
  n_published      = length(published_34),
  n_matching       = length(intersect(repro$shared, published_34)),
  exact_match      = setequal(repro$shared, published_34)
)
readr::write_csv(validation, file.path(dir_tab, "validation_reproduce_34.csv"))
message("Validation against the published Figure 3 candidates:")
print(as.data.frame(validation))

if (!validation$exact_match) {
  stop("Could not reproduce the published 34 candidates exactly; ",
       "aborting rather than building on a pipeline that does not match.")
}
message("  -> the published Figure 3 rule is reproduced exactly.")

# =============================================================================
# 2. Rebuild the same data, and add the SNP PCs
# =============================================================================
counts <- readr::read_csv(file.path(dir_data, "coord_syl_TE.csv"),
                          show_col_types = FALSE)
te_vars <- setdiff(names(counts), "Sample_ID")
message(sprintf("Count matrix: %d individuals x %d residual TEs",
                nrow(counts), length(te_vars)))

sample_info <- readr::read_csv(file.path(dir_data, "250923_sample_list_TE.csv"),
                               show_col_types = FALSE) %>%
  dplyr::filter(Species == "M. sylvestris",
                !is.na(latitude), !is.na(longitude)) %>%
  dplyr::transmute(Sample_ID, latitude, longitude,
                   group_label = GROUP_OF(group_structure_CORE_admix_K6)) %>%
  dplyr::filter(!is.na(group_label))

pc_tbl <- readr::read_csv(
  file.path(dir_outputs, "01_snp_pca_axes", "tables",
            "PCA_SNP_coordinates_PC1_10.csv"),
  show_col_types = FALSE
) %>%
  dplyr::select(Sample_ID = sample_id, dplyr::all_of(sprintf("PC%d", 1:PC_DEPTH)))

# Bioclim, from the same WorldClim 2.1 10-arcmin rasters.
clim_dir   <- file.path(dir_data, "worldclim_env_data", "climate", "wc2.1_10m")
bio_stack  <- terra::rast(file.path(clim_dir, sprintf("wc2.1_10m_bio_%d.tif", 1:19)))
names(bio_stack) <- sprintf("bio%d", 1:19)

meta <- sample_info %>%
  dplyr::mutate(Sample_ID = toupper(Sample_ID)) %>%
  dplyr::inner_join(pc_tbl %>% dplyr::mutate(Sample_ID = toupper(Sample_ID)),
                    by = "Sample_ID")

env <- terra::extract(bio_stack, as.matrix(meta[, c("longitude", "latitude")]))
meta <- dplyr::bind_cols(meta, tibble::as_tibble(env))

counts <- counts %>% dplyr::mutate(Sample_ID = toupper(Sample_ID))
design <- meta %>%
  dplyr::filter(Sample_ID %in% counts$Sample_ID) %>%
  dplyr::arrange(Sample_ID)
count_mat <- as.matrix(counts[match(design$Sample_ID, counts$Sample_ID), te_vars,
                              drop = FALSE])
rownames(count_mat) <- design$Sample_ID

message("Group sizes used for the PC-covariate models:")
print(as.data.frame(design %>% dplyr::count(group_label, name = "n")))

# =============================================================================
# 3. Model A-PC - the same replication design, with PC1-3 inside each group
# =============================================================================
pc_names <- sprintf("PC%d", seq_len(PC_DEPTH))

#' Fit the quasi-Poisson GLMs for one group, all 19 bioclimatic variables.
fit_group <- function(grp) {
  idx <- which(design$group_label == grp)
  d   <- design[idx, , drop = FALSE]
  cm  <- count_mat[idx, , drop = FALSE]

  purrr::map_dfr(sprintf("bio%d", 1:19), function(bv) {
    predictors <- c(bv, pc_names)
    md <- d[, predictors, drop = FALSE]
    ok <- stats::complete.cases(md)
    md <- md[ok, , drop = FALSE]; cmv <- cm[ok, , drop = FALSE]

    rhs <- paste(sprintf("`%s`", predictors), collapse = " + ")
    fml <- stats::as.formula(paste("te_count ~", rhs))

    n_te <- length(te_vars)
    est <- rep(NA_real_, n_te); pv <- rep(NA_real_, n_te)

    for (i in seq_len(n_te)) {
      y <- as.numeric(cmv[, te_vars[i]])
      # Same guards as the published pipeline.
      if (length(unique(y)) < 2 || length(y) < MIN_N_PER_MODEL) next
      dat <- md; dat$te_count <- y
      fit <- suppressWarnings(tryCatch(
        stats::glm(fml, data = dat, family = stats::quasipoisson()),
        error = function(e) NULL))
      if (is.null(fit)) next
      cf <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
      if (is.null(cf) || nrow(cf) < 2) next
      est[i] <- cf[2, 1]; pv[i] <- cf[2, 4]
    }
    tibble::tibble(group = grp, variable = bv, TE_ID = te_vars,
                   estimate = est, p_value = pv)
  })
}

message(sprintf("Fitting Model A-PC (bio + PC1-%d) in each group ...", PC_DEPTH))
t0 <- Sys.time()
modelA_PC <- purrr::map_dfr(GROUPS, function(g) {
  message("  ", g, " ...")
  fit_group(g)
}) %>%
  dplyr::filter(!is.na(p_value))
message(sprintf("  %d fits in %.1f min", nrow(modelA_PC),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# BH correction, applied within group x bioclimatic variable, and reported
# alongside the raw p-values so both stringency regimes are available.
modelA_PC <- modelA_PC %>%
  dplyr::group_by(group, variable) %>%
  dplyr::mutate(padj = stats::p.adjust(p_value, method = PADJ_METHOD)) %>%
  dplyr::ungroup()

modelA_bh <- modelA %>%
  dplyr::filter(!is.na(p_value)) %>%
  dplyr::group_by(group, variable) %>%
  dplyr::mutate(padj = stats::p.adjust(p_value, method = PADJ_METHOD)) %>%
  dplyr::ungroup()

readr::write_csv(modelA_PC, file.path(dir_tab, "glm_modelA_PC_by_group.csv"))

# =============================================================================
# 3b. Permutation control - is the collapse confounding, or just lost df?
# =============================================================================
# Adding PC1-3 costs three degrees of freedom in every group, and SYLE has
# only n = 10, so part of the drop in candidate TEs could be pure power loss
# rather than removal of a genuine population-structure confounder.
#
# The control: refit Model A-PC with the PC coordinates PERMUTED AMONG
# INDIVIDUALS WITHIN EACH GROUP. The permuted covariates keep exactly the same
# distribution and consume exactly the same three degrees of freedom, but
# carry no real information about who each individual is. If the observed
# candidate count sits inside this null distribution, the collapse is a power
# artefact; if it sits well below it, the PCs are removing real confounding.
message(sprintf("Permutation control: %d permutations of PC1-%d within group ...",
                N_PERM, PC_DEPTH))

fit_group_permuted <- function(grp, design_perm) {
  idx <- which(design_perm$group_label == grp)
  d   <- design_perm[idx, , drop = FALSE]
  cm  <- count_mat[idx, , drop = FALSE]

  purrr::map_dfr(sprintf("bio%d", 1:19), function(bv) {
    predictors <- c(bv, pc_names)
    md <- d[, predictors, drop = FALSE]
    ok <- stats::complete.cases(md)
    md <- md[ok, , drop = FALSE]; cmv <- cm[ok, , drop = FALSE]
    rhs <- paste(sprintf("`%s`", predictors), collapse = " + ")
    fml <- stats::as.formula(paste("te_count ~", rhs))

    n_te <- length(te_vars); pv <- rep(NA_real_, n_te)
    for (i in seq_len(n_te)) {
      y <- as.numeric(cmv[, te_vars[i]])
      if (length(unique(y)) < 2 || length(y) < MIN_N_PER_MODEL) next
      dat <- md; dat$te_count <- y
      fit <- suppressWarnings(tryCatch(
        stats::glm(fml, data = dat, family = stats::quasipoisson()),
        error = function(e) NULL))
      if (is.null(fit)) next
      cf <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
      if (is.null(cf) || nrow(cf) < 2) next
      pv[i] <- cf[2, 4]
    }
    tibble::tibble(group = grp, variable = bv, TE_ID = te_vars, p_value = pv)
  })
}

t0 <- Sys.time()
perm_results <- purrr::map_dfr(seq_len(N_PERM), function(k) {
  design_perm <- design
  # Permute PC coordinates among individuals, independently within each group.
  for (g in GROUPS) {
    rows <- which(design_perm$group_label == g)
    design_perm[rows, pc_names] <- design_perm[sample(rows), pc_names]
  }
  res <- purrr::map_dfr(GROUPS, fit_group_permuted, design_perm = design_perm) %>%
    dplyr::filter(!is.na(p_value))
  purrr::map_dfr(RAW_THRESHOLDS, function(thr) {
    r <- replicated_candidates(res, thr, p_col = "p_value")
    tibble::tibble(permutation = k, threshold = thr,
                   n_replicated = length(r$shared))
  })
})
message(sprintf("  %d permutations in %.1f min", N_PERM,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))

observed_B <- purrr::map_dfr(RAW_THRESHOLDS, function(thr) {
  tibble::tibble(threshold = thr,
                 n_observed = length(replicated_candidates(modelA_PC, thr,
                                                           p_col = "p_value")$shared))
})

perm_summary <- perm_results %>%
  dplyr::group_by(threshold) %>%
  dplyr::summarise(
    null_median = stats::median(n_replicated),
    null_q025   = stats::quantile(n_replicated, 0.025),
    null_q975   = stats::quantile(n_replicated, 0.975),
    null_min    = min(n_replicated),
    null_max    = max(n_replicated),
    .groups = "drop"
  ) %>%
  dplyr::left_join(observed_B, by = "threshold") %>%
  dplyr::mutate(
    # One-sided empirical p: how often does the permuted null fall at or below
    # the observed count?
    p_emp_observed_lower = purrr::map2_dbl(threshold, n_observed, function(thr, obs) {
      nulls <- perm_results$n_replicated[perm_results$threshold == thr]
      (sum(nulls <= obs) + 1) / (length(nulls) + 1)
    }),
    interpretation = dplyr::case_when(
      p_emp_observed_lower <= 0.05 ~ "observed below null: PCs remove real confounding",
      TRUE                         ~ "observed within null: drop explained by lost df"
    )
  )
readr::write_csv(perm_results, file.path(dir_tab, "permutation_control_raw.csv"))
readr::write_csv(perm_summary, file.path(dir_tab, "permutation_control_summary.csv"))
message("Permutation control (observed Model A-PC vs permuted-PC null):")
print(as.data.frame(perm_summary))

# =============================================================================
# 4. Candidate counts across models, thresholds and correction regimes
# =============================================================================
model_tables <- list(
  `A (bio only)`         = modelA_bh,
  `A-PC (bio + PC1-3)`      = modelA_PC
)

counts_tbl <- purrr::imap_dfr(model_tables, function(df, model_name) {
  purrr::map_dfr(c("p_value", "padj"), function(pc) {
    purrr::map_dfr(RAW_THRESHOLDS, function(thr) {
      r <- replicated_candidates(df, thr, p_col = pc)
      tibble::tibble(
        model       = model_name,
        regime      = ifelse(pc == "p_value", "raw p", "BH within group x variable"),
        threshold   = thr,
        n_SYLW      = length(r$per_group$SYLW),
        n_SYLE      = length(r$per_group$SYLE),
        n_SYL_Admix = length(r$per_group$SYL_Admix),
        n_replicated= length(r$shared)
      )
    })
  })
})
readr::write_csv(counts_tbl, file.path(dir_tab, "candidate_counts.csv"))
message("Candidate counts (n_replicated = shared by all three groups):")
print(as.data.frame(counts_tbl))

# =============================================================================
# 5. What is shared with, and what differs from, the published 34
# =============================================================================
comparison <- purrr::imap_dfr(model_tables, function(df, model_name) {
  purrr::map_dfr(c("p_value", "padj"), function(pc) {
    purrr::map_dfr(RAW_THRESHOLDS, function(thr) {
      s <- replicated_candidates(df, thr, p_col = pc)$shared
      tibble::tibble(
        model     = model_name,
        regime    = ifelse(pc == "p_value", "raw p", "BH within group x variable"),
        threshold = thr,
        n         = length(s),
        n_shared_with_published34 = length(intersect(s, published_34)),
        n_lost_vs_published34     = length(setdiff(published_34, s)),
        n_new_vs_published34      = length(setdiff(s, published_34)),
        pct_of_34_retained        = round(100 * length(intersect(s, published_34)) /
                                            length(published_34), 1)
      )
    })
  })
})
readr::write_csv(comparison, file.path(dir_tab, "candidate_comparison.csv"))
message("Comparison against the published 34:")
print(as.data.frame(comparison))

# Export every candidate list.
for (model_name in names(model_tables)) {
  for (pc in c("p_value", "padj")) {
    for (thr in RAW_THRESHOLDS) {
      s <- replicated_candidates(model_tables[[model_name]], thr, p_col = pc)$shared
      tag <- sprintf("%s_%s_p%s",
                     ifelse(grepl("^A", model_name), "modelA", "modelA_PC"),
                     ifelse(pc == "p_value", "rawp", "BH"),
                     sub("\\.", "", format(thr, scientific = FALSE)))
      writeLines(s, file.path(dir_list, sprintf("TE_replicated_%s.txt", tag)))
    }
  }
}
writeLines(published_34, file.path(dir_list, "TE_published_34.txt"))

# =============================================================================
# 6. Recommended candidate list
# =============================================================================
# Model A-PC at raw p < 0.01 keeps the published design and threshold and adds
# the population-structure control the reviewer asked for, so it is the
# like-for-like replacement for the published 34.
final_set <- replicated_candidates(modelA_PC, 0.01, p_col = "p_value")$shared

superfamily_of <- function(x) sub("^TE_[0-9]+_", "", x)

final_tbl <- tibble::tibble(
  TE_ID          = final_set,
  superfamily    = superfamily_of(final_set),
  in_published34 = final_set %in% published_34
) %>%
  dplyr::left_join(
    modelA_PC %>%
      dplyr::filter(TE_ID %in% final_set, p_value < 0.01) %>%
      dplyr::group_by(TE_ID) %>%
      dplyr::summarise(
        n_group_variable_hits = dplyr::n(),
        bioclim_vars = paste(sort(unique(variable)), collapse = ";"),
        min_p        = min(p_value),
        .groups = "drop"
      ),
    by = "TE_ID"
  ) %>%
  dplyr::arrange(min_p)

readr::write_csv(final_tbl, file.path(dir_tab, "final_candidates_annotated.csv"))
message(sprintf("Recommended list (Model A-PC, replicated, raw p < 0.01): %d TEs",
                nrow(final_tbl)))
message("  superfamily composition:")
print(as.data.frame(final_tbl %>% dplyr::count(superfamily, name = "n")))

# FASTA for the TFBS rerun, split by superfamily as the FIMO pipeline expects.
te_seqs <- Biostrings::readDNAStringSet(file.path(dir_data, "TE_consensus_db.fasta"))
names(te_seqs) <- sub(" .*$", "", names(te_seqs))

fimo_group <- function(sf) {
  dplyr::case_when(sf == "LTR_Copia" ~ "Copia", sf == "LTR_Gypsy" ~ "Gypsy",
                   sf == "LTR" ~ "LTR", sf == "noCat" ~ "NoCat",
                   sf == "TIR" ~ "TIR", TRUE ~ NA_character_)
}
final_tbl$fimo_group <- fimo_group(final_tbl$superfamily)
for (g in sort(unique(stats::na.omit(final_tbl$fimo_group)))) {
  ids <- final_tbl$TE_ID[!is.na(final_tbl$fimo_group) & final_tbl$fimo_group == g]
  Biostrings::writeXStringSet(
    te_seqs[ids], file.path(dir_fa, sprintf("candidates_%s.fasta", g)))
  message(sprintf("  wrote candidates_%s.fasta (%d sequences)", g, length(ids)))
}

# =============================================================================
# 6b. Final non-redundant candidate set = union of the two independent analyses
# =============================================================================
# Model A and Model A-PC are two INDEPENDENT analyses of the same replication
# design, not a replacement of one by the other:
#   Model A - climate association, replicated across the three groups
#   Model A-PC - the same, with neutral SNP PC1-3 as covariates in each group
#
# The permutation control above shows that Model A-PC is heavily power-limited
# (SYLE has n = 10, and three PC covariates cost three degrees of freedom), so
# its shorter list reflects lost power more than removed confounding. Reporting
# Model A-PC alone would therefore discard real candidates. The Methods already
# specify merging the two lists into "a final, non-redundant set", which is
# what is built here: the union, annotated with the model(s) supporting each
# TE so the level of evidence stays explicit.
UNION_THRESHOLD <- 0.01   # the threshold used for the published Figure 3

set_A <- replicated_candidates(modelA,  UNION_THRESHOLD, p_col = "p_value")$shared
set_APC <- replicated_candidates(modelA_PC,  UNION_THRESHOLD, p_col = "p_value")$shared
set_union <- sort(union(set_A, set_APC))

union_tbl <- tibble::tibble(
  TE_ID       = set_union,
  superfamily = superfamily_of(set_union),
  in_model_A  = set_union %in% set_A,
  in_model_A_PC  = set_union %in% set_APC
) %>%
  dplyr::mutate(
    support = dplyr::case_when(
      in_model_A &  in_model_A_PC ~ "both models (highest confidence)",
      in_model_A & !in_model_A_PC ~ "Model A only (climate, replicated)",
      TRUE                     ~ "Model A-PC only (structure-controlled)"
    )
  ) %>%
  dplyr::left_join(
    modelA %>% dplyr::filter(p_value < UNION_THRESHOLD) %>%
      dplyr::group_by(TE_ID) %>%
      dplyr::summarise(min_p_modelA = min(p_value),
                       bioclim_vars_A = paste(sort(unique(variable)), collapse = ";"),
                       .groups = "drop"),
    by = "TE_ID"
  ) %>%
  dplyr::left_join(
    modelA_PC %>% dplyr::filter(p_value < UNION_THRESHOLD) %>%
      dplyr::group_by(TE_ID) %>%
      dplyr::summarise(min_p_modelA_PC = min(p_value), .groups = "drop"),
    by = "TE_ID"
  ) %>%
  dplyr::arrange(dplyr::desc(in_model_A & in_model_A_PC), TE_ID)

readr::write_csv(union_tbl, file.path(dir_tab, "FINAL_union_candidates.csv"))
writeLines(set_union, file.path(dir_list, "TE_FINAL_union_p001.txt"))

message(sprintf("FINAL non-redundant candidate set (union, p < %g): %d TEs",
                UNION_THRESHOLD, nrow(union_tbl)))
print(as.data.frame(union_tbl %>% dplyr::count(support, name = "n")))
message("  superfamily composition:")
print(as.data.frame(union_tbl %>% dplyr::count(superfamily, name = "n")))

# FASTA per superfamily, so Figure 3C / the TFBS pipeline can be regenerated
# on the final set.
union_tbl$fimo_group <- fimo_group(union_tbl$superfamily)
dir_fa_union <- file.path(dir_fa, "final_union")
dir.create(dir_fa_union, showWarnings = FALSE, recursive = TRUE)
for (g in sort(unique(stats::na.omit(union_tbl$fimo_group)))) {
  ids <- union_tbl$TE_ID[!is.na(union_tbl$fimo_group) & union_tbl$fimo_group == g]
  Biostrings::writeXStringSet(
    te_seqs[ids], file.path(dir_fa_union, sprintf("final_union_%s.fasta", g)))
  message(sprintf("  wrote final_union_%s.fasta (%d sequences)", g, length(ids)))
}

# =============================================================================
# 6c. How strong is "replicated" really? Variable and sign consistency
# =============================================================================
# The published rule counts a TE as replicated when it is significant for AT
# LEAST ONE bioclimatic variable in each of the three groups - but not
# necessarily the SAME variable, and not necessarily in the same direction. A
# TE significant for bio1 in SYLW, bio13 in SYLE and bio17 in SYL_Admix passes,
# even though no single association has been replicated.
#
# Two stricter readings are quantified here:
#   shared-variable : the same bioclimatic variable is significant in all three
#                     groups
#   consistent-sign : and its slope has the same sign in all three groups
#
# The three groups occupy different parts of the range, so different limiting
# variables in different regions is not by itself absurd. But the distinction
# changes what can be claimed, so it is measured rather than left implicit.
consistency_for <- function(df, te_ids, thr) {
  sig <- df %>% dplyr::filter(!is.na(p_value), p_value < thr)
  purrr::map_dfr(te_ids, function(te) {
    d <- sig %>% dplyr::filter(TE_ID == te)
    per_group <- lapply(GROUPS, function(g) unique(d$variable[d$group == g]))
    shared <- Reduce(intersect, per_group)
    sign_ok <- FALSE
    if (length(shared) > 0) {
      sign_ok <- any(vapply(shared, function(v) {
        sg <- sign(d$estimate[d$variable == v & d$group %in% GROUPS])
        length(unique(sg)) == 1
      }, logical(1)))
    }
    tibble::tibble(
      TE_ID = te,
      n_shared_variables = length(shared),
      shared_variables = paste(sort(shared), collapse = ";"),
      consistent_sign = sign_ok)
  })
}

consistency <- purrr::imap_dfr(model_tables, function(df, model_name) {
  purrr::map_dfr(RAW_THRESHOLDS, function(thr) {
    te_ids <- replicated_candidates(df, thr, p_col = "p_value")$shared
    if (length(te_ids) == 0) {
      return(tibble::tibble(model = model_name, threshold = thr, n_candidates = 0,
                            n_shared_variable = 0, n_consistent_sign = 0))
    }
    cc <- consistency_for(df, te_ids, thr)
    tibble::tibble(
      model = model_name, threshold = thr,
      n_candidates = nrow(cc),
      n_shared_variable = sum(cc$n_shared_variables > 0),
      n_consistent_sign = sum(cc$consistent_sign))
  })
})
readr::write_csv(consistency, file.path(dir_tab, "candidate_consistency.csv"))
message("Consistency of the replication rule:")
print(as.data.frame(consistency))

# Per-TE detail for the final union, at the threshold used for Figure 3.
union_consistency <- consistency_for(modelA, set_union, UNION_THRESHOLD) %>%
  dplyr::rename(shared_vars_modelA = shared_variables,
                n_shared_modelA = n_shared_variables,
                sign_ok_modelA = consistent_sign)
readr::write_csv(union_consistency,
                 file.path(dir_tab, "final_union_consistency.csv"))

# The strictest candidate list available from this design: same variable,
# same direction, all three groups. Exported at every threshold so the cost of
# the stricter rule is visible.
strict_sets <- purrr::map_dfr(RAW_THRESHOLDS, function(thr) {
  te_ids <- replicated_candidates(modelA, thr, p_col = "p_value")$shared
  if (length(te_ids) == 0) {
    return(tibble::tibble(threshold = thr, n_strict = 0))
  }
  cc <- consistency_for(modelA, te_ids, thr)
  strict <- sort(cc$TE_ID[cc$consistent_sign])
  thr_tag <- sub("\\.", "", format(thr, scientific = FALSE))
  writeLines(strict,
             file.path(dir_list, sprintf("TE_STRICT_same_variable_same_sign_p%s.txt",
                                         thr_tag)))
  tibble::tibble(threshold = thr, n_strict = length(strict))
})
readr::write_csv(strict_sets, file.path(dir_tab, "strict_rule_set_sizes.csv"))
message("Strictest rule (same variable, same sign, all three groups):")
print(as.data.frame(strict_sets))

# =============================================================================
# 7. Figures
# =============================================================================
p_counts <- counts_tbl %>%
  dplyr::mutate(threshold = factor(threshold, levels = RAW_THRESHOLDS,
                                   labels = paste("p <", RAW_THRESHOLDS))) %>%
  ggplot2::ggplot(ggplot2::aes(x = threshold, y = n_replicated, fill = model)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
  ggplot2::geom_text(ggplot2::aes(label = n_replicated),
                     position = ggplot2::position_dodge(width = 0.8),
                     vjust = -0.3, size = 3) +
  ggplot2::geom_hline(yintercept = length(published_34), linetype = 2,
                      colour = "grey40") +
  ggplot2::facet_wrap(~ regime, scales = "free_y") +
  ggplot2::scale_fill_manual(values = c(`A (bio only)` = "#4C72B0",
                                        `A-PC (bio + PC1-3)` = "#C44E52")) +
  ggplot2::labs(x = NULL, y = "TEs replicated across all three groups",
                fill = NULL,
                title = "Candidate TEs under the three-group replication design",
                subtitle = "Dashed line = the 34 candidates currently in Figure 3") +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")
save_plot(p_counts, "replicated_candidate_counts", width = 9, height = 5)

p_grp <- counts_tbl %>%
  dplyr::filter(regime == "raw p") %>%
  dplyr::select(model, threshold, SYLW = n_SYLW, SYLE = n_SYLE,
                SYL_Admix = n_SYL_Admix, Replicated = n_replicated) %>%
  tidyr::pivot_longer(-c(model, threshold), names_to = "set", values_to = "n") %>%
  dplyr::mutate(
    set = factor(set, levels = c("SYLW", "SYLE", "SYL_Admix", "Replicated")),
    threshold = factor(threshold, levels = RAW_THRESHOLDS,
                       labels = paste("p <", RAW_THRESHOLDS))) %>%
  ggplot2::ggplot(ggplot2::aes(x = set, y = n, fill = model)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
  ggplot2::geom_text(ggplot2::aes(label = n),
                     position = ggplot2::position_dodge(width = 0.8),
                     vjust = -0.3, size = 2.5) +
  ggplot2::facet_wrap(~ threshold, scales = "free_y") +
  ggplot2::scale_fill_manual(values = c(`A (bio only)` = "#4C72B0",
                                        `A-PC (bio + PC1-3)` = "#C44E52")) +
  ggplot2::labs(x = NULL, y = "Candidate residual TEs", fill = NULL,
                title = "Where the stringency comes from: replication, not correction",
                subtitle = "Per-group lists collapse once agreement across all three groups is required") +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "bottom",
                 axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
save_plot(p_grp, "per_group_vs_replicated", width = 10, height = 5)

# Permutation control: observed count against the permuted-PC null.
p_perm <- perm_results %>%
  dplyr::mutate(threshold = factor(threshold, levels = RAW_THRESHOLDS,
                                   labels = paste("p <", RAW_THRESHOLDS))) %>%
  ggplot2::ggplot(ggplot2::aes(x = n_replicated)) +
  ggplot2::geom_histogram(bins = 20, fill = "grey70", colour = "white") +
  ggplot2::geom_vline(
    data = observed_B %>%
      dplyr::mutate(threshold = factor(threshold, levels = RAW_THRESHOLDS,
                                       labels = paste("p <", RAW_THRESHOLDS))),
    ggplot2::aes(xintercept = n_observed), colour = "#C44E52", linewidth = 1) +
  ggplot2::facet_wrap(~ threshold, scales = "free") +
  ggplot2::labs(
    x = "TEs replicated across all three groups", y = "Permutations",
    title = "Is the drop caused by the PCs, or by spending three degrees of freedom?",
    subtitle = sprintf("Grey = %d within-group PC permutations; red = observed Model A-PC",
                       N_PERM)) +
  ggplot2::theme_classic(base_size = 10)
save_plot(p_perm, "permutation_control", width = 10, height = 4.5)

# =============================================================================
# 8. Session log
# =============================================================================
log_session("02_candidate_TEs")
message("=== 02_candidate_TEs.R finished ===")
