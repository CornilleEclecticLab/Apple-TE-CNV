# =============================================================================
# Script:  05_tfbs_enrichment.R
# Purpose: Reviewer 2, comment 6 (R2.6) - supply the missing background for
#          the TFBS enrichment analysis, and regenerate panel Fig. 3C on the
#          final candidate set.
#
#          The published analysis ran FIMO (PlantTFDB v5.0 Malus domestica
#          motifs) on the 34 climate-associated residual TEs, split by TE
#          superfamily, scored enrichment as the median -log10(q-value) per
#          TE x TF-family pair, and called a pair enriched above 1.3. There
#          was no background group, so the reviewer is right that the result
#          may only show that certain TE families contain TFBSs at all,
#          independently of climate.
#
#          Two things changed since the first version of this script. The
#          candidate set is now the final non-redundant union of 38 TEs from
#          02_candidate_TEs.R, not the 34 and not the
#          723 that the COMBINED-only model produced. And "non-climate-
#          associated" is now defined symmetrically with how the candidates
#          were defined: a background TE must not reach p < 0.01 for ANY
#          bioclimatic variable in ANY of the three M. sylvestris groups,
#          under either Model A or Model A-PC.
#
#            PART 1 (runs now)
#              - assembles the 38-TE foreground and writes its FASTAs
#              - draws a superfamily-matched background from the 550 TEs that
#                show no climate association anywhere
#              - writes a ready-to-run FIMO driver covering BOTH sets, so
#                foreground and background are scanned in the same run with
#                the same parameters
#              - reports the legacy 34-TE foreground for reference
#
#            PART 2 (runs automatically once the FIMO output exists)
#              - applies the identical enrichment pipeline to both sets
#              - compares them per TF family and overall, controlling for
#                sequence length
#              - verifies or replaces the TF families named in the Discussion
#
# NAMING NOTE: in the manuscript, "Model A" is the climate GLM
# (TE count ~ bioclimatic variable) and "Model B" is the SPATIAL model
# (TE count ~ latitude + longitude, negative binomial), which produces the
# vector-strength (VS) values shown in Fig. 3A-B. To avoid a clash with that
# established terminology, the PC-covariate version of the climate model is
# called "Model A-PC" throughout this folder, never "Model B".
#
# Author:  Antho (analysis run with Claude Code)
# Date:    2026-08-24
# Project: TE_Project / review (reviewer response, Genome Biology submission)
#
# Inputs:
#   outputs/02_candidate_TEs/tables/
#       FINAL_union_candidates.csv          the 38 final candidates
#       glm_modelA_PC_by_group.csv          Model A-PC, per group
#   data/modelA_fig3_raw_p.csv              Model A, per group
#   data/Mdo_TF_list.txt                    motif_id -> TF family
#   data/Mdo_TF_binding_motifs.meme         PlantTFDB v5.0 motifs
#   data/TE_consensus_db.fasta              TE consensus sequences
#   data/fimo_modelA_outputs/fimo_*.tsv     legacy FIMO output for the 34
#   outputs/05_tfbs_enrichment/fimo_out/{foreground,background}/
#       fimo_*.tsv                          PART 2 only - from the FIMO driver
#
# Outputs (all under outputs/05_tfbs_enrichment/):
#   fasta/foreground/*.fasta  fasta/background/*.fasta
#       per superfamily, plus all_candidates_38.fasta / all_background_304.fasta
#   run_fimo.sh / run_fimo.slurm
#   tables/background_TE_selection.csv      both sets, with sequence lengths
#   tables/background_availability.csv      matching feasibility
#   tables/legacy34_TF_family_summary.csv   published foreground, for reference
#   tables/enrichment_comparison_by_TF_family.csv   (PART 2)
#   tables/enrichment_overall_tests.csv             (PART 2)
#   tables/manuscript_TF_family_verdict.csv         (PART 2)
#   figures/*.pdf + *.svg
#   ../../logs/05_tfbs_enrichment_*_sessionInfo.txt
# =============================================================================

# --- Setup -------------------------------------------------------------------
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(Biostrings)
  library(patchwork)
  library(svglite)
})

# The seed is set in 00_setup.R; restated because the background draw below is
# the one stochastic step here and must be reproducible on its own.
set.seed(20260824)

# --- Parameters --------------------------------------------------------------
BACKGROUND_MULTIPLIER <- 10     # desired background size; capped if infeasible
CLIMATE_EXCL_P        <- 0.01   # a background TE must not reach this anywhere
Q_THRESHOLD           <- 0.05   # FIMO q-value cutoff, as in the published run
ENRICH_THRESHOLD      <- 1.3    # median -log10(q) above which a pair is enriched

dir_out    <- file.path(dir_outputs, "05_tfbs_enrichment")
dir_tab    <- file.path(dir_out, "tables")
dir_fig    <- file.path(dir_out, "figures")
dir_fa_fg  <- file.path(dir_out, "fasta", "foreground")
dir_fa_bg  <- file.path(dir_out, "fasta", "background")
dir_fimo   <- file.path(dir_out, "fimo_out")
for (d in c(dir_out, dir_tab, dir_fig, dir_fa_fg, dir_fa_bg,
            file.path(dir_fimo, "foreground"), file.path(dir_fimo, "background"))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

save_plot <- function(p, name, width = 7, height = 5) {
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".pdf")), p,
                  width = width, height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".svg")), p,
                  width = width, height = height, device = svglite::svglite)
}

message("=== 05_tfbs_enrichment.R ===")

# =============================================================================
# 1. Shared definitions
# =============================================================================
# FIMO was run separately per TE superfamily; these are the five groups.
FIMO_GROUP_OF <- function(te_id) {
  suffix <- sub("^TE_[0-9]+_", "", te_id)
  dplyr::case_when(
    suffix == "LTR_Copia" ~ "Copia",
    suffix == "LTR_Gypsy" ~ "Gypsy",
    suffix == "LTR"       ~ "LTR",
    suffix == "noCat"     ~ "NoCat",
    suffix == "TIR"       ~ "TIR",
    TRUE                  ~ NA_character_
  )
}

tf_family <- readr::read_tsv(file.path(dir_data, "Mdo_TF_list.txt"),
                             show_col_types = FALSE)

#' Apply the published enrichment pipeline to a set of FIMO output files.
#' Keeps hits with q < Q_THRESHOLD, scores them as -log10(q), attaches the TF
#' family, drops unannotated motifs, then takes the median score per
#' TF family x TE pair - identical to the published code.
run_enrichment_pipeline <- function(files) {
  files <- files[file.exists(files)]
  if (length(files) == 0) return(NULL)

  raw <- purrr::imap_dfr(files, function(path, sf) {
    readr::read_tsv(path, show_col_types = FALSE, comment = "#",
                    progress = FALSE) %>%
      dplyr::mutate(TE_type = sf)
  })

  hits <- raw %>%
    dplyr::filter(!is.na(`q-value`), `q-value` < Q_THRESHOLD) %>%
    dplyr::mutate(neglog10_q = -log10(`q-value`)) %>%
    dplyr::left_join(tf_family, by = c("motif_id" = "TF_ID")) %>%
    dplyr::filter(!is.na(Family))

  pairs <- hits %>%
    dplyr::group_by(Family, sequence_name, TE_type) %>%
    dplyr::summarise(median_neglog10_q = stats::median(neglog10_q),
                     n_motif_hits = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(enriched = median_neglog10_q > ENRICH_THRESHOLD)

  list(hits = hits, pairs = pairs)
}

# =============================================================================
# 2. PART 1a - the 38-TE foreground
# =============================================================================
final_union <- readr::read_csv(
  file.path(dir_outputs, "02_candidate_TEs", "tables",
            "FINAL_union_candidates.csv"),
  show_col_types = FALSE
)
fg_te_ids <- sort(unique(final_union$TE_ID))

fg_composition <- tibble::tibble(TE_ID = fg_te_ids,
                                 fimo_group = FIMO_GROUP_OF(fg_te_ids)) %>%
  dplyr::filter(!is.na(fimo_group)) %>%
  dplyr::count(fimo_group, name = "n_foreground")

message(sprintf("Foreground: %d final candidate TEs", length(fg_te_ids)))
print(as.data.frame(fg_composition))

# =============================================================================
# 3. PART 1b - background pool and superfamily-matched draw
# =============================================================================
# Symmetric with the candidate definition: a background TE must not reach
# p < CLIMATE_EXCL_P for any bioclimatic variable, in any of the three groups,
# under either the climate model or its PC-covariate version.
modelA <- readr::read_csv(file.path(dir_data, "modelA_fig3_raw_p.csv"),
                          show_col_types = FALSE)
modelA_PC <- readr::read_csv(
  file.path(dir_outputs, "02_candidate_TEs", "tables",
            "glm_modelA_PC_by_group.csv"),
  show_col_types = FALSE
)

te_universe <- sort(unique(modelA$TE_ID))
climate_any <- union(
  unique(modelA$TE_ID[!is.na(modelA$p_value) & modelA$p_value < CLIMATE_EXCL_P]),
  unique(modelA_PC$TE_ID[!is.na(modelA_PC$p_value) &
                           modelA_PC$p_value < CLIMATE_EXCL_P])
)
background_pool <- setdiff(te_universe, union(climate_any, fg_te_ids))

message(sprintf("TE universe: %d | climate-associated anywhere (p < %g): %d | background pool: %d",
                length(te_universe), CLIMATE_EXCL_P, length(climate_any),
                length(background_pool)))

pool_tbl <- tibble::tibble(TE_ID = background_pool,
                           fimo_group = FIMO_GROUP_OF(background_pool)) %>%
  dplyr::filter(!is.na(fimo_group))

# The requested multiplier is capped at the largest value the pool can supply
# while preserving the foreground superfamily composition exactly, rather than
# failing or silently distorting the matching.
availability <- fg_composition %>%
  dplyr::left_join(pool_tbl %>% dplyr::count(fimo_group, name = "n_available"),
                   by = "fimo_group") %>%
  dplyr::mutate(n_available = tidyr::replace_na(n_available, 0L),
                max_multiplier = floor(n_available / n_foreground))

feasible_multiplier <- min(BACKGROUND_MULTIPLIER, min(availability$max_multiplier))
availability <- availability %>%
  dplyr::mutate(multiplier_used = feasible_multiplier,
                n_required = n_foreground * feasible_multiplier)
readr::write_csv(availability, file.path(dir_tab, "background_availability.csv"))

message("Background availability by superfamily:")
print(as.data.frame(availability))
if (feasible_multiplier < BACKGROUND_MULTIPLIER) {
  message(sprintf("NOTE: requested %dx is not available; using %dx (%d background TEs). ",
                  BACKGROUND_MULTIPLIER, feasible_multiplier,
                  feasible_multiplier * length(fg_te_ids)),
          "Limiting superfamily: ",
          paste(availability$fimo_group[availability$max_multiplier ==
                                          min(availability$max_multiplier)],
                collapse = ", "))
}

background_te <- purrr::map_dfr(seq_len(nrow(availability)), function(i) {
  g <- availability$fimo_group[i]
  k <- availability$n_required[i]
  candidates <- sort(pool_tbl$TE_ID[pool_tbl$fimo_group == g])
  tibble::tibble(TE_ID = sample(candidates, size = k, replace = FALSE),
                 fimo_group = g)
}) %>%
  dplyr::arrange(fimo_group, TE_ID)

message(sprintf("Background selected: %d TEs (%dx the %d foreground TEs)",
                nrow(background_te), feasible_multiplier, length(fg_te_ids)))

# =============================================================================
# 4. PART 1c - sequences, length balance, FASTA export
# =============================================================================
te_seqs <- Biostrings::readDNAStringSet(file.path(dir_data, "TE_consensus_db.fasta"))
names(te_seqs) <- sub(" .*$", "", names(te_seqs))

missing_seq <- setdiff(c(fg_te_ids, background_te$TE_ID), names(te_seqs))
if (length(missing_seq) > 0) {
  stop("TEs absent from the consensus FASTA: ",
       paste(utils::head(missing_seq, 10), collapse = ", "))
}

# A longer consensus gives FIMO more chances to score a good hit, so sequence
# length is recorded and tested rather than assumed away.
selection <- dplyr::bind_rows(
  tibble::tibble(TE_ID = fg_te_ids, fimo_group = FIMO_GROUP_OF(fg_te_ids),
                 group = "climate_associated"),
  background_te %>% dplyr::mutate(group = "background")
) %>%
  dplyr::filter(!is.na(fimo_group)) %>%
  dplyr::mutate(seq_length = as.integer(BiocGenerics::width(te_seqs[TE_ID])))

readr::write_csv(selection, file.path(dir_tab, "background_TE_selection.csv"))

len_test <- stats::wilcox.test(seq_length ~ group, data = selection, exact = FALSE)
len_summary <- selection %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(n = dplyr::n(), median_length = stats::median(seq_length),
                   mean_length = mean(seq_length), .groups = "drop") %>%
  dplyr::mutate(wilcoxon_p_length = len_test$p.value)
readr::write_csv(len_summary, file.path(dir_tab, "sequence_length_balance.csv"))
message("Sequence length balance between groups:")
print(as.data.frame(len_summary))

write_set_fasta <- function(ids_tbl, outdir, prefix) {
  for (g in sort(unique(ids_tbl$fimo_group))) {
    ids <- ids_tbl$TE_ID[ids_tbl$fimo_group == g]
    out <- file.path(outdir, sprintf("%s_%s.fasta", prefix, g))
    Biostrings::writeXStringSet(te_seqs[ids], filepath = out)
    message(sprintf("  wrote %s (%d sequences)", basename(out), length(ids)))
  }
}
write_set_fasta(selection %>% dplyr::filter(group == "climate_associated"),
                dir_fa_fg, "foreground")
write_set_fasta(selection %>% dplyr::filter(group == "background"),
                dir_fa_bg, "background")

# Also write each set as a single combined FASTA. FIMO is run per superfamily
# in the driver below, matching the published pipeline, but a single file is
# handier for a one-off scan or for depositing with the manuscript.
write_combined <- function(ids, path) {
  Biostrings::writeXStringSet(te_seqs[ids], filepath = path)
  message(sprintf("  wrote %s (%d sequences)", basename(path), length(ids)))
}
write_combined(selection$TE_ID[selection$group == "climate_associated"],
               file.path(dir_fa_fg, "all_candidates_38.fasta"))
write_combined(selection$TE_ID[selection$group == "background"],
               file.path(dir_fa_bg, "all_background_304.fasta"))

# =============================================================================
# 5. PART 1d - FIMO driver covering both sets
# =============================================================================
# Foreground and background are scanned in the same run, with the same motif
# file and FIMO defaults, so the comparison cannot be attributed to a
# difference in how the two sets were processed.
fimo_sh <- c(
  "#!/usr/bin/env bash",
  "# ---------------------------------------------------------------------",
  "# FIMO run for reviewer comment R2.6 and Fig. 3C.",
  "# Generated by scripts/05_tfbs_enrichment.R - do not edit by hand.",
  "#",
  "# Scans BOTH the final candidate TEs and the superfamily-matched",
  "# background, with DEFAULT FIMO parameters and the same PlantTFDB v5.0",
  "# M. domestica motif file used for the published analysis.",
  "#",
  "# Requires the MEME Suite on PATH (fimo --version should work).",
  "# Run from the review/ project root:",
  "#   bash outputs/05_tfbs_enrichment/run_fimo.sh",
  "# Then re-run: Rscript scripts/05_tfbs_enrichment.R",
  "# ---------------------------------------------------------------------",
  "set -euo pipefail",
  "",
  'HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
  'MOTIFS="$HERE/../../data/Mdo_TF_binding_motifs.meme"',
  "",
  '# MEME Suite 5.5.7 was built locally into $HOME/opt/meme; fall back to it',
  '# so the driver works without editing any shell profile.',
  'if ! command -v fimo >/dev/null 2>&1; then',
  '  if [ -x "$HOME/opt/meme/bin/fimo" ]; then',
  '    export PATH="$HOME/opt/meme/bin:$PATH"',
  '  else',
  '    echo "ERROR: fimo not found on PATH and not in $HOME/opt/meme/bin."; exit 1',
  '  fi',
  'fi',
  'echo "using $(command -v fimo) ($(fimo --version))"',
  "",
  "for SET in foreground background; do",
  '  for SF in Copia Gypsy LTR NoCat TIR; do',
  '    FA="$HERE/fasta/${SET}/${SET}_${SF}.fasta"',
  '    [ -f "$FA" ] || { echo "skip ${SET} ${SF} (no FASTA)"; continue; }',
  '    echo ">>> FIMO on ${SET} ${SF} ..."',
  '    fimo --oc "$HERE/fimo_out/${SET}/raw_${SF}" "$MOTIFS" "$FA"',
  '    cp "$HERE/fimo_out/${SET}/raw_${SF}/fimo.tsv" "$HERE/fimo_out/${SET}/fimo_${SF}.tsv"',
  "  done",
  "done",
  "",
  'echo "Done. Now re-run: Rscript scripts/05_tfbs_enrichment.R"'
)
sh_path <- file.path(dir_out, "run_fimo.sh")
writeLines(fimo_sh, sh_path)
Sys.chmod(sh_path, mode = "0755")

writeLines(c(
  "#!/bin/bash",
  "#SBATCH --job-name=fimo_R2.6",
  "#SBATCH --time=06:00:00",
  "#SBATCH --cpus-per-task=4",
  "#SBATCH --mem=16G",
  "#SBATCH --output=fimo_R26_%j.out",
  "# ---------------------------------------------------------------------",
  "# Slurm wrapper for the FIMO run (R2.6 / Fig. 3C), for jubail.",
  "# Adjust the module line to whatever the cluster provides, then:",
  "#   sbatch run_fimo.slurm",
  "# Copy outputs/05_tfbs_enrichment/fimo_out/ back to the local",
  "# review/ folder before re-running script 03.",
  "# ---------------------------------------------------------------------",
  "set -euo pipefail",
  "",
  "module load meme || echo 'adjust this module load line for your cluster'",
  "",
  'bash "$(dirname "$0")/run_fimo.sh"'
), file.path(dir_out, "run_fimo.slurm"))

message("FIMO driver written: ", sh_path)

# =============================================================================
# 6. PART 1e - legacy 34-TE foreground, for reference
# =============================================================================
# The published FIMO output still exists and covers the 34 candidates that
# were in Fig. 3C. It is summarised here so the TF families currently named
# in the Discussion can be checked against it even before the new run, and so
# the old and new foregrounds can be compared afterwards.
legacy_files <- c(
  Copia = file.path(dir_data, "fimo_modelA_outputs", "fimo_copia.tsv"),
  Gypsy = file.path(dir_data, "fimo_modelA_outputs", "fimo_gypsy.tsv"),
  LTR   = file.path(dir_data, "fimo_modelA_outputs", "fimo_LTR.tsv"),
  NoCat = file.path(dir_data, "fimo_modelA_outputs", "fimo_nocat.tsv"),
  TIR   = file.path(dir_data, "fimo_modelA_outputs", "fimo_TIR.tsv")
)
legacy <- run_enrichment_pipeline(legacy_files)
if (!is.null(legacy)) {
  legacy_ids <- unique(legacy$hits$sequence_name)
  legacy_summary <- legacy$pairs %>%
    dplyr::group_by(Family) %>%
    dplyr::summarise(n_TE_with_pair = dplyr::n(),
                     n_TE_enriched = sum(enriched),
                     prop_TE_enriched = sum(enriched) / length(legacy_ids),
                     median_of_medians = stats::median(median_neglog10_q),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_TE_enriched))
  readr::write_csv(legacy_summary,
                   file.path(dir_tab, "legacy34_TF_family_summary.csv"))
  message(sprintf("Legacy foreground (published Fig. 3C): %d TEs, %d TF families",
                  length(legacy_ids), nrow(legacy_summary)))
}

# =============================================================================
# 7. PART 2 - the comparison, once FIMO output exists
# =============================================================================
sf_levels <- sort(unique(selection$fimo_group))
fg_files <- stats::setNames(
  file.path(dir_fimo, "foreground", sprintf("fimo_%s.tsv", sf_levels)), sf_levels)
bg_files <- stats::setNames(
  file.path(dir_fimo, "background", sprintf("fimo_%s.tsv", sf_levels)), sf_levels)

if (!all(file.exists(c(fg_files, bg_files)))) {
  missing <- c(fg_files, bg_files)[!file.exists(c(fg_files, bg_files))]
  message("")
  message("---------------------------------------------------------------")
  message("PART 1 complete. FIMO output not found, so PART 2 (the ")
  message("statistical comparison) has been skipped.")
  message("")
  message("Next step:")
  message("    bash ", sh_path)
  message("then re-run this script; PART 2 will execute automatically.")
  message("Missing: ", paste(basename(missing), collapse = ", "))
  message("---------------------------------------------------------------")
  log_session("05_tfbs_enrichment")
  message("=== 05_tfbs_enrichment.R finished (PART 1 only) ===")
  quit(save = "no", status = 0)
}

message("FIMO output found - running PART 2 ...")

fg <- run_enrichment_pipeline(fg_files)
bg <- run_enrichment_pipeline(bg_files)

n_fg_total <- sum(selection$group == "climate_associated")
n_bg_total <- sum(selection$group == "background")

message(sprintf("Foreground: %d TEs with hits / %d; background: %d / %d",
                dplyr::n_distinct(fg$hits$sequence_name), n_fg_total,
                dplyr::n_distinct(bg$hits$sequence_name), n_bg_total))

# --- 7a. Per-TF-family comparison -------------------------------------------
all_families <- sort(union(unique(fg$pairs$Family), unique(bg$pairs$Family)))
pairs_all <- dplyr::bind_rows(
  fg$pairs %>% dplyr::mutate(group = "climate_associated"),
  bg$pairs %>% dplyr::mutate(group = "background")
)

comparison <- purrr::map_dfr(all_families, function(fam) {
  f_en <- sum(fg$pairs$Family == fam & fg$pairs$enriched)
  b_en <- sum(bg$pairs$Family == fam & bg$pairs$enriched)

  # Fisher's exact test: enriched vs not, foreground vs background. TEs with no
  # significant hit still count in the denominators.
  ct <- matrix(c(f_en, n_fg_total - f_en, b_en, n_bg_total - b_en), nrow = 2)
  ft <- stats::fisher.test(ct)

  f_sc <- fg$pairs$median_neglog10_q[fg$pairs$Family == fam]
  b_sc <- bg$pairs$median_neglog10_q[bg$pairs$Family == fam]
  wt_p <- if (length(f_sc) >= 3 && length(b_sc) >= 3) {
    stats::wilcox.test(f_sc, b_sc, exact = FALSE)$p.value
  } else NA_real_

  tibble::tibble(
    Family = fam,
    n_climate_enriched = f_en, n_climate_total = n_fg_total,
    prop_climate_enriched = f_en / n_fg_total,
    n_bg_enriched = b_en, n_bg_total = n_bg_total,
    prop_bg_enriched = b_en / n_bg_total,
    odds_ratio = unname(ft$estimate), fisher_p = ft$p.value,
    median_score_climate = if (length(f_sc)) stats::median(f_sc) else NA_real_,
    median_score_bg = if (length(b_sc)) stats::median(b_sc) else NA_real_,
    wilcoxon_p_score = wt_p
  )
}) %>%
  dplyr::mutate(fisher_padj = stats::p.adjust(fisher_p, method = "BH"),
                wilcoxon_padj = stats::p.adjust(wilcoxon_p_score, method = "BH")) %>%
  dplyr::arrange(fisher_padj, dplyr::desc(odds_ratio))

readr::write_csv(comparison,
                 file.path(dir_tab, "enrichment_comparison_by_TF_family.csv"))

sig_families <- comparison %>%
  dplyr::filter(fisher_padj < 0.05,
                prop_climate_enriched > prop_bg_enriched) %>%
  dplyr::pull(Family)
message(sprintf("TF families enriched in candidates vs background (BH < 0.05): %s",
                if (length(sig_families)) paste(sig_families, collapse = ", ") else "none"))

# --- 7b. Overall tests, controlling for sequence length ----------------------
per_te <- selection %>%
  dplyr::select(TE_ID, group, seq_length, fimo_group) %>%
  dplyr::left_join(
    pairs_all %>%
      dplyr::group_by(sequence_name) %>%
      dplyr::summarise(n_families_enriched = sum(enriched),
                       max_score = max(median_neglog10_q), .groups = "drop"),
    by = c("TE_ID" = "sequence_name")
  ) %>%
  dplyr::mutate(n_families_enriched = tidyr::replace_na(n_families_enriched, 0),
                max_score = tidyr::replace_na(max_score, 0),
                group = factor(group,
                               levels = c("background", "climate_associated")))

w_fam <- stats::wilcox.test(n_families_enriched ~ group, data = per_te, exact = FALSE)
w_sc  <- stats::wilcox.test(max_score ~ group, data = per_te, exact = FALSE)
m_len <- stats::glm(n_families_enriched ~ group + log(seq_length),
                    data = per_te, family = stats::quasipoisson())
m_coef <- summary(m_len)$coefficients

overall <- tibble::tibble(
  test = c("Wilcoxon: enriched TF families per TE",
           "Wilcoxon: max median -log10(q) per TE",
           "quasi-Poisson GLM: group effect adjusted for log(sequence length)"),
  statistic = c(unname(w_fam$statistic), unname(w_sc$statistic),
                unname(m_coef["groupclimate_associated", 3])),
  p_value = c(w_fam$p.value, w_sc$p.value,
              unname(m_coef["groupclimate_associated", 4])),
  effect = c(
    sprintf("median %0.1f (candidates) vs %0.1f (background)",
            stats::median(per_te$n_families_enriched[per_te$group == "climate_associated"]),
            stats::median(per_te$n_families_enriched[per_te$group == "background"])),
    sprintf("median %0.2f (candidates) vs %0.2f (background)",
            stats::median(per_te$max_score[per_te$group == "climate_associated"]),
            stats::median(per_te$max_score[per_te$group == "background"])),
    sprintf("log rate ratio = %0.3f", unname(m_coef["groupclimate_associated", 1]))
  )
)
readr::write_csv(overall, file.path(dir_tab, "enrichment_overall_tests.csv"))
readr::write_csv(per_te, file.path(dir_tab, "per_TE_enrichment_summary.csv"))
message("Overall comparison:")
print(as.data.frame(overall))

# --- 7c. Verdict on the TF families named in the Discussion ------------------
manuscript_families <- c("C2H2", "MYB", "NAC", "WRKY", "ERF", "bZIP",
                         "MIKC_MADS", "M-type_MADS", "BBR-BPC")

verdict <- comparison %>%
  dplyr::filter(Family %in% manuscript_families) %>%
  dplyr::mutate(named_in_manuscript = TRUE,
                verdict = dplyr::case_when(
                  fisher_padj < 0.05 & prop_climate_enriched > prop_bg_enriched ~ "SUPPORTED",
                  fisher_padj < 0.05 & prop_climate_enriched < prop_bg_enriched ~ "REVERSED (depleted)",
                  TRUE ~ "NOT SUPPORTED"))

missing_from_ms <- comparison %>%
  dplyr::filter(fisher_padj < 0.05, prop_climate_enriched > prop_bg_enriched,
                !Family %in% manuscript_families) %>%
  dplyr::mutate(named_in_manuscript = FALSE, verdict = "SHOULD BE ADDED")

verdict_full <- dplyr::bind_rows(verdict, missing_from_ms) %>%
  dplyr::select(Family, named_in_manuscript, verdict, prop_climate_enriched,
                prop_bg_enriched, odds_ratio, fisher_p, fisher_padj) %>%
  dplyr::arrange(dplyr::desc(named_in_manuscript), fisher_padj)
readr::write_csv(verdict_full,
                 file.path(dir_tab, "manuscript_TF_family_verdict.csv"))
message("Verdict on the TF families named in the Discussion:")
print(as.data.frame(verdict_full))

# =============================================================================
# 8. Figures
# =============================================================================
p_vol <- comparison %>%
  dplyr::mutate(log2_or = log2(pmax(pmin(odds_ratio, 1e3), 1e-3)),
                signif = fisher_padj < 0.05,
                label = ifelse(fisher_padj < 0.05, Family, NA_character_)) %>%
  ggplot2::ggplot(ggplot2::aes(x = log2_or, y = -log10(fisher_padj))) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  ggplot2::geom_point(ggplot2::aes(colour = signif), size = 2, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "grey65")) +
  ggplot2::geom_text(ggplot2::aes(label = label), size = 2.6, vjust = -0.8,
                     na.rm = TRUE) +
  ggplot2::labs(x = "log2 odds ratio (candidates vs background)",
                y = "-log10 BH-adjusted Fisher p",
                title = "TF-family binding-site enrichment against a matched background",
                subtitle = sprintf("%d candidate TEs vs %d non-climate-associated TEs",
                                   n_fg_total, n_bg_total)) +
  ggplot2::theme_classic(base_size = 10)
save_plot(p_vol, "TF_family_enrichment_volcano", width = 8, height = 6)

top_fams <- comparison %>% dplyr::filter(fisher_padj < 0.05) %>%
  dplyr::arrange(fisher_padj) %>% utils::head(25) %>% dplyr::pull(Family)
if (length(top_fams) > 0) {
  p_prop <- comparison %>%
    dplyr::filter(Family %in% top_fams) %>%
    dplyr::select(Family, Candidates = prop_climate_enriched,
                  Background = prop_bg_enriched) %>%
    tidyr::pivot_longer(-Family, names_to = "group", values_to = "prop") %>%
    ggplot2::ggplot(ggplot2::aes(x = stats::reorder(Family, prop), y = prop,
                                 fill = group)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(Candidates = "#C44E52",
                                          Background = "#4C72B0")) +
    ggplot2::labs(x = NULL, y = "Proportion of TEs with an enriched binding site",
                  fill = NULL,
                  title = "TF families differing between candidates and background",
                  subtitle = "Fisher's exact test, BH-adjusted p < 0.05") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(legend.position = "bottom")
  save_plot(p_prop, "TF_family_enrichment_proportions", width = 8, height = 7)
}

# Fig. 3C regenerated on the final candidate set.
p_heat <- fg$pairs %>%
  ggplot2::ggplot(ggplot2::aes(x = sequence_name,
                               y = forcats::fct_reorder(Family, median_neglog10_q,
                                                        max, na.rm = TRUE))) +
  ggplot2::geom_tile(ggplot2::aes(fill = median_neglog10_q), colour = "black",
                     linewidth = 0.1) +
  ggplot2::scale_fill_viridis_c(name = "Median -log10(q)") +
  ggplot2::facet_grid(. ~ TE_type, scales = "free_x", space = "free_x") +
  ggplot2::labs(x = "TE", y = "TF family",
                title = "Fig. 3C regenerated on the final candidate set") +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,
                                                     hjust = 1, size = 6),
                 axis.text.y = ggplot2::element_text(size = 7))
save_plot(p_heat, "Fig3C_heatmap_final_candidates", width = 12, height = 7)

p_te <- per_te %>%
  ggplot2::ggplot(ggplot2::aes(x = group, y = n_families_enriched, fill = group)) +
  ggplot2::geom_boxplot(outlier.size = 0.7, width = 0.6, show.legend = FALSE) +
  ggplot2::scale_fill_manual(values = c(background = "#4C72B0",
                                        climate_associated = "#C44E52")) +
  ggplot2::labs(x = NULL, y = "Enriched TF families per TE",
                title = "Breadth of TFBS enrichment per TE",
                subtitle = sprintf("Wilcoxon p = %.3g", w_fam$p.value)) +
  ggplot2::theme_classic(base_size = 10)
save_plot(p_te, "per_TE_enrichment_breadth", width = 5, height = 5)

# =============================================================================
# 8b. Figure 3C, redrawn
# =============================================================================
# The published panel was a 29 x 34 tile grid of median -log10(q) with 6 pt
# labels: it shows which motifs are present, but presence is exactly what the
# reviewer objected to, since a motif can be common in every TE regardless of
# climate. The redraw puts the COMPARISON first and keeps the per-TE detail as
# support, restricted to the families that actually differ.

# --- 8b-i. Headline: candidates against the matched background --------------
# A dumbbell reads the paired comparison directly - one row per TF family, the
# two proportions joined - which a heatmap cannot show at all.
dumb <- comparison %>%
  dplyr::filter(prop_climate_enriched > 0 | prop_bg_enriched > 0) %>%
  dplyr::mutate(
    signif = fisher_padj < 0.05,
    direction = dplyr::case_when(
      !signif ~ "not significant",
      prop_climate_enriched > prop_bg_enriched ~ "enriched in candidates",
      TRUE ~ "depleted in candidates"),
    # an infinite odds ratio means the family is absent from the background;
    # printing "999" would read as a real estimate
    or_txt = dplyr::case_when(
      !is.finite(odds_ratio) ~ "OR inf",
      odds_ratio >= 100      ~ "OR >100",
      TRUE                   ~ sprintf("OR %.1f", odds_ratio)),
    label = ifelse(signif, paste0(Family, "  (", or_txt, ")"), Family)) %>%
  dplyr::arrange(prop_climate_enriched - prop_bg_enriched)
dumb$label <- factor(dumb$label, levels = dumb$label)

p_dumb <- ggplot2::ggplot(dumb) +
  ggplot2::geom_segment(
    ggplot2::aes(y = label, yend = label,
                 x = prop_bg_enriched * 100, xend = prop_climate_enriched * 100,
                 colour = direction),
    linewidth = 1.1, lineend = "round") +
  ggplot2::geom_point(ggplot2::aes(y = label, x = prop_bg_enriched * 100,
                                   shape = "Background"), size = 2.4,
                      colour = "grey35", fill = "white", stroke = 0.7) +
  ggplot2::geom_point(ggplot2::aes(y = label, x = prop_climate_enriched * 100,
                                   colour = direction, shape = "Candidates"),
                      size = 2.9) +
  ggplot2::scale_shape_manual(values = c(Background = 21, Candidates = 19),
                              name = NULL) +
  ggplot2::scale_colour_manual(
    values = c(`enriched in candidates` = "#C44E52",
               `depleted in candidates` = "#4C72B0",
               `not significant` = "grey72"), name = NULL) +
  ggplot2::labs(
    x = "TEs carrying an enriched binding site (%)", y = NULL,
    title = "Transcription-factor binding sites: candidates vs a matched background",
    subtitle = sprintf(
      "%d climate-associated TEs vs %d non-climate-associated TEs, matched by superfamily; Fisher's exact test, BH-adjusted",
      n_fg_total, n_bg_total),
    caption = "Odds ratios shown for families differing at BH-adjusted p < 0.05.") +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "top",
                 axis.text.y = ggplot2::element_text(size = 8),
                 plot.title = ggplot2::element_text(face = "bold", size = 12),
                 plot.subtitle = ggplot2::element_text(size = 8.5, colour = "grey30"),
                 plot.caption = ggplot2::element_text(size = 7.5, colour = "grey40",
                                                      hjust = 0),
                 plot.background = ggplot2::element_rect(fill = "white",
                                                         colour = NA))
save_plot(p_dumb, "Fig3C_enrichment_dumbbell", width = 8, height = 8)

# --- 8b-ii. Per-TE detail, restricted to the families that differ -----------
# Only the significantly different families are kept, so the grid is small
# enough to label properly. TEs are ordered by superfamily, then by how many
# of those families they carry.
focus_fams <- comparison %>%
  dplyr::filter(fisher_padj < 0.05) %>%
  dplyr::arrange(dplyr::desc(odds_ratio)) %>%
  dplyr::pull(Family)

if (length(focus_fams) > 0) {
  fg_sf <- selection %>%
    dplyr::filter(group == "climate_associated") %>%
    dplyr::select(TE_ID, fimo_group)

  heat <- tidyr::expand_grid(
    TE_ID = fg_sf$TE_ID, Family = focus_fams) %>%
    dplyr::left_join(
      fg$pairs %>% dplyr::select(TE_ID = sequence_name, Family,
                                 median_neglog10_q, enriched),
      by = c("TE_ID", "Family")) %>%
    dplyr::left_join(fg_sf, by = "TE_ID")

  te_order <- heat %>%
    dplyr::group_by(TE_ID, fimo_group) %>%
    dplyr::summarise(n_en = sum(enriched, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(fimo_group, dplyr::desc(n_en)) %>%
    dplyr::pull(TE_ID)

  heat <- heat %>%
    dplyr::mutate(
      TE_ID = factor(TE_ID, levels = te_order),
      Family = factor(Family, levels = rev(focus_fams)),
      # short label: TE_124385_LTR -> 124385
      te_short = sub("^TE_0*", "", sub("_.*$", "", sub("^TE_", "", TE_ID))))

  p_heat <- ggplot2::ggplot(heat,
                            ggplot2::aes(x = TE_ID, y = Family)) +
    ggplot2::geom_tile(ggplot2::aes(fill = median_neglog10_q),
                       colour = "white", linewidth = 0.4) +
    # a ring marks the pairs that clear the enrichment threshold
    ggplot2::geom_point(data = heat %>% dplyr::filter(enriched %in% TRUE),
                        shape = 21, size = 1.5, colour = "grey10",
                        fill = NA, stroke = 0.45) +
    ggplot2::scale_fill_gradientn(
      colours = c("#FFF7EC", "#FDD49E", "#FC8D59", "#B30000"),
      na.value = "grey93",
      name = expression(Median~-log[10]*(q))) +
    ggplot2::facet_grid(. ~ fimo_group, scales = "free_x", space = "free_x") +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = "Binding-site strength per candidate TE",
      subtitle = paste0("Only the ", length(focus_fams),
                        " TF families that differ from the background are shown; ",
                        "grey = no significant motif hit; ring = median -log10(q) > 1.3")) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1,
                                          size = 6.5),
      axis.text.y = ggplot2::element_text(size = 9),
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 8.5, colour = "grey30"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  save_plot(p_heat, "Fig3C_heatmap_clean", width = 11, height = 4.6)

  p_c <- p_dumb / p_heat + patchwork::plot_layout(heights = c(1.5, 1))
  save_plot(p_c, "Fig3C_composite", width = 11.5, height = 12)
}

# =============================================================================
# 9. Session log
# =============================================================================
log_session("05_tfbs_enrichment")
message("=== 05_tfbs_enrichment.R finished (PART 1 + PART 2) ===")
