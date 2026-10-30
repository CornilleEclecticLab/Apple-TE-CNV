# =============================================================================
# Script:  13_tfbs_single_run.R
# Purpose: Redo the TFBS comparison with ONE FIMO run over all sequences,
#          instead of one run per TE superfamily.
#
#          WHY THIS MATTERS. FIMO's q-value is not a property of a
#          motif-sequence pair on its own: it is an FDR estimated WITHIN the
#          run, using the p-value distribution of that run (FIMO reports the
#          Storey pi_0 it estimates each time). The published pipeline, and
#          05_tfbs_enrichment.R which reproduced it, scanned each superfamily
#          separately - so the 38 candidates were split across five small runs
#          (4, 11, 13, 6 and 4 sequences) while the background was split across
#          five much larger ones (32-104 sequences), each with its own
#          calibration. Observed pi_0 ranged from 0.84 to 1.0 across runs.
#
#          Comparing q-values across separately calibrated runs is not sound,
#          and the direction of the resulting bias is not obvious a priori: a
#          set whose p-value distribution is richer in strong hits gets more
#          generous q-values, which is exactly the quantity being compared.
#
#          The fix is simply to scan every sequence in a single FIMO run, so
#          all hits share one calibration, and only then split by group. This
#          script does that and reports how much the conclusions move, which is
#          the number that decides whether the published panel needs redoing.
#
# NAMING NOTE: "Model A" is the climate GLM, "Model B" the spatial model that
# yields VS, and "Model A-PC" the PC-covariate climate model.
#
# Author:  Antho (analysis run with Claude Code)
# Date:    2026-09-01
# Project: TE_Project / final_pipeline (Genome Biology revision)
#
# Inputs:
#   outputs/05_tfbs_enrichment/tables/background_TE_selection.csv
#   outputs/05_tfbs_enrichment/tables/enrichment_comparison_by_TF_family.csv
#   outputs/12_tfbs_robustness/tables/full_pool_CMH.csv
#   data/TE_consensus_db.fasta, data/Mdo_TF_binding_motifs.meme, data/Mdo_TF_list.txt
#
# Outputs (all under outputs/13_tfbs_single_run/):
#   fasta/all_sequences.fasta            everything scanned together
#   fimo_out/fimo.tsv                    the single run
#   tables/enrichment_single_run.csv     comparison with comparable q-values
#   tables/per_superfamily_breakdown.csv direction within each superfamily
#   tables/split_vs_single_run.csv       how much the artefact mattered
#   figures/*.pdf + *.svg
#   ../../logs/13_tfbs_single_run_*_sessionInfo.txt
# =============================================================================

# --- Setup -------------------------------------------------------------------
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(Biostrings)
  library(svglite)
})

set.seed(20260824)

# --- Parameters --------------------------------------------------------------
Q_THRESHOLD   <- 0.05
ENRICH_THRESH <- 1.3
MEME_BIN      <- path.expand("~/opt/meme/bin")
SF            <- c("Copia", "Gypsy", "LTR", "NoCat", "TIR")
# A family counts as consistently enriched only if it moves the same way in at
# least this many superfamilies; the stratified test alone can average over a
# reversal (see 12_tfbs_robustness.R).
MIN_CONSISTENT_SF <- 4

dir_out  <- file.path(dir_outputs, "13_tfbs_single_run")
dir_tab  <- file.path(dir_out, "tables")
dir_fig  <- file.path(dir_out, "figures")
dir_fa   <- file.path(dir_out, "fasta")
dir_fimo <- file.path(dir_out, "fimo_out")
for (d in c(dir_out, dir_tab, dir_fig, dir_fa, dir_fimo)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

save_plot <- function(p, name, width = 7, height = 5) {
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".pdf")), p, width = width,
                  height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".svg")), p, width = width,
                  height = height, device = svglite::svglite)
}

message("=== 13_tfbs_single_run.R ===")

# =============================================================================
# 1. Assemble every sequence into one FASTA
# =============================================================================
sel <- readr::read_csv(file.path(dir_outputs, "05_tfbs_enrichment", "tables",
                                 "background_TE_selection.csv"),
                       show_col_types = FALSE)
fg_ids <- sel$TE_ID[sel$group == "climate_associated"]

modelA_raw <- readr::read_csv(file.path(dir_data, "modelA_fig3_raw_p.csv"),
                              show_col_types = FALSE)
modelAPC <- readr::read_csv(file.path(dir_outputs, "02_candidate_TEs", "tables",
                                      "glm_modelA_PC_by_group.csv"),
                            show_col_types = FALSE)
climate_any <- union(
  unique(modelA_raw$TE_ID[!is.na(modelA_raw$p_value) & modelA_raw$p_value < 0.01]),
  unique(modelAPC$TE_ID[!is.na(modelAPC$p_value) & modelAPC$p_value < 0.01]))

sfam <- function(x) {
  sx <- sub("^TE_[0-9]+_", "", x)
  dplyr::case_when(sx == "LTR_Copia" ~ "Copia", sx == "LTR_Gypsy" ~ "Gypsy",
                   sx == "LTR" ~ "LTR", sx == "noCat" ~ "NoCat",
                   sx == "TIR" ~ "TIR", TRUE ~ NA_character_)
}
bg_ids <- setdiff(sort(unique(modelA_raw$TE_ID)), union(climate_any, fg_ids))
bg_ids <- bg_ids[!is.na(sfam(bg_ids))]

groups <- tibble::tibble(
  TE_ID = c(fg_ids, bg_ids),
  group = c(rep("climate_associated", length(fg_ids)),
            rep("background", length(bg_ids))),
  superfamily = sfam(c(fg_ids, bg_ids)))
readr::write_csv(groups, file.path(dir_tab, "sequence_groups.csv"))
message(sprintf("Scanning together: %d candidates + %d background = %d sequences",
                length(fg_ids), length(bg_ids), nrow(groups)))

fa_all <- file.path(dir_fa, "all_sequences.fasta")
if (!file.exists(fa_all)) {
  seqs <- Biostrings::readDNAStringSet(file.path(dir_data, "TE_consensus_db.fasta"))
  names(seqs) <- sub(" .*$", "", names(seqs))
  Biostrings::writeXStringSet(seqs[groups$TE_ID], fa_all)
}

# =============================================================================
# 2. One FIMO run
# =============================================================================
fimo_bin <- file.path(MEME_BIN, "fimo")
motifs <- file.path(dir_data, "Mdo_TF_binding_motifs.meme")
tsv <- file.path(dir_fimo, "fimo.tsv")
if (!file.exists(tsv)) {
  stopifnot(file.exists(fimo_bin))
  message("Running FIMO once over all sequences (this takes a while) ...")
  od <- file.path(dir_fimo, "raw")
  system2(fimo_bin, c("--oc", shQuote(od), shQuote(motifs), shQuote(fa_all)),
          stdout = NULL, stderr = NULL)
  stopifnot(file.exists(file.path(od, "fimo.tsv")))
  file.copy(file.path(od, "fimo.tsv"), tsv, overwrite = TRUE)
}

tf_family <- readr::read_tsv(file.path(dir_data, "Mdo_TF_list.txt"),
                             show_col_types = FALSE)
pairs <- readr::read_tsv(tsv, show_col_types = FALSE, comment = "#",
                         progress = FALSE) %>%
  dplyr::filter(!is.na(`q-value`), `q-value` < Q_THRESHOLD) %>%
  dplyr::left_join(tf_family, by = c("motif_id" = "TF_ID")) %>%
  dplyr::filter(!is.na(Family)) %>%
  dplyr::group_by(TE_ID = sequence_name, Family) %>%
  dplyr::summarise(med = stats::median(-log10(`q-value`)), .groups = "drop") %>%
  dplyr::mutate(enriched = med > ENRICH_THRESH)

# =============================================================================
# 3. Enrichment, now on one calibration
# =============================================================================
fams <- sort(unique(pairs$Family))
n_fg <- length(fg_ids); n_bg <- length(bg_ids)

single <- purrr::map_dfr(fams, function(f) {
  a <- sum(pairs$Family == f & pairs$enriched & pairs$TE_ID %in% fg_ids)
  b <- sum(pairs$Family == f & pairs$enriched & pairs$TE_ID %in% bg_ids)
  ft <- stats::fisher.test(matrix(c(a, n_fg - a, b, n_bg - b), nrow = 2))
  tibble::tibble(Family = f, n_cand = a, n_bg = b,
                 pct_cand = 100 * a / n_fg, pct_bg = 100 * b / n_bg,
                 odds_ratio = unname(ft$estimate), p = ft$p.value)
}) %>%
  dplyr::mutate(padj = stats::p.adjust(p, method = "BH")) %>%
  dplyr::arrange(padj)

# Direction within each superfamily, so a reversal cannot hide inside an average.
breakdown <- purrr::map_dfr(fams, function(f) {
  purrr::map_dfr(SF, function(s) {
    fi <- groups$TE_ID[groups$group == "climate_associated" &
                         groups$superfamily == s]
    bi <- groups$TE_ID[groups$group == "background" & groups$superfamily == s]
    a <- sum(pairs$Family == f & pairs$enriched & pairs$TE_ID %in% fi)
    b <- sum(pairs$Family == f & pairs$enriched & pairs$TE_ID %in% bi)
    tibble::tibble(Family = f, superfamily = s,
                   n_cand = a, n_cand_total = length(fi),
                   n_bg = b, n_bg_total = length(bi),
                   pct_cand = 100 * a / max(1, length(fi)),
                   pct_bg = 100 * b / max(1, length(bi)))
  })
}) %>%
  dplyr::mutate(direction = dplyr::case_when(
    pct_cand > pct_bg ~ "up", pct_cand < pct_bg ~ "down", TRUE ~ "flat"))
readr::write_csv(breakdown, file.path(dir_tab, "per_superfamily_breakdown.csv"))

consistency <- breakdown %>%
  dplyr::group_by(Family) %>%
  dplyr::summarise(n_up = sum(direction == "up"),
                   n_down = sum(direction == "down"),
                   n_flat = sum(direction == "flat"), .groups = "drop")

single <- single %>%
  dplyr::left_join(consistency, by = "Family") %>%
  dplyr::mutate(
    consistent = dplyr::case_when(
      odds_ratio > 1 & n_up >= MIN_CONSISTENT_SF ~ TRUE,
      odds_ratio < 1 & n_down >= MIN_CONSISTENT_SF ~ TRUE,
      TRUE ~ FALSE),
    call = dplyr::case_when(
      padj >= 0.05 ~ "not significant",
      !consistent ~ "significant but direction reverses across superfamilies",
      odds_ratio > 1 ~ "enriched",
      TRUE ~ "depleted"))
readr::write_csv(single, file.path(dir_tab, "enrichment_single_run.csv"))

message("Single-run result:")
print(as.data.frame(single %>% dplyr::filter(padj < 0.05) %>%
  dplyr::transmute(Family, cand = sprintf("%d/%d", n_cand, n_fg),
                   bg = sprintf("%d/%d", n_bg, n_bg),
                   OR = round(odds_ratio, 2), padj = signif(padj, 2),
                   up = n_up, down = n_down, call)))

# =============================================================================
# 4. How much did the split-run artefact matter?
# =============================================================================
split <- readr::read_csv(file.path(dir_outputs, "12_tfbs_robustness", "tables",
                                   "full_pool_CMH.csv"), show_col_types = FALSE)
cmp <- single %>%
  dplyr::select(Family, single_OR = odds_ratio, single_padj = padj, call) %>%
  dplyr::full_join(split %>% dplyr::select(Family, split_OR = common_OR,
                                           split_padj = cmh_padj),
                   by = "Family") %>%
  dplyr::mutate(
    sig_single = !is.na(single_padj) & single_padj < 0.05,
    sig_split  = !is.na(split_padj) & split_padj < 0.05,
    agreement = dplyr::case_when(
      sig_single & sig_split ~ "significant in both",
      sig_single & !sig_split ~ "only with one FIMO run",
      !sig_single & sig_split ~ "only with per-superfamily runs",
      TRUE ~ "neither")) %>%
  dplyr::arrange(dplyr::desc(sig_single | sig_split), single_padj)
readr::write_csv(cmp, file.path(dir_tab, "split_vs_single_run.csv"))
message("Split-run versus single-run:")
print(as.data.frame(cmp %>% dplyr::filter(sig_single | sig_split) %>%
  dplyr::transmute(Family, single_OR = round(single_OR, 2),
                   split_OR = round(split_OR, 2), agreement)))

# =============================================================================
# 5. Figure
# =============================================================================
p <- single %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::mutate(Family = stats::reorder(Family, pct_cand - pct_bg)) %>%
  dplyr::select(Family, call, Candidates = pct_cand, Background = pct_bg) %>%
  tidyr::pivot_longer(c(Candidates, Background), names_to = "set",
                      values_to = "pct") %>%
  ggplot2::ggplot(ggplot2::aes(x = pct, y = Family, fill = set)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75),
                    width = 0.68) +
  ggplot2::facet_grid(call ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_fill_manual(values = c(Candidates = "#C44E52",
                                        Background = "#4C72B0"), name = NULL) +
  ggplot2::labs(
    x = "TEs carrying an enriched binding site (%)", y = NULL,
    title = "TFBS enrichment with a single FIMO calibration",
    subtitle = sprintf("%d candidates vs %d background TEs scanned in one run",
                       n_fg, n_bg)) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "top",
                 strip.text.y = ggplot2::element_text(angle = 0, size = 7.5),
                 plot.background = ggplot2::element_rect(fill = "white",
                                                         colour = NA))
save_plot(p, "enrichment_single_run", width = 9, height = 6.5)

# =============================================================================
# 6. Session log
# =============================================================================
log_session("13_tfbs_single_run")
message("=== 13_tfbs_single_run.R finished ===")
