# =============================================================================
# Script:  12_tfbs_robustness.R
# Purpose: Stress-test the TFBS enrichment result of 05_tfbs_enrichment.R
#          against three objections that the headline comparison alone cannot
#          answer.
#
#          OBJECTION 1 - "the groups are not the same size".
#          The comparison is 38 candidates against 304 background TEs. Fisher's
#          exact test conditions on the margins, so unequal group sizes change
#          the PRECISION of the estimate, not its expectation: the odds ratio
#          is not biased by the 8:1 ratio. What unequal sizes do change is how
#          coarse the candidate side is - one TE is 2.6% there and 0.33% in the
#          background - so a family can reach significance on very few
#          sequences. This script therefore re-runs the whole comparison on
#          BALANCED draws (38 vs 38, superfamily-matched, many times) and
#          reports how often each family survives.
#
#          OBJECTION 2 - "why not random sequences?".
#          Random or shuffled sequences answer a different question: whether
#          these TEs carry more motifs than their own base composition
#          predicts. That is a weaker test, and almost guaranteed to be
#          positive, because TE consensus sequences are AT-rich and repetitive.
#          The reviewer asked specifically for the comparison against
#          NON-CLIMATE-ASSOCIATED TEs, which is what the main analysis does and
#          is the stronger control: it holds TE-ness constant and varies only
#          the climate association. A compositional control is still worth
#          having as a complement, so this script adds one - each candidate is
#          dinucleotide-shuffled (preserving base and dinucleotide composition)
#          and rescanned with the same FIMO settings.
#
#          OBJECTION 3 - "where do the background sequences come from?".
#          They are real TEs from the same consensus library: the 550 residual
#          TEs reaching p < 0.01 for NO bioclimatic variable in ANY of the three
#          M. sylvestris groups, under either Model A or Model A-PC, sampled to
#          match the candidate superfamily composition exactly. The provenance
#          table is written out here so it can go straight into a supplement.
#
# NAMING NOTE: "Model A" is the climate GLM, "Model B" the spatial model that
# yields VS, and "Model A-PC" the PC-covariate climate model.
#
# Author:  Antho (analysis run with Claude Code)
# Date:    2026-08-25
# Project: TE_Project / final_pipeline (Genome Biology revision)
#
# Inputs:
#   outputs/05_tfbs_enrichment/tables/background_TE_selection.csv
#   outputs/05_tfbs_enrichment/tables/enrichment_comparison_by_TF_family.csv
#   outputs/05_tfbs_enrichment/fimo_out/{foreground,background}/fimo_*.tsv
#   outputs/05_tfbs_enrichment/fasta/foreground/*.fasta
#   data/Mdo_TF_list.txt, data/Mdo_TF_binding_motifs.meme
#
# Outputs (all under outputs/12_tfbs_robustness/):
#   tables/enrichment_raw_counts.csv        counts, not just proportions
#   tables/balanced_subsampling.csv         38 vs 38, how often each family holds
#   tables/shuffled_control.csv             candidates vs their own shuffles
#   tables/background_provenance.csv        where the 304 come from
#   tables/full_pool_CMH.csv                every eligible background TE, CMH
#   figures/*.pdf + *.svg
#   ../../logs/12_tfbs_robustness_*_sessionInfo.txt
# =============================================================================

# --- Setup -------------------------------------------------------------------
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(Biostrings)
  library(svglite)
})

set.seed(20260824)

# --- Parameters --------------------------------------------------------------
N_SUBSAMPLE   <- 1000    # balanced draws
Q_THRESHOLD   <- 0.05
ENRICH_THRESH <- 1.3
SHUFFLE_COPIES<- 5       # shuffled replicates per candidate sequence
MEME_BIN      <- path.expand("~/opt/meme/bin")
MEME_LIBEXEC  <- path.expand("~/opt/meme/libexec/meme-5.5.7")

dir_src <- file.path(dir_outputs, "05_tfbs_enrichment")
dir_out <- file.path(dir_outputs, "12_tfbs_robustness")
dir_tab <- file.path(dir_out, "tables")
dir_fig <- file.path(dir_out, "figures")
dir_sh  <- file.path(dir_out, "fasta_shuffled")
dir_shf <- file.path(dir_out, "fimo_shuffled")
dir_extra <- file.path(dir_out, "fasta_extra_background")
dir_extraf <- file.path(dir_out, "fimo_extra_background")
for (d in c(dir_out, dir_tab, dir_fig, dir_sh, dir_shf, dir_extra,
            dir_extraf)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

save_plot <- function(p, name, width = 7, height = 5) {
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".pdf")), p, width = width,
                  height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".svg")), p, width = width,
                  height = height, device = svglite::svglite)
}

shuffler <- file.path(MEME_LIBEXEC, "fasta-shuffle-letters")
fimo_bin <- file.path(MEME_BIN, "fimo")
motifs   <- file.path(dir_data, "Mdo_TF_binding_motifs.meme")

message("=== 12_tfbs_robustness.R ===")

# =============================================================================
# 1. Reload the enrichment pipeline exactly as the main script defines it
# =============================================================================
tf_family <- readr::read_tsv(file.path(dir_data, "Mdo_TF_list.txt"),
                             show_col_types = FALSE)
SF <- c("Copia", "Gypsy", "LTR", "NoCat", "TIR")

read_pairs <- function(dir_fimo) {
  f <- file.path(dir_fimo, sprintf("fimo_%s.tsv", SF))
  f <- f[file.exists(f)]
  if (!length(f)) return(NULL)
  purrr::map_dfr(f, function(p)
    readr::read_tsv(p, show_col_types = FALSE, comment = "#", progress = FALSE)) %>%
    dplyr::filter(!is.na(`q-value`), `q-value` < Q_THRESHOLD) %>%
    dplyr::mutate(neglog10_q = -log10(`q-value`)) %>%
    dplyr::left_join(tf_family, by = c("motif_id" = "TF_ID")) %>%
    dplyr::filter(!is.na(Family)) %>%
    dplyr::group_by(TE_ID = sequence_name, Family) %>%
    dplyr::summarise(med = stats::median(neglog10_q), .groups = "drop") %>%
    dplyr::mutate(enriched = med > ENRICH_THRESH)
}

fg_pairs <- read_pairs(file.path(dir_src, "fimo_out", "foreground"))
bg_pairs <- read_pairs(file.path(dir_src, "fimo_out", "background"))
stopifnot(!is.null(fg_pairs), !is.null(bg_pairs))

selection <- readr::read_csv(file.path(dir_src, "tables",
                                       "background_TE_selection.csv"),
                             show_col_types = FALSE)
fg_ids <- selection$TE_ID[selection$group == "climate_associated"]
bg_ids <- selection$TE_ID[selection$group == "background"]

# Every residual TE that is climate-associated nowhere: the full pool the 304
# were drawn from.
modelA_raw <- readr::read_csv(file.path(dir_data, "modelA_fig3_raw_p.csv"),
                              show_col_types = FALSE)
modelAPC <- readr::read_csv(file.path(dir_outputs, "02_candidate_TEs", "tables",
                                      "glm_modelA_PC_by_group.csv"),
                            show_col_types = FALSE)
climate_any <- union(
  unique(modelA_raw$TE_ID[!is.na(modelA_raw$p_value) & modelA_raw$p_value < 0.01]),
  unique(modelAPC$TE_ID[!is.na(modelAPC$p_value) & modelAPC$p_value < 0.01]))
eligible_pool <- setdiff(sort(unique(modelA_raw$TE_ID)),
                         union(climate_any, fg_ids))
message(sprintf("Candidates: %d | background: %d", length(fg_ids), length(bg_ids)))

# =============================================================================
# 2. Objection 3 - provenance of the background, written out explicitly
# =============================================================================
prov <- selection %>%
  dplyr::count(group, fimo_group, name = "n") %>%
  tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0) %>%
  dplyr::rename(superfamily = fimo_group) %>%
  dplyr::mutate(ratio = round(background / climate_associated, 1))
readr::write_csv(prov, file.path(dir_tab, "background_provenance.csv"))
message("Background composition (matched by superfamily):")
print(as.data.frame(prov))

# =============================================================================
# 3. Objection 1a - report counts, not only proportions
# =============================================================================
fams <- sort(union(unique(fg_pairs$Family), unique(bg_pairs$Family)))

count_table <- function(pairs, ids) {
  vapply(fams, function(f) sum(pairs$Family == f & pairs$enriched &
                                 pairs$TE_ID %in% ids), integer(1))
}
k_fg <- count_table(fg_pairs, fg_ids)
k_bg <- count_table(bg_pairs, bg_ids)

raw <- tibble::tibble(
  Family = fams,
  n_candidates_with_site = k_fg, n_candidates = length(fg_ids),
  n_background_with_site = k_bg, n_background = length(bg_ids),
  pct_candidates = round(100 * k_fg / length(fg_ids), 1),
  pct_background = round(100 * k_bg / length(bg_ids), 1)) %>%
  dplyr::arrange(dplyr::desc(pct_candidates - pct_background))
readr::write_csv(raw, file.path(dir_tab, "enrichment_raw_counts.csv"))

main <- readr::read_csv(file.path(dir_src, "tables",
                                  "enrichment_comparison_by_TF_family.csv"),
                        show_col_types = FALSE)
sig_fams <- main$Family[main$fisher_padj < 0.05]
fragile <- raw$Family[raw$Family %in% sig_fams & raw$n_candidates_with_site < 5]
message(sprintf("Significant families resting on fewer than 5 candidate TEs: %s",
                if (length(fragile)) paste(fragile, collapse = ", ") else "none"))

# =============================================================================
# 4. Objection 1b - balanced subsampling, 38 vs 38
# =============================================================================
# Each draw keeps the superfamily composition of the candidates, so the only
# thing that changes relative to the main analysis is that the background side
# now has exactly the same number of sequences as the candidate side.
sel_bg <- selection %>% dplyr::filter(group == "background")
need <- selection %>% dplyr::filter(group == "climate_associated") %>%
  dplyr::count(fimo_group, name = "k")

message(sprintf("Balanced subsampling: %d draws of %d background TEs ...",
                N_SUBSAMPLE, length(fg_ids)))
hits <- matrix(FALSE, nrow = N_SUBSAMPLE, ncol = length(fams),
               dimnames = list(NULL, fams))
dir_hits <- hits

for (i in seq_len(N_SUBSAMPLE)) {
  draw <- unlist(lapply(seq_len(nrow(need)), function(j) {
    pool <- sel_bg$TE_ID[sel_bg$fimo_group == need$fimo_group[j]]
    sample(pool, need$k[j], replace = FALSE)
  }))
  kb <- count_table(bg_pairs, draw)
  for (m in seq_along(fams)) {
    ct <- matrix(c(k_fg[m], length(fg_ids) - k_fg[m],
                   kb[m], length(draw) - kb[m]), nrow = 2)
    p <- stats::fisher.test(ct)$p.value
    hits[i, m] <- p < 0.05
    dir_hits[i, m] <- k_fg[m] / length(fg_ids) > kb[m] / length(draw)
  }
}

balanced <- tibble::tibble(
  Family = fams,
  n_candidates_with_site = k_fg,
  prop_draws_significant = round(colMeans(hits), 3),
  prop_draws_same_direction = round(colMeans(dir_hits), 3),
  significant_in_main = fams %in% sig_fams) %>%
  dplyr::arrange(dplyr::desc(prop_draws_significant))
readr::write_csv(balanced, file.path(dir_tab, "balanced_subsampling.csv"))
message("Balanced 38 vs 38 (unadjusted p < 0.05 per draw):")
print(as.data.frame(balanced %>% dplyr::filter(significant_in_main)))

# =============================================================================
# 5. Objection 2 - dinucleotide-shuffled compositional control
# =============================================================================
# fasta-shuffle-letters -kmer 2 preserves each sequence's dinucleotide
# composition while destroying its motif content, which isolates "more sites
# than base composition predicts" from "more sites than other TEs".

if (file.exists(shuffler) && file.exists(fimo_bin)) {
  for (sf in SF) {
    fa <- file.path(dir_src, "fasta", "foreground",
                    sprintf("foreground_%s.fasta", sf))
    if (!file.exists(fa)) next
    out <- file.path(dir_sh, sprintf("shuffled_%s.fasta", sf))
    if (!file.exists(out)) {
      system2(shuffler, c("-kmer", "2", "-copies", SHUFFLE_COPIES,
                          "-seed", "20260824", shQuote(fa)), stdout = out)
    }
    od <- file.path(dir_shf, sprintf("raw_%s", sf))
    tsv <- file.path(dir_shf, sprintf("fimo_%s.tsv", sf))
    if (!file.exists(tsv)) {
      message("  FIMO on shuffled ", sf, " ...")
      system2(fimo_bin, c("--oc", shQuote(od), shQuote(motifs), shQuote(out)),
              stdout = NULL, stderr = NULL)
      if (file.exists(file.path(od, "fimo.tsv")))
        file.copy(file.path(od, "fimo.tsv"), tsv, overwrite = TRUE)
    }
  }

  sh_pairs <- read_pairs(dir_shf)
  if (!is.null(sh_pairs)) {
    sh_ids <- unique(sh_pairs$TE_ID)
    n_sh <- length(unique(sub("_shuf_[0-9]+$", "", sh_ids)))
    n_sh_total <- length(sh_ids)
    k_sh <- count_table(sh_pairs, sh_ids)

    shuffled <- tibble::tibble(
      Family = fams,
      pct_candidates = round(100 * k_fg / length(fg_ids), 1),
      pct_shuffled = round(100 * k_sh / n_sh_total, 1),
      pct_background = round(100 * k_bg / length(bg_ids), 1)) %>%
      dplyr::mutate(
        fisher_p_vs_shuffled = purrr::map2_dbl(k_fg, k_sh, function(a, b)
          stats::fisher.test(matrix(c(a, length(fg_ids) - a,
                                      b, n_sh_total - b), nrow = 2))$p.value)) %>%
      dplyr::mutate(padj_vs_shuffled = stats::p.adjust(fisher_p_vs_shuffled,
                                                       method = "BH")) %>%
      dplyr::arrange(padj_vs_shuffled)
    readr::write_csv(shuffled, file.path(dir_tab, "shuffled_control.csv"))
    message(sprintf("Dinucleotide-shuffled control: %d shuffled sequences", n_sh_total))
    message("Families exceeding their own shuffled sequences (BH < 0.05): ",
            paste(shuffled$Family[shuffled$padj_vs_shuffled < 0.05], collapse = ", "))
  }
} else {
  message("MEME utilities not found; the shuffled control was skipped.")
  shuffled <- NULL
}

# =============================================================================
# 5b. Objection 1c - use EVERY eligible background TE, not a matched subsample
# =============================================================================
# The 8x matched draw of 304 was capped by Copia (33 eligible TEs for the 4
# candidates), which left 246 eligible TEs unused. Nothing forces a subsample:
# the whole eligible pool can be used if the superfamily imbalance is handled
# in the TEST rather than in the sampling. A Cochran-Mantel-Haenszel test does
# exactly that - it estimates a common odds ratio across superfamily strata,
# so composition is controlled without discarding data.
extra_ids <- setdiff(eligible_pool, bg_ids)
message(sprintf("Full eligible background pool: %d TEs (%d already scanned, %d to add)",
                length(eligible_pool), length(bg_ids), length(extra_ids)))

if (length(extra_ids) > 0 && file.exists(fimo_bin)) {
  seqs_all <- Biostrings::readDNAStringSet(file.path(dir_data, "TE_consensus_db.fasta"))
  names(seqs_all) <- sub(" .*$", "", names(seqs_all))
  sfam <- function(x) {
    sx <- sub("^TE_[0-9]+_", "", x)
    dplyr::case_when(sx == "LTR_Copia" ~ "Copia", sx == "LTR_Gypsy" ~ "Gypsy",
                     sx == "LTR" ~ "LTR", sx == "noCat" ~ "NoCat",
                     sx == "TIR" ~ "TIR", TRUE ~ NA_character_)
  }
  for (sf in SF) {
    ids <- extra_ids[sfam(extra_ids) %in% sf]
    if (!length(ids)) next
    fa <- file.path(dir_extra, sprintf("extra_%s.fasta", sf))
    if (!file.exists(fa)) Biostrings::writeXStringSet(seqs_all[ids], fa)
    tsv <- file.path(dir_extraf, sprintf("fimo_%s.tsv", sf))
    if (!file.exists(tsv)) {
      message("  FIMO on the remaining ", sf, " background (", length(ids), ") ...")
      od <- file.path(dir_extraf, sprintf("raw_%s", sf))
      system2(fimo_bin, c("--oc", shQuote(od), shQuote(motifs), shQuote(fa)),
              stdout = NULL, stderr = NULL)
      if (file.exists(file.path(od, "fimo.tsv")))
        file.copy(file.path(od, "fimo.tsv"), tsv, overwrite = TRUE)
    }
  }
  extra_pairs <- read_pairs(dir_extraf)
  full_pairs <- dplyr::bind_rows(bg_pairs, extra_pairs)
  full_ids <- union(bg_ids, extra_ids)
  full_ids <- full_ids[!is.na(sfam(full_ids))]

  fg_sf   <- sfam(fg_ids)
  full_sf <- sfam(full_ids)

  cmh <- purrr::map_dfr(fams, function(f) {
    # one 2x2 table per superfamily stratum
    arr <- array(0L, dim = c(2, 2, length(SF)))
    for (k in seq_along(SF)) {
      fi <- fg_ids[fg_sf == SF[k]]
      bi <- full_ids[full_sf == SF[k]]
      a <- sum(fg_pairs$Family == f & fg_pairs$enriched & fg_pairs$TE_ID %in% fi)
      b <- sum(full_pairs$Family == f & full_pairs$enriched &
                 full_pairs$TE_ID %in% bi)
      arr[, , k] <- matrix(c(a, length(fi) - a, b, length(bi) - b), nrow = 2)
    }
    keep <- apply(arr, 3, function(m) sum(m) > 0 && all(colSums(m) > 0))
    if (sum(keep) < 1) return(NULL)
    tst <- try(stats::mantelhaen.test(arr[, , keep, drop = FALSE],
                                      exact = FALSE, correct = TRUE),
               silent = TRUE)
    if (inherits(tst, "try-error")) return(NULL)
    a_tot <- sum(arr[1, 1, keep]); b_tot <- sum(arr[1, 2, keep])
    tibble::tibble(Family = f,
                   n_candidates_with_site = a_tot,
                   n_background_with_site = b_tot,
                   n_background_total = length(full_ids),
                   common_OR = unname(tst$estimate),
                   cmh_p = tst$p.value)
  }) %>%
    dplyr::mutate(cmh_padj = stats::p.adjust(cmh_p, method = "BH")) %>%
    dplyr::arrange(cmh_padj)
  readr::write_csv(cmh, file.path(dir_tab, "full_pool_CMH.csv"))
  message(sprintf("CMH against the full pool of %d background TEs, enriched at BH < 0.05: %s",
                  length(full_ids),
                  paste(cmh$Family[cmh$cmh_padj < 0.05 & cmh$common_OR > 1],
                        collapse = ", ")))
  message(sprintf("  depleted at BH < 0.05: %s",
                  paste(cmh$Family[cmh$cmh_padj < 0.05 & cmh$common_OR < 1],
                        collapse = ", ")))
} else {
  cmh <- NULL
}

# =============================================================================
# 6. Figures
# =============================================================================
p_bal <- balanced %>%
  dplyr::filter(significant_in_main) %>%
  dplyr::mutate(Family = stats::reorder(Family, prop_draws_significant),
                fragile = n_candidates_with_site < 5) %>%
  ggplot2::ggplot(ggplot2::aes(x = prop_draws_significant, y = Family,
                               fill = fragile)) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::geom_vline(xintercept = 0.5, linetype = 2, colour = "grey40") +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%d/38", n_candidates_with_site)),
                     hjust = -0.15, size = 3) +
  ggplot2::scale_fill_manual(values = c(`FALSE` = "#4C72B0", `TRUE` = "#C44E52"),
                             labels = c("5 or more candidate TEs",
                                        "fewer than 5 candidate TEs"),
                             name = NULL) +
  ggplot2::scale_x_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25)) +
  ggplot2::labs(
    x = "Fraction of balanced draws still significant (p < 0.05)", y = NULL,
    title = "Does the result survive equal group sizes?",
    subtitle = sprintf("%d draws of 38 background TEs, superfamily-matched; labels give the candidate count",
                       N_SUBSAMPLE)) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "top",
                 plot.background = ggplot2::element_rect(fill = "white",
                                                         colour = NA))
save_plot(p_bal, "balanced_subsampling", width = 8, height = 5.5)

if (!is.null(shuffled)) {
  p_sh <- shuffled %>%
    dplyr::filter(Family %in% sig_fams) %>%
    dplyr::select(Family, Candidates = pct_candidates,
                  `Background TEs` = pct_background,
                  `Shuffled candidates` = pct_shuffled) %>%
    tidyr::pivot_longer(-Family, names_to = "set", values_to = "pct") %>%
    dplyr::mutate(set = factor(set, levels = c("Candidates", "Background TEs",
                                               "Shuffled candidates"))) %>%
    ggplot2::ggplot(ggplot2::aes(x = pct, y = stats::reorder(Family, pct),
                                 fill = set)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                      width = 0.72) +
    ggplot2::scale_fill_manual(values = c(Candidates = "#C44E52",
                                          `Background TEs` = "#4C72B0",
                                          `Shuffled candidates` = "grey65"),
                               name = NULL) +
    ggplot2::labs(x = "Sequences carrying an enriched binding site (%)", y = NULL,
                  title = "Two controls, two different questions",
                  subtitle = "Background TEs hold TE-ness constant; shuffled sequences hold base composition constant") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(legend.position = "top",
                   plot.background = ggplot2::element_rect(fill = "white",
                                                           colour = NA))
  save_plot(p_sh, "shuffled_vs_background_control", width = 8.5, height = 6)
}

# =============================================================================
# 7. Session log
# =============================================================================
log_session("12_tfbs_robustness")
message("=== 12_tfbs_robustness.R finished ===")
