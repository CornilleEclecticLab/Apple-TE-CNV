# =============================================================================
# Script:  14_tfbs_final.R
# Purpose: The definitive TFBS analysis for Fig. 3C and reviewer comment R2.6.
#          It supersedes the enrichment part of 05_tfbs_enrichment.R, which
#          scanned the two groups in separately calibrated FIMO runs.
#
#          THE DESIGN, AND WHY IT IS THIS ONE. FIMO's q-value is an FDR
#          estimated within a run, from that run's p-value distribution. It is
#          therefore not a property of a motif-sequence pair: the same pair
#          gets a different q-value depending on what else was scanned with it.
#          Scanning the 38 candidates in five small runs (4-13 sequences) and
#          the background in five much larger ones (32-104) made the two sides
#          incomparable, and the effect was not subtle - GRAS went from 9/38 to
#          29/38 and bHLH from 9/38 to 0/38 when only the calibration changed.
#
#          Pooling everything into one run fixes that but abandons the
#          per-superfamily structure used throughout the study. The design kept
#          here does both: ONE RUN PER SUPERFAMILY, with the candidates and
#          their background inside the same run. The historical structure is
#          preserved, and within each run the two groups share one calibration,
#          which is all the comparison requires.
#
#          WHAT IT FINDS. No overall difference in motif content, no individual
#          motif enriched (0 of 119 distinct matrices), and one TF family -
#          C3H - marginally enriched (6/38 vs 15/523, BH p = 0.05). That family
#          signal rests on essentially a single position weight matrix
#          (MP00113, assigned by PlantTFDB to two apple paralogues), so it is
#          reported as a lead rather than a conclusion.
#
#          This is a negative result, and it is the direct answer to the
#          reviewer's point: motif content in these TEs reflects what TEs carry
#          generally, not a climate-associated signature.
#
# NAMING NOTE: "Model A" is the climate GLM, "Model B" the spatial model that
# yields VS, and "Model A-PC" the PC-covariate climate model.
#
# Author:  Antho (analysis run with Claude Code)
# Date:    2026-09-01
# Project: TE_Project / final_pipeline (Genome Biology revision)
#
# Inputs:
#   outputs/13_tfbs_single_run/tables/sequence_groups.csv   the 38 + 523
#   outputs/13_tfbs_single_run/fasta/all_sequences.fasta
#   data/Mdo_TF_binding_motifs.meme, data/Mdo_TF_list.txt
#
# Outputs (all under outputs/14_tfbs_final/):
#   fasta/<SF>.fasta, fimo_out/fimo_<SF>.tsv   one run per superfamily
#   tables/enrichment_by_TF_family.csv         the headline test
#   tables/enrichment_by_motif.csv             per-matrix test, PWM-deduplicated
#   tables/overall_comparison.csv              global motif content
#   tables/manuscript_family_verdict.csv       the eight families named in the text
#   figures/Fig3C_final*.pdf|.svg
#   ../../logs/14_tfbs_final_*_sessionInfo.txt
# =============================================================================

# --- Setup -------------------------------------------------------------------
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(Biostrings)
  library(patchwork)
  library(svglite)
})

set.seed(20260824)

# --- Parameters --------------------------------------------------------------
SF            <- c("Copia", "Gypsy", "LTR", "NoCat", "TIR")
Q_THRESHOLD   <- 0.05
ENRICH_THRESH <- 1.3
MEME_BIN      <- path.expand("~/opt/meme/bin")
MANUSCRIPT_FAMS <- c("C2H2", "MYB", "NAC", "WRKY", "ERF", "bZIP",
                     "MIKC_MADS", "M-type_MADS", "BBR-BPC")

dir_out  <- file.path(dir_outputs, "14_tfbs_final")
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

message("=== 14_tfbs_final.R ===")

# =============================================================================
# 1. One FASTA per superfamily, both groups inside
# =============================================================================
groups <- readr::read_csv(file.path(dir_outputs, "13_tfbs_single_run", "tables",
                                    "sequence_groups.csv"),
                          show_col_types = FALSE)
fg_ids <- groups$TE_ID[groups$group == "climate_associated"]
bg_ids <- groups$TE_ID[groups$group == "background"]
message(sprintf("Candidates %d | background %d", length(fg_ids), length(bg_ids)))

seqs <- Biostrings::readDNAStringSet(
  file.path(dir_outputs, "13_tfbs_single_run", "fasta", "all_sequences.fasta"))
names(seqs) <- sub(" .*$", "", names(seqs))

for (s in SF) {
  fa <- file.path(dir_fa, sprintf("%s.fasta", s))
  if (!file.exists(fa)) {
    ids <- groups$TE_ID[groups$superfamily == s]
    Biostrings::writeXStringSet(seqs[ids], fa)
  }
}
message("Per-superfamily composition:")
print(as.data.frame(groups %>% dplyr::count(superfamily, group) %>%
  tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0)))

# =============================================================================
# 2. FIMO, one run per superfamily
# =============================================================================
fimo_bin <- file.path(MEME_BIN, "fimo")
motifs <- file.path(dir_data, "Mdo_TF_binding_motifs.meme")
for (s in SF) {
  tsv <- file.path(dir_fimo, sprintf("fimo_%s.tsv", s))
  if (file.exists(tsv)) next
  stopifnot(file.exists(fimo_bin))
  message("  FIMO on ", s, " (candidates + background together) ...")
  od <- file.path(dir_fimo, sprintf("raw_%s", s))
  system2(fimo_bin, c("--oc", shQuote(od), shQuote(motifs),
                      shQuote(file.path(dir_fa, sprintf("%s.fasta", s)))),
          stdout = NULL, stderr = NULL)
  if (file.exists(file.path(od, "fimo.tsv")))
    file.copy(file.path(od, "fimo.tsv"), tsv, overwrite = TRUE)
}

tf_family <- readr::read_tsv(file.path(dir_data, "Mdo_TF_list.txt"),
                             show_col_types = FALSE)
raw <- purrr::map_dfr(SF, function(s)
  readr::read_tsv(file.path(dir_fimo, sprintf("fimo_%s.tsv", s)),
                  show_col_types = FALSE, comment = "#", progress = FALSE)) %>%
  dplyr::left_join(tf_family, by = c("motif_id" = "TF_ID")) %>%
  dplyr::filter(!is.na(Family), !is.na(`q-value`), `q-value` < Q_THRESHOLD)

message(sprintf("Significant hits with a family annotation: %d | %d genes | %d distinct matrices",
                nrow(raw), dplyr::n_distinct(raw$motif_id),
                dplyr::n_distinct(raw$motif_alt_id)))

# =============================================================================
# 3. Tests
# =============================================================================
# PlantTFDB assigns the same position weight matrix to several paralogues, so
# hits are collapsed to distinct matrices before anything is counted; otherwise
# one matrix is tallied several times and large families are inflated.
pairs <- raw %>%
  dplyr::group_by(TE_ID = sequence_name, Family, PWM = motif_alt_id) %>%
  dplyr::summarise(med = stats::median(-log10(`q-value`)), .groups = "drop") %>%
  dplyr::filter(med > ENRICH_THRESH)

fisher_by <- function(df, key) {
  ks <- sort(unique(df[[key]]))
  purrr::map_dfr(ks, function(k) {
    a <- dplyr::n_distinct(df$TE_ID[df[[key]] == k & df$TE_ID %in% fg_ids])
    b <- dplyr::n_distinct(df$TE_ID[df[[key]] == k & df$TE_ID %in% bg_ids])
    ft <- stats::fisher.test(matrix(c(a, length(fg_ids) - a,
                                      b, length(bg_ids) - b), nrow = 2))
    tibble::tibble(key = k, n_cand = a, n_bg = b,
                   pct_cand = 100 * a / length(fg_ids),
                   pct_bg = 100 * b / length(bg_ids),
                   odds_ratio = unname(ft$estimate), p = ft$p.value)
  }) %>% dplyr::mutate(padj = stats::p.adjust(p, method = "BH"))
}

fam_res <- fisher_by(pairs %>% dplyr::distinct(TE_ID, Family), "Family") %>%
  dplyr::rename(Family = key) %>% dplyr::arrange(padj)
readr::write_csv(fam_res, file.path(dir_tab, "enrichment_by_TF_family.csv"))

motif_res <- fisher_by(pairs %>% dplyr::distinct(TE_ID, PWM), "PWM") %>%
  dplyr::rename(PWM = key) %>% dplyr::arrange(padj)
readr::write_csv(motif_res, file.path(dir_tab, "enrichment_by_motif.csv"))

message(sprintf("TF families enriched (BH < 0.05): %s",
                paste(fam_res$Family[fam_res$padj < 0.05], collapse = ", ")))
message(sprintf("Individual matrices enriched (BH < 0.05): %d of %d",
                sum(motif_res$padj < 0.05), nrow(motif_res)))

# Overall motif content, the question behind "is there any enrichment at all".
per_te <- tibble::tibble(TE_ID = c(fg_ids, bg_ids),
                         group = c(rep("candidates", length(fg_ids)),
                                   rep("background", length(bg_ids)))) %>%
  dplyr::left_join(pairs %>% dplyr::count(TE_ID, name = "n_pairs"),
                   by = "TE_ID") %>%
  dplyr::mutate(n_pairs = tidyr::replace_na(n_pairs, 0L))
w <- stats::wilcox.test(n_pairs ~ group, data = per_te, exact = FALSE)
overall <- tibble::tibble(
  median_candidates = stats::median(per_te$n_pairs[per_te$group == "candidates"]),
  median_background = stats::median(per_te$n_pairs[per_te$group == "background"]),
  wilcoxon_p = w$p.value)
readr::write_csv(overall, file.path(dir_tab, "overall_comparison.csv"))
message("Overall motif content:")
print(as.data.frame(overall))

verdict <- fam_res %>%
  dplyr::filter(Family %in% MANUSCRIPT_FAMS) %>%
  dplyr::mutate(verdict = dplyr::case_when(
    padj < 0.05 & odds_ratio > 1 ~ "SUPPORTED",
    padj < 0.05 & odds_ratio < 1 ~ "REVERSED",
    TRUE ~ "NOT SUPPORTED")) %>%
  dplyr::select(Family, n_cand, n_bg, pct_cand, pct_bg, odds_ratio, padj, verdict)
readr::write_csv(verdict, file.path(dir_tab, "manuscript_family_verdict.csv"))
message("The eight families named in the Discussion:")
print(as.data.frame(verdict %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3)))))

# =============================================================================
# 4. Figure 3C
# =============================================================================
plot_df <- fam_res %>%
  dplyr::filter(n_cand > 0 | n_bg > 0) %>%
  dplyr::mutate(
    signif = padj < 0.05,
    label = ifelse(signif, sprintf("%s  (OR %.1f, p = %.2f)", Family,
                                   odds_ratio, padj), Family)) %>%
  dplyr::arrange(pct_cand - pct_bg)
plot_df$label <- factor(plot_df$label, levels = plot_df$label)

p_main <- ggplot2::ggplot(plot_df) +
  ggplot2::geom_segment(
    ggplot2::aes(y = label, yend = label, x = pct_bg, xend = pct_cand,
                 colour = signif), linewidth = 1.05, lineend = "round") +
  ggplot2::geom_point(ggplot2::aes(y = label, x = pct_bg, shape = "Background"),
                      size = 2.3, colour = "grey30", fill = "white", stroke = 0.7) +
  ggplot2::geom_point(ggplot2::aes(y = label, x = pct_cand, colour = signif,
                                   shape = "Candidates"), size = 2.8) +
  ggplot2::scale_shape_manual(values = c(Background = 21, Candidates = 19),
                              name = NULL) +
  ggplot2::scale_colour_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "grey72"),
                               guide = "none") +
  ggplot2::labs(
    x = "TEs carrying the motif (%)", y = NULL,
    title = "No TF-family signature distinguishes climate-associated TEs",
    subtitle = sprintf(
      "%d candidates vs %d non-climate-associated TEs, scanned together within each superfamily;\noverall motif content does not differ (median %d vs %d families per TE, Wilcoxon p = %.2f)",
      length(fg_ids), length(bg_ids), overall$median_candidates,
      overall$median_background, overall$wilcoxon_p),
    caption = "Fisher's exact test per TF family, Benjamini-Hochberg adjusted. Motifs collapsed to distinct matrices before counting.") +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "top",
                 axis.text.y = ggplot2::element_text(size = 7.5),
                 plot.title = ggplot2::element_text(face = "bold", size = 12),
                 plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30"),
                 plot.caption = ggplot2::element_text(size = 7, colour = "grey45",
                                                      hjust = 0),
                 plot.background = ggplot2::element_rect(fill = "white",
                                                         colour = NA))
save_plot(p_main, "Fig3C_final_comparison", width = 8.5, height = 8)

# Descriptive companion: which motifs the candidates carry at all. This is the
# content of the published panel, kept without the enrichment claim.
top_fams <- fam_res %>% dplyr::filter(n_cand >= 3) %>% dplyr::pull(Family)
heat <- tidyr::expand_grid(TE_ID = fg_ids, Family = top_fams) %>%
  dplyr::left_join(pairs %>% dplyr::distinct(TE_ID, Family) %>%
                     dplyr::mutate(present = TRUE),
                   by = c("TE_ID", "Family")) %>%
  dplyr::mutate(present = tidyr::replace_na(present, FALSE)) %>%
  dplyr::left_join(groups %>% dplyr::select(TE_ID, superfamily), by = "TE_ID")

p_heat <- ggplot2::ggplot(heat, ggplot2::aes(x = TE_ID, y = Family,
                                             fill = present)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
  ggplot2::scale_fill_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "grey92"),
                             labels = c("absent", "present"), name = NULL) +
  ggplot2::facet_grid(. ~ superfamily, scales = "free_x", space = "free_x") +
  ggplot2::labs(x = NULL, y = NULL,
                title = "Motif content of the 38 candidates",
                subtitle = "Descriptive only: none of these families is enriched relative to the background") +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1,
                                        size = 5.5),
    strip.text = ggplot2::element_text(face = "bold", size = 8),
    legend.position = "top",
    plot.title = ggplot2::element_text(face = "bold", size = 11),
    plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30"),
    plot.background = ggplot2::element_rect(fill = "white", colour = NA))
save_plot(p_heat, "Fig3C_final_motif_content", width = 11, height = 4.4)

save_plot(p_main / p_heat + patchwork::plot_layout(heights = c(1.6, 1)),
          "Fig3C_final_composite", width = 11, height = 12)

# =============================================================================
# 5. Session log
# =============================================================================
log_session("14_tfbs_final")
message("=== 14_tfbs_final.R finished ===")
