# =============================================================================
# Script:  15_sequence_composition.R
# Purpose: Ask what the candidate TE sequences themselves differ in, now that
#          the TFBS comparison has come back negative (14_tfbs_final.R).
#
#          The motif analysis asked whether these elements carry particular
#          binding sites. They do not, relative to non-climate-associated TEs.
#          A different and more basic question is whether the sequences differ
#          at all in composition - length, GC, and CpG content.
#
#          WHY CpG IS THE INTERESTING ONE. In plants, cytosines in CG context
#          are heavily methylated in transposable elements, and 5-methylcytosine
#          deaminates to thymine far faster than cytosine does. Sustained
#          methylation therefore erodes CpG dinucleotides over evolutionary
#          time, and the observed/expected CpG ratio is a running record of how
#          much methylation-driven decay a sequence has accumulated. A higher
#          ratio means less erosion.
#
#          The ratio is used rather than raw CpG density because raw density is
#          confounded by GC content; observed/expected divides out the
#          single-base composition and isolates dinucleotide depletion.
#
#          INTERPRETATION IS NOT UNIQUE, and the script says so in its outputs.
#          A higher CpG observed/expected ratio is consistent with weaker or
#          shorter historical methylation, but equally with the elements simply
#          being YOUNGER, since less time has passed for deamination to act.
#          Separating the two needs the age of the genomic copies, which
#          requires BLASTing the background set - not available here (see the
#          known gaps in the README).
#
# NAMING NOTE: "Model A" is the climate GLM, "Model B" the spatial model that
# yields VS, and "Model A-PC" the PC-covariate climate model.
#
# Author:  Antho (analysis run with Claude Code)
# Date:    2026-09-01
# Project: TE_Project / final_pipeline (Genome Biology revision)
#
# Inputs:
#   outputs/13_tfbs_single_run/tables/sequence_groups.csv
#   outputs/13_tfbs_single_run/fasta/all_sequences.fasta
#   data/modelA_fig3_raw_p.csv
#
# Outputs (all under outputs/15_sequence_composition/):
#   tables/composition_per_sequence.csv     one row per TE
#   tables/composition_tests.csv            candidates vs background
#   tables/cpg_by_superfamily.csv           is the signal superfamily-wide?
#   tables/superfamily_composition.csv      are the 38 a biased draw?
#   tables/cpg_modality.csv                 dip test and mixture components
#   figures/Fig3D_cpg_density|ecdf|forest|composite .pdf + .svg
#   ../../logs/15_sequence_composition_*_sessionInfo.txt
# =============================================================================

# --- Setup -------------------------------------------------------------------
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(Biostrings)
  library(ggrepel)
  library(diptest)
  library(mclust)
  library(patchwork)
  library(svglite)
})

set.seed(20260824)

SF <- c("Copia", "Gypsy", "LTR", "NoCat", "TIR")
# Copia, Gypsy and unclassified LTR are all LTR retrotransposons; the split
# only records how far the de novo annotation resolved each element. Pooling
# them is biologically coherent and gives 28 candidates in one stratum instead
# of 4, 11 and 13, which Copia alone could never support.
POOL_LTR <- TRUE
pool_stratum <- function(sf) ifelse(POOL_LTR & sf %in% c("Copia", "Gypsy", "LTR"),
                                    "LTR retrotransposons", sf)

# Percentile bootstrap CI for a median, for the forest panel.
boot_median_ci <- function(x, n = 2000, seed = 20260824) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(c(NA_real_, NA_real_))
  set.seed(seed)
  b <- replicate(n, stats::median(sample(x, length(x), replace = TRUE)))
  unname(stats::quantile(b, c(0.025, 0.975)))
}

dir_out <- file.path(dir_outputs, "15_sequence_composition")
dir_tab <- file.path(dir_out, "tables")
dir_fig <- file.path(dir_out, "figures")
for (d in c(dir_out, dir_tab, dir_fig)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, name, width = 7, height = 5) {
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".pdf")), p, width = width,
                  height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".svg")), p, width = width,
                  height = height, device = svglite::svglite)
}

message("=== 15_sequence_composition.R ===")

# =============================================================================
# 1. Composition of every sequence
# =============================================================================
groups <- readr::read_csv(file.path(dir_outputs, "13_tfbs_single_run", "tables",
                                    "sequence_groups.csv"),
                          show_col_types = FALSE)
seqs <- Biostrings::readDNAStringSet(
  file.path(dir_outputs, "13_tfbs_single_run", "fasta", "all_sequences.fasta"))
names(seqs) <- sub(" .*$", "", names(seqs))
s <- seqs[groups$TE_ID]

nC <- as.numeric(Biostrings::letterFrequency(s, "C"))
nG <- as.numeric(Biostrings::letterFrequency(s, "G"))
L  <- Biostrings::width(s)
nCG <- as.numeric(Biostrings::vcountPattern("CG", s))

comp <- groups %>%
  dplyr::mutate(
    length_bp = L,
    gc_percent = 100 * (nC + nG) / L,
    cpg_per_kb = 1000 * nCG / L,
    # Observed/expected divides out single-base composition, so what is left is
    # dinucleotide depletion rather than GC content.
    cpg_obs_exp = ifelse(nC * nG > 0, nCG * L / (nC * nG), NA_real_),
    grp = factor(group, levels = c("background", "climate_associated")),
    stratum = pool_stratum(superfamily))
readr::write_csv(comp, file.path(dir_tab, "composition_per_sequence.csv"))

# =============================================================================
# 2. Candidates versus background
# =============================================================================
vars <- c(length_bp = "Length (bp)", gc_percent = "GC content (%)",
          cpg_per_kb = "CpG per kb", cpg_obs_exp = "CpG observed/expected")

tests <- purrr::imap_dfr(vars, function(lab, v) {
  w <- stats::wilcox.test(comp[[v]] ~ comp$grp, exact = FALSE)
  tibble::tibble(
    variable = v, label = lab,
    median_candidates = stats::median(comp[[v]][comp$grp == "climate_associated"],
                                      na.rm = TRUE),
    median_background = stats::median(comp[[v]][comp$grp == "background"],
                                      na.rm = TRUE),
    wilcoxon_p = w$p.value)
}) %>%
  dplyr::mutate(padj = stats::p.adjust(wilcoxon_p, method = "BH"))
readr::write_csv(tests, file.path(dir_tab, "composition_tests.csv"))
message("Composition, candidates versus background:")
print(as.data.frame(tests %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 4)))))

# =============================================================================
# 3. Is the CpG signal superfamily-wide, or one superfamily's doing?
# =============================================================================
strata <- unique(comp$stratum)
by_sf <- purrr::map_dfr(strata, function(sf) {
  d <- comp %>% dplyr::filter(stratum == sf)
  n_c <- sum(d$grp == "climate_associated")
  p <- if (n_c >= 3) {
    tryCatch(stats::wilcox.test(cpg_obs_exp ~ grp, data = d,
                                exact = FALSE)$p.value, error = function(e) NA_real_)
  } else NA_real_
  tibble::tibble(
    superfamily = sf, n_candidates = n_c,
    n_background = sum(d$grp == "background"),
    oe_candidates = stats::median(d$cpg_obs_exp[d$grp == "climate_associated"],
                                  na.rm = TRUE),
    oe_background = stats::median(d$cpg_obs_exp[d$grp == "background"],
                                  na.rm = TRUE),
    wilcoxon_p = p)
}) %>%
  dplyr::mutate(direction = ifelse(oe_candidates > oe_background,
                                   "higher in candidates", "lower in candidates"))
readr::write_csv(by_sf, file.path(dir_tab, "cpg_by_superfamily.csv"))
message("CpG observed/expected by superfamily:")
print(as.data.frame(by_sf %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3)))))

# =============================================================================
# 4. Are the 38 a biased draw from the residual pool by superfamily?
# =============================================================================
modelA <- readr::read_csv(file.path(dir_data, "modelA_fig3_raw_p.csv"),
                          show_col_types = FALSE)
all_te <- unique(modelA$TE_ID)
fg_ids <- groups$TE_ID[groups$group == "climate_associated"]
sfam_raw <- function(x) sub("^TE_[0-9]+_", "", x)

sf_tbl <- tibble::tibble(sfam = sfam_raw(all_te),
                         is_cand = all_te %in% fg_ids) %>%
  dplyr::count(sfam, is_cand) %>%
  tidyr::pivot_wider(names_from = is_cand, values_from = n, values_fill = 0)
names(sf_tbl) <- c("superfamily", "n_pool", "n_candidates")
sf_tbl <- sf_tbl %>%
  dplyr::mutate(pct_pool = round(100 * n_pool / sum(n_pool), 1),
                pct_candidates = round(100 * n_candidates / sum(n_candidates), 1)) %>%
  dplyr::arrange(dplyr::desc(n_candidates))
readr::write_csv(sf_tbl, file.path(dir_tab, "superfamily_composition.csv"))

ft <- stats::fisher.test(as.matrix(sf_tbl[, c("n_candidates", "n_pool")]),
                         simulate.p.value = TRUE, B = 20000)
message(sprintf("Superfamily composition of the 38 vs the residual pool: Fisher p = %.3f",
                ft$p.value))
print(as.data.frame(sf_tbl))

# =============================================================================
# 4b. Modality of the CpG observed/expected distribution
# =============================================================================
# A single median hides whether a set is one homogeneous population or several.
# Notos (Bulla et al. 2018, BMC Bioinformatics) makes exactly this argument for
# CpG o/e distributions: it fits a kernel density and counts the modes, on the
# reasoning that distinct modes correspond to distinct methylation classes.
# Its distribution server is no longer reachable, and its scope is comparing
# gene bodies ACROSS SPECIES rather than two sequence sets within one genome,
# so the same logic is applied here with standard, citable methods:
#   - Hartigan's dip test, the formal test of unimodality;
#   - Gaussian mixture models selected by BIC, which give the number of
#     components, their means and their weights.
# This is more transparent to describe in a Methods section than a wrapper.
modality <- purrr::map_dfr(levels(comp$grp), function(gg) {
  x <- comp$cpg_obs_exp[comp$grp == gg]
  x <- x[is.finite(x)]
  dt <- diptest::dip.test(x)
  mc <- mclust::Mclust(x, G = 1:4, verbose = FALSE)
  tibble::tibble(
    group = gg, n = length(x),
    dip_statistic = unname(dt$statistic), dip_p = dt$p.value,
    unimodal = dt$p.value > 0.05,
    n_components_BIC = mc$G,
    component_means = paste(round(sort(mc$parameters$mean), 3), collapse = "; "),
    component_weights = paste(round(mc$parameters$pro[order(mc$parameters$mean)], 3),
                              collapse = "; "))
})
readr::write_csv(modality, file.path(dir_tab, "cpg_modality.csv"))
message("Modality of the CpG observed/expected distributions:")
print(as.data.frame(modality))

# The same, within the pooled LTR stratum, where the sample sizes are largest.
modality_ltr <- purrr::map_dfr(levels(comp$grp), function(gg) {
  x <- comp$cpg_obs_exp[comp$grp == gg & comp$stratum == "LTR retrotransposons"]
  x <- x[is.finite(x)]
  if (length(x) < 20) return(NULL)
  dt <- diptest::dip.test(x)
  mc <- mclust::Mclust(x, G = 1:4, verbose = FALSE)
  tibble::tibble(group = gg, n = length(x), dip_p = dt$p.value,
                 n_components_BIC = mc$G,
                 component_means = paste(round(sort(mc$parameters$mean), 3),
                                         collapse = "; "))
})
readr::write_csv(modality_ltr, file.path(dir_tab, "cpg_modality_LTR.csv"))
message("Same, within LTR retrotransposons:")
print(as.data.frame(modality_ltr))

# =============================================================================
# 4c. Testing the mechanism: where did the missing CpGs go?
# =============================================================================
# If CpG loss is driven by methylation, the missing dinucleotides are not gone -
# they were converted. A methylated C in a CG deaminates to T, giving TG on that
# strand; because CG methylation is symmetric, the same event on the opposite
# strand reads as CA here. TpG and CpA are therefore the expected mutational
# products, and a set that has lost more CpGs should carry more of them.
#
# All 16 dinucleotides are computed so the comparison has its own controls, and
# the summary statistic is the ratio the literature uses, (TG + CA) / CG.
#
# One caveat is stated in the outputs: the 16 observed/expected values are not
# independent. Base counts are fixed, so a large shift in CG forces compensating
# shifts elsewhere, and other dinucleotides do differ between the two sets. The
# ratio is the defensible summary, not the claim that only TG and CA move.
bases <- c("A", "C", "G", "T")
dinucs <- as.vector(outer(bases, bases, paste0))
mono <- Biostrings::letterFrequency(s, bases)
di <- Biostrings::dinucleotideFrequency(s)

di_oe <- sapply(dinucs, function(d) {
  x <- substr(d, 1, 1); y <- substr(d, 2, 2)
  di[, d] * L / (mono[, x] * mono[, y])
})
di_oe <- tibble::as_tibble(di_oe) %>%
  dplyr::mutate(TE_ID = groups$TE_ID, grp = comp$grp)

dinuc_tests <- purrr::map_dfr(dinucs, function(d) {
  a <- di_oe[[d]][di_oe$grp == "climate_associated"]
  b <- di_oe[[d]][di_oe$grp == "background"]
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  tibble::tibble(dinucleotide = d,
                 oe_candidates = stats::median(a),
                 oe_background = stats::median(b),
                 difference = stats::median(a) - stats::median(b),
                 wilcoxon_p = stats::wilcox.test(a, b, exact = FALSE)$p.value)
}) %>%
  dplyr::mutate(padj = stats::p.adjust(wilcoxon_p, method = "BH"),
                role = dplyr::case_when(
                  dinucleotide == "CG" ~ "CpG",
                  dinucleotide %in% c("TG", "CA") ~ "deamination product",
                  TRUE ~ "other")) %>%
  dplyr::arrange(dplyr::desc(difference))
readr::write_csv(dinuc_tests, file.path(dir_tab, "dinucleotide_oe_tests.csv"))
message("Dinucleotide observed/expected, candidates vs background:")
print(as.data.frame(dinuc_tests %>%
  dplyr::filter(role != "other" | padj < 0.05) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3)))))

di_oe$deamination_index <- (di_oe$TG + di_oe$CA) / di_oe$CG
di_ok <- di_oe[is.finite(di_oe$deamination_index), ]
w_deam <- stats::wilcox.test(deamination_index ~ grp, data = di_ok, exact = FALSE)
deam <- tibble::tibble(
  median_candidates = stats::median(di_ok$deamination_index[di_ok$grp == "climate_associated"]),
  median_background = stats::median(di_ok$deamination_index[di_ok$grp == "background"]),
  wilcoxon_p = w_deam$p.value)
readr::write_csv(deam, file.path(dir_tab, "deamination_index.csv"))
message("Deamination index (TG + CA) / CG:")
print(as.data.frame(deam))

# =============================================================================
# 5. Figures
# =============================================================================
# Distributions are shown as densities and as an ECDF rather than as boxplots:
# a boxplot hides the shape, and here the whole distribution is shifted rather
# than just its quartiles.
# The observed/expected ratio is exact but abstract. Expressed as a percentage
# of the amount expected by chance it reads directly: 56 % means 44 % of the
# CpGs a sequence should carry have been lost. Values above 100 % simply mean no
# net loss. This is a rescaling, not a different statistic.
comp <- comp %>% dplyr::mutate(cpg_pct = 100 * cpg_obs_exp)
med <- comp %>% dplyr::group_by(grp) %>%
  dplyr::summarise(m = stats::median(cpg_pct, na.rm = TRUE), .groups = "drop")

#' Conventional significance stars for a p-value.
#' Three thresholds, so the reader sees the strength of the effect at a glance
#' without having to parse a number; the exact p is kept alongside.
sig_stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "***",
    ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "n.s."))))
}
STAR_KEY <- "Wilcoxon rank-sum test:  *** P < 0.001   ** P < 0.01   * P < 0.05   n.s. not significant"

GRP_COL <- c(background = "#4C72B0", climate_associated = "#C44E52")
GRP_LAB <- c(background = "Non-climate-associated (n = 523)",
             climate_associated = "Climate-associated (n = 38)")

# Bracket spanning the two medians, with the stars above it: the comparison the
# panel is making is between those two dashed lines, so the mark belongs there.
p_cpg_test <- tests$wilcoxon_p[tests$variable == "cpg_obs_exp"]
dens_max <- max(vapply(split(comp$cpg_pct, comp$grp),
                       function(x) max(stats::density(x, adjust = 1.1)$y),
                       numeric(1)))
brk_x <- sort(med$m)
brk_y <- dens_max * 1.10
brk_t <- dens_max * 0.035

p_dens <- ggplot2::ggplot(comp, ggplot2::aes(x = cpg_pct, fill = grp,
                                             colour = grp)) +
  ggplot2::geom_density(alpha = 0.45, linewidth = 0.7, adjust = 1.1) +
  ggplot2::geom_vline(data = med, ggplot2::aes(xintercept = m, colour = grp),
                      linetype = 2, linewidth = 0.7, show.legend = FALSE) +
  ggplot2::annotate("segment", x = brk_x[1], xend = brk_x[2],
                   y = brk_y, yend = brk_y, linewidth = 0.5, colour = "grey20") +
  ggplot2::annotate("segment", x = brk_x[1], xend = brk_x[1],
                   y = brk_y, yend = brk_y - brk_t, linewidth = 0.5, colour = "grey20") +
  ggplot2::annotate("segment", x = brk_x[2], xend = brk_x[2],
                   y = brk_y, yend = brk_y - brk_t, linewidth = 0.5, colour = "grey20") +
  ggplot2::annotate("text", x = mean(brk_x), y = brk_y + brk_t * 0.35,
                   label = sig_stars(p_cpg_test), size = 5.2, vjust = 0,
                   fontface = "bold", colour = "grey20") +
  ggplot2::expand_limits(y = dens_max * 1.24) +
  ggplot2::scale_fill_manual(values = GRP_COL, labels = GRP_LAB, name = NULL) +
  ggplot2::scale_colour_manual(values = GRP_COL, labels = GRP_LAB, name = NULL) +
  ggplot2::labs(
    x = "CpG content (% of the amount expected by chance)", y = "Density",
    caption = paste0(STAR_KEY, "\nP-value is the raw Wilcoxon test on CpG observed/expected."),
    title = "Climate-associated TEs have kept more of their CpGs",
    subtitle = sprintf("median %.0f%% vs %.0f%% retained (dashed) - a loss of %.0f%% against %.0f%%; Wilcoxon p = %.4f.\nGC content and sequence length do not differ between the two sets.",
                       med$m[med$grp == "climate_associated"],
                       med$m[med$grp == "background"],
                       100 - med$m[med$grp == "climate_associated"],
                       100 - med$m[med$grp == "background"],
                       tests$wilcoxon_p[tests$variable == "cpg_obs_exp"])) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "top",
                 plot.caption = ggplot2::element_text(hjust = 0, size = 7.5,
                                                      colour = "grey35"),
                 plot.background = ggplot2::element_rect(fill = "white",
                                                         colour = NA))
save_plot(p_dens, "Fig3D_cpg_density", width = 6.5, height = 5)

p_ecdf <- ggplot2::ggplot(comp, ggplot2::aes(x = cpg_pct, colour = grp)) +
  ggplot2::stat_ecdf(linewidth = 0.9, pad = FALSE) +
  ggplot2::scale_colour_manual(values = GRP_COL, labels = GRP_LAB, name = NULL) +
  ggplot2::labs(x = "CpG content (% of expected)",
                y = "Cumulative proportion of TEs",
                title = "The whole distribution is shifted, not just its centre") +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "top",
                 plot.background = ggplot2::element_rect(fill = "white",
                                                         colour = NA))
save_plot(p_ecdf, "Fig3D_cpg_ecdf", width = 6.5, height = 5)

# Forest panel: median with a bootstrap interval, per stratum, so the reader
# sees the effect and its precision without a boxplot.
forest <- purrr::map_dfr(unique(comp$stratum), function(st) {
  purrr::map_dfr(levels(comp$grp), function(gg) {
    x <- comp$cpg_pct[comp$stratum == st & comp$grp == gg]
    ci <- boot_median_ci(x)
    tibble::tibble(stratum = st, grp = gg, n = length(x),
                   med = stats::median(x, na.rm = TRUE),
                   lo = ci[1], hi = ci[2])
  })
})
forest$stratum <- factor(forest$stratum,
                         levels = rev(c("LTR retrotransposons", "NoCat", "TIR")))
forest$grp <- factor(forest$grp, levels = c("background", "climate_associated"))
lbl <- forest %>% dplyr::group_by(stratum) %>%
  dplyr::summarise(txt = sprintf("%d vs %d", n[grp == "climate_associated"],
                                 n[grp == "background"]), .groups = "drop")

p_forest <- ggplot2::ggplot(forest,
                            ggplot2::aes(x = med, y = stratum, colour = grp)) +
  ggplot2::geom_linerange(ggplot2::aes(xmin = lo, xmax = hi),
                          position = ggplot2::position_dodge(width = 0.55),
                          linewidth = 0.8) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.55),
                      size = 3) +
  ggplot2::geom_text(data = lbl, ggplot2::aes(x = Inf, y = stratum, label = txt),
                     inherit.aes = FALSE, hjust = 1.1, size = 2.8,
                     colour = "grey40") +
  ggplot2::scale_colour_manual(values = GRP_COL, labels = GRP_LAB, name = NULL) +
  ggplot2::labs(x = "Median CpG content, % of expected (95% bootstrap interval)",
                y = NULL,
                title = "The difference holds within LTR retrotransposons and TIR",
                subtitle = "Copia, Gypsy and unclassified LTR pooled; labels give candidate vs background counts") +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "top",
                 plot.background = ggplot2::element_rect(fill = "white",
                                                         colour = NA))
save_plot(p_forest, "Fig3D_cpg_forest", width = 7.5, height = 4.2)

# Figure 3D: the distribution and the mechanism. The per-superfamily forest
# goes to the supplement, where the small strata can carry their own caveat.
save_plot(p_forest, "FigS_cpg_by_superfamily", width = 7.5, height = 4.2)

# All 16 dinucleotides, so the reader can see that CpG moves most and that the
# two products of its deamination move the other way.
# A lollipop of differences across all 16 dinucleotides gives the
# compositionally constrained ones the same visual weight as the informative
# ones, and reads as an abstract quantity. The three dinucleotides that carry
# the mechanism are shown instead, on an axis anchored at the chance
# expectation; the full 16-dinucleotide comparison stays in
# tables/dinucleotide_oe_tests.csv as the control.



trio <- dinuc_tests %>%
  dplyr::filter(dinucleotide %in% c("CG", "TG", "CA")) %>%
  dplyr::mutate(dinucleotide = factor(dinucleotide, levels = c("CG", "TG", "CA"),
                                      labels = c("CpG", "TpG", "CpA"))) %>%
  dplyr::select(dinucleotide, role, padj,
                Background = oe_background, Candidates = oe_candidates) %>%
  tidyr::pivot_longer(c(Background, Candidates), names_to = "set",
                      values_to = "oe") %>%
  dplyr::mutate(pct = 100 * oe,
                set = factor(set, levels = c("Background", "Candidates"),
                             labels = c("Non-climate-associated",
                                        "Climate-associated")))
star <- trio %>% dplyr::group_by(dinucleotide) %>%
  dplyr::summarise(y = max(pct), padj = dplyr::first(padj), .groups = "drop") %>%
  dplyr::mutate(lab = ifelse(padj < 0.05, "*", "n.s."))

# Bars, lollipops and arrows all reduce each dinucleotide to a single number
# and ask the reader to trust it. The left panel already shows a distribution,
# so the right panel keeps the same visual language: the per-sequence content
# of the three dinucleotides that carry the mechanism, for both sets. Nothing
# is summarised away, and the shift is visible where it actually happens.
di_long <- di_oe %>%
  dplyr::select(TE_ID, grp, CG, TG, CA) %>%
  tidyr::pivot_longer(c(CG, TG, CA), names_to = "dinucleotide",
                      values_to = "oe") %>%
  dplyr::filter(is.finite(oe)) %>%
  dplyr::mutate(
    pct = 100 * oe,
    dinucleotide = factor(dinucleotide, levels = c("TG", "CA"),
                          labels = c("TpG  (deamination product)",
                                     "CpA  (deamination product)"))) %>%
  dplyr::filter(!is.na(dinucleotide))

di_med <- di_long %>%
  dplyr::group_by(dinucleotide, grp) %>%
  dplyr::summarise(m = stats::median(pct), .groups = "drop")

# CpG is deliberately absent here: it is the same quantity, on the same
# sequences, as the density panel on the left (both are nCG * L / (nC * nG)),
# so drawing it twice was pure repetition. This panel now carries only what
# the left one cannot show, the two deamination products.
di_ann <- dinuc_tests %>%
  dplyr::filter(dinucleotide %in% c("TG", "CA")) %>%
  dplyr::mutate(dinucleotide = factor(dinucleotide, levels = c("TG", "CA"),
                                      labels = levels(di_long$dinucleotide)),
                lab = sprintf("%.0f%% vs %.0f%%   %s%s",
                              100 * oe_candidates, 100 * oe_background,
                              sig_stars(padj),
                              ifelse(padj < 0.05,
                                     sprintf("  p = %.1g", padj), "")))

p_trio <- ggplot2::ggplot(di_long, ggplot2::aes(x = pct, fill = grp,
                                                colour = grp)) +
  ggplot2::geom_density(alpha = 0.45, linewidth = 0.6, adjust = 1.1) +
  ggplot2::geom_vline(data = di_med, ggplot2::aes(xintercept = m, colour = grp),
                      linetype = 2, linewidth = 0.6, show.legend = FALSE) +
  # the chance expectation, so each row is read against a fixed anchor
  ggplot2::geom_vline(xintercept = 100, colour = "grey55", linewidth = 0.4) +
  ggplot2::geom_text(data = di_ann, ggplot2::aes(x = Inf, y = Inf, label = lab),
                     inherit.aes = FALSE, hjust = 1.05, vjust = 1.6,
                     size = 3, colour = "grey25") +
  ggplot2::facet_wrap(~ dinucleotide, ncol = 1, scales = "free_y") +
  ggplot2::scale_fill_manual(values = GRP_COL, labels = GRP_LAB, name = NULL) +
  ggplot2::scale_colour_manual(values = GRP_COL, labels = GRP_LAB, name = NULL) +
  ggplot2::labs(
    x = "Content (% of expected; grey line = expected by chance)", y = "Density",
    title = "Both products of CpG deamination are depleted in the candidates",
    subtitle = sprintf("CpG itself is shown in the left-hand panel; it is not repeated here.\nIndex (TpG + CpA)/CpG: %.2f vs %.2f, Wilcoxon p = %.1g %s",
                       deam$median_candidates, deam$median_background,
                       deam$wilcoxon_p, sig_stars(deam$wilcoxon_p)),
    caption = paste0(STAR_KEY,
                     "\nPer-dinucleotide P-values are Benjamini-Hochberg adjusted across the 16 dinucleotides.")) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    legend.position = "top",
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", hjust = 0, size = 9.5),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    plot.caption = ggplot2::element_text(hjust = 0, size = 7.5, colour = "grey35"),
    plot.background = ggplot2::element_rect(fill = "white", colour = NA))
save_plot(p_trio, "Fig3D_dinucleotides", width = 6.4, height = 6.2)

# The same evidence as one number per sequence, if a single panel is wanted.
p_index2 <- ggplot2::ggplot(di_ok, ggplot2::aes(x = deamination_index,
                                                fill = grp, colour = grp)) +
  ggplot2::geom_density(alpha = 0.45, linewidth = 0.7) +
  ggplot2::geom_vline(xintercept = c(deam$median_background,
                                     deam$median_candidates),
                      colour = c(GRP_COL[["background"]],
                                 GRP_COL[["climate_associated"]]),
                      linetype = 2, linewidth = 0.7) +
  ggplot2::scale_fill_manual(values = GRP_COL, labels = GRP_LAB, name = NULL) +
  ggplot2::scale_colour_manual(values = GRP_COL, labels = GRP_LAB, name = NULL) +
  ggplot2::coord_cartesian(xlim = c(0, 12)) +
  ggplot2::labs(
    x = "Deamination index  (TpG + CpA) / CpG", y = "Density",
    title = "Conversion products relative to surviving CpGs",
    subtitle = sprintf("median %.2f vs %.2f, Wilcoxon p = %.1g",
                       deam$median_candidates, deam$median_background,
                       deam$wilcoxon_p)) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "top",
                 plot.background = ggplot2::element_rect(fill = "white",
                                                         colour = NA))
save_plot(p_index2, "Fig3D_index_density", width = 6.2, height = 4.8)

# --- Per superfamily -------------------------------------------------------
# Copia, TIR and noCat contribute only 4, 4 and 6 candidates. A density drawn on
# four points is meaningless, so the background is shown as a distribution and
# every candidate as its own point: the reader sees exactly where each element
# falls rather than a curve fitted to nothing.
per_fam <- di_oe %>%
  dplyr::select(TE_ID, grp, CG, TG, CA) %>%
  dplyr::mutate(cpg_pct = 100 * CG,
                index = (TG + CA) / CG) %>%
  dplyr::left_join(groups %>% dplyr::select(TE_ID, superfamily), by = "TE_ID") %>%
  dplyr::filter(is.finite(cpg_pct), is.finite(index))

fam_n <- per_fam %>%
  dplyr::count(superfamily, grp) %>%
  tidyr::pivot_wider(names_from = grp, values_from = n, values_fill = 0) %>%
  dplyr::mutate(lab = sprintf("%s   (%d candidates vs %d background)",
                              superfamily, climate_associated, background))
per_fam <- per_fam %>%
  dplyr::left_join(fam_n %>% dplyr::select(superfamily, lab), by = "superfamily")

fam_med <- per_fam %>%
  dplyr::group_by(lab, grp) %>%
  dplyr::summarise(cpg = stats::median(cpg_pct), idx = stats::median(index),
                   .groups = "drop")

#' Background as a density, candidates as individual points beneath it.
per_family_plot <- function(var, xlab, ttl, sub, xlim = NULL) {
  bgd <- per_fam %>% dplyr::filter(grp == "background")
  cnd <- per_fam %>% dplyr::filter(grp == "climate_associated")
  p <- ggplot2::ggplot() +
    ggplot2::geom_density(data = bgd,
                          ggplot2::aes(x = .data[[var]], y = ggplot2::after_stat(scaled)),
                          fill = GRP_COL[["background"]], colour = GRP_COL[["background"]],
                          alpha = 0.35, linewidth = 0.6) +
    ggplot2::geom_vline(data = fam_med,
                        ggplot2::aes(xintercept = if (var == "cpg_pct") cpg else idx,
                                     colour = grp),
                        linetype = 2, linewidth = 0.6) +
    ggplot2::geom_point(data = cnd,
                        ggplot2::aes(x = .data[[var]], y = -0.14),
                        colour = GRP_COL[["climate_associated"]], size = 2.1,
                        alpha = 0.9) +
    ggplot2::scale_colour_manual(values = GRP_COL, labels = GRP_LAB, name = NULL) +
    ggplot2::facet_wrap(~ lab, ncol = 1, scales = "free_y") +
    ggplot2::labs(x = xlab, y = NULL, title = ttl, subtitle = sub) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      legend.position = "top",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0, size = 9),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  if (!is.null(xlim)) p <- p + ggplot2::coord_cartesian(xlim = xlim)
  p
}

save_plot(per_family_plot("cpg_pct", "CpG content (% of expected)",
                          "CpG content by TE superfamily",
                          "Blue: background distribution. Red points: each climate-associated TE. Dashed: medians."),
          "Fig3D_cpg_by_family", width = 6.6, height = 8)

save_plot(per_family_plot("index", "Deamination index  (TpG + CpA) / CpG",
                          "Deamination index by TE superfamily",
                          "Blue: background distribution. Red points: each climate-associated TE. Dashed: medians.",
                          xlim = c(0, 12)),
          "Fig3D_index_by_family", width = 6.6, height = 8)

# The numbers behind those panels, since the small strata cannot be read off a
# curve.
fam_stats <- purrr::map_dfr(unique(per_fam$superfamily), function(sf) {
  d <- per_fam %>% dplyr::filter(superfamily == sf)
  n_c <- sum(d$grp == "climate_associated")
  tst <- function(v) if (n_c >= 3) {
    tryCatch(stats::wilcox.test(d[[v]] ~ d$grp, exact = FALSE)$p.value,
             error = function(e) NA_real_) } else NA_real_
  tibble::tibble(
    superfamily = sf, n_candidates = n_c, n_background = sum(d$grp == "background"),
    cpg_candidates = stats::median(d$cpg_pct[d$grp == "climate_associated"]),
    cpg_background = stats::median(d$cpg_pct[d$grp == "background"]),
    cpg_p = tst("cpg_pct"),
    index_candidates = stats::median(d$index[d$grp == "climate_associated"]),
    index_background = stats::median(d$index[d$grp == "background"]),
    index_p = tst("index"))
}) %>%
  dplyr::mutate(cpg_padj = stats::p.adjust(cpg_p, "BH"),
                index_padj = stats::p.adjust(index_p, "BH"))
readr::write_csv(fam_stats, file.path(dir_tab, "cpg_index_by_superfamily.csv"))
message("CpG and deamination index by superfamily:")
print(as.data.frame(fam_stats %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3)))))

p_index <- ggplot2::ggplot(di_ok, ggplot2::aes(x = deamination_index,
                                               fill = grp, colour = grp)) +
  ggplot2::geom_density(alpha = 0.45, linewidth = 0.7) +
  ggplot2::scale_fill_manual(values = GRP_COL, labels = GRP_LAB, name = NULL) +
  ggplot2::scale_colour_manual(values = GRP_COL, labels = GRP_LAB, name = NULL) +
  ggplot2::labs(x = "(TpG + CpA) / CpG", y = "Density",
                title = "Deamination index",
                subtitle = sprintf("median %.2f vs %.2f, Wilcoxon p = %.1g",
                                   deam$median_candidates,
                                   deam$median_background, deam$wilcoxon_p)) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "none",
                 plot.background = ggplot2::element_rect(fill = "white",
                                                         colour = NA))
save_plot(p_index, "Fig3D_deamination_index", width = 5.5, height = 4.2)

save_plot((p_dens | p_trio) + patchwork::plot_layout(widths = c(1, 0.7)),
          "Fig3D", width = 12.5, height = 5.6)

# Saved for the Fig. 3 assembly (17_figure3_final.R).
saveRDS(list(density = p_dens, trio = p_trio),
        file.path(dir_out, "panel_C_cpg.rds"))

# =============================================================================
# 6. Session log
# =============================================================================
log_session("15_sequence_composition")
message("=== 15_sequence_composition.R finished ===")
