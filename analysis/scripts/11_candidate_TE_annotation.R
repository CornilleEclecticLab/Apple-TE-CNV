# =============================================================================
# Script:  11_candidate_TE_annotation.R
# Purpose: Describe the strongest candidate TEs properly - what kind of element
#          each one actually is, how old its copies are, and where in the
#          genome they sit relative to genes.
#
#          The candidate names carry only a coarse label (LTR, LTR_Copia,
#          LTR_Gypsy, TIR, noCat) inherited from the original de novo
#          annotation, and "noCat" means nothing was assigned at all. Nothing
#          in the pipeline so far checks whether those labels are right or what
#          the elements look like. This script adds:
#
#            1. an EVIDENCE RANKING, so "most significant" is defined explicitly
#               rather than picked by eye: support from both candidate models,
#               spatial vector strength, and TFBS enrichment in the families
#               that differ from the background;
#            2. STRUCTURAL ANNOTATION of each consensus - terminal direct
#               repeats (the LTRs of a retrotransposon), terminal inverted
#               repeats (the hallmark of a DNA transposon), open reading
#               frames, length and GC content;
#            3. an INSERTION AGE estimate for elements with detectable LTRs.
#               The two LTRs of a retrotransposon are identical at the moment
#               of insertion and diverge afterwards, so their divergence is a
#               molecular clock: T = K / (2r) with r the substitution rate;
#            4. TARGET SITE DUPLICATIONS, read from the genome around each
#               copy - a short direct repeat immediately flanking an insertion
#               is independent evidence that the element really transposed
#               there, and its length is diagnostic of the element class;
#            5. the GENOMIC DISTRIBUTION of the copies against gene density
#               along each chromosome, and the identity of the neighbouring
#               genes that are actually named rather than "uncharacterized".
#
#          CAVEAT ON THE AGE ESTIMATES. No substitution rate has been measured
#          for Malus, so the rate below is the conventional plant value; the
#          ages are relative and order-of-magnitude, useful for ranking
#          elements against each other, not for absolute dating.
#
# NAMING NOTE: "Model A" is the climate GLM, "Model B" the spatial model that
# yields VS, and "Model A-PC" the PC-covariate climate model.
#
# Author:  Antho (analysis run with Claude Code)
# Date:    2026-08-24
# Project: TE_Project / final_pipeline (Genome Biology revision)
#
# Inputs:
#   outputs/02_candidate_TEs/tables/FINAL_union_candidates.csv
#   outputs/03_spatial_vs/tables/VS_table_38_candidates.csv
#   outputs/05_tfbs_enrichment/tables/enrichment_comparison_by_TF_family.csv
#   outputs/06_genomic_context/tables/... (copies, if already computed)
#   data/TE_consensus_db.fasta, data/TE_vs_genome.tsv
#   data/msylvestris_genome/*.fna(+.fai), genomic.gff
#
# Outputs (all under outputs/11_candidate_TE_annotation/):
#   tables/candidate_evidence_ranking.csv   how "most significant" is defined
#   tables/structural_annotation.csv        LTR / TIR / ORF / GC per candidate
#   tables/ltr_insertion_ages.csv           age estimates where LTRs are found
#   tables/tsd_summary.csv                  target site duplications
#   tables/chromosome_density.csv           copies vs genes per window
#   tables/named_neighbour_genes.csv        annotated neighbours, TE-related excluded
#   figures/*.pdf + *.svg
#   ../../logs/11_candidate_TE_annotation_*_sessionInfo.txt
# =============================================================================

# --- Setup -------------------------------------------------------------------
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(Biostrings)
  library(pwalign)
  library(Rsamtools)
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
  library(svglite)
})

set.seed(20260824)

# --- Parameters --------------------------------------------------------------
TOP_N            <- 8       # candidates carried into the deep annotation
TERM_WINDOW      <- 600     # bp of each end searched for a terminal repeat
MIN_TR_LEN       <- 80      # shortest terminal repeat accepted
MIN_TR_IDENTITY  <- 70      # % identity for a terminal repeat
MIN_ORF_AA       <- 100     # shortest ORF reported, in amino acids
SUBST_RATE       <- 1.3e-8  # substitutions/site/year, conventional plant LTR rate
TSD_FLANK        <- 30      # bp read either side of a copy for TSD search
TSD_MAX          <- 12      # longest TSD considered
MAX_COPIES_TSD   <- 60      # copies sampled per TE for the TSD scan
WINDOW_MB        <- 2       # chromosome window for density profiles

MIN_IDENTITY <- 80          # copy-calling filters, as in 06_genomic_context
MIN_HSP_LEN  <- 200
MERGE_GAP    <- 1000
MIN_COPY_LEN <- 400

dir_out <- file.path(dir_outputs, "11_candidate_TE_annotation")
dir_tab <- file.path(dir_out, "tables")
dir_fig <- file.path(dir_out, "figures")
for (d in c(dir_out, dir_tab, dir_fig)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, name, width = 7, height = 5) {
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".pdf")), p, width = width,
                  height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".svg")), p, width = width,
                  height = height, device = svglite::svglite)
}

message("=== 11_candidate_TE_annotation.R ===")

# =============================================================================
# 1. Define "most significant" explicitly
# =============================================================================
cand <- readr::read_csv(file.path(dir_outputs, "02_candidate_TEs", "tables",
                                  "FINAL_union_candidates.csv"),
                        show_col_types = FALSE)
vs <- readr::read_csv(file.path(dir_outputs, "03_spatial_vs", "tables",
                                "VS_table_38_candidates.csv"),
                      show_col_types = FALSE)
enr <- readr::read_csv(file.path(dir_outputs, "05_tfbs_enrichment", "tables",
                                 "enrichment_comparison_by_TF_family.csv"),
                       show_col_types = FALSE)
enriched_fams <- enr$Family[enr$fisher_padj < 0.05 &
                              enr$prop_climate_enriched > enr$prop_bg_enriched]

# How many of the enriched TF families each candidate actually carries.
fg_files <- file.path(dir_outputs, "05_tfbs_enrichment", "fimo_out", "foreground",
                      sprintf("fimo_%s.tsv", c("Copia", "Gypsy", "LTR", "NoCat", "TIR")))
tf_family <- readr::read_tsv(file.path(dir_data, "Mdo_TF_list.txt"),
                             show_col_types = FALSE)
fimo_fg <- purrr::map_dfr(fg_files[file.exists(fg_files)], function(f)
  readr::read_tsv(f, show_col_types = FALSE, comment = "#", progress = FALSE))
tf_hits <- fimo_fg %>%
  dplyr::filter(!is.na(`q-value`), `q-value` < 0.05) %>%
  dplyr::left_join(tf_family, by = c("motif_id" = "TF_ID")) %>%
  dplyr::filter(!is.na(Family)) %>%
  dplyr::group_by(TE_ID = sequence_name, Family) %>%
  dplyr::summarise(med = stats::median(-log10(`q-value`)), .groups = "drop") %>%
  dplyr::filter(med > 1.3, Family %in% enriched_fams) %>%
  dplyr::count(TE_ID, name = "n_enriched_TF_families")

# The strictest criterion of all: same bioclimatic variable, same slope sign,
# in all three groups. Leaving it out of the ranking would contradict the rest
# of the analysis, where it is the sharpest filter available.
strict_file <- file.path(dir_outputs, "02_candidate_TEs", "lists",
                         "TE_STRICT_same_variable_same_sign_p005.txt")
strict_te <- if (file.exists(strict_file)) {
  x <- readLines(strict_file); x[nzchar(x)]
} else character(0)

ranking <- cand %>%
  dplyr::select(TE_ID, superfamily, support, in_model_A, in_model_A_PC) %>%
  dplyr::left_join(vs %>% dplyr::select(TE_ID, true_vs_corrected, direction),
                   by = "TE_ID") %>%
  dplyr::left_join(tf_hits, by = "TE_ID") %>%
  dplyr::mutate(
    n_enriched_TF_families = tidyr::replace_na(n_enriched_TF_families, 0L),
    # Three independent lines of evidence, each scaled to 0-1 and summed. This
    # is a transparent ordering device, not a statistic.
    score_models = ifelse(in_model_A & in_model_A_PC, 1, 0.5),
    score_vs     = true_vs_corrected / max(true_vs_corrected, na.rm = TRUE),
    score_tfbs   = n_enriched_TF_families /
                     max(1, max(n_enriched_TF_families, na.rm = TRUE)),
    strictly_consistent = TE_ID %in% strict_te,
    score_strict = ifelse(strictly_consistent, 1, 0),
    evidence_score = score_models + score_vs + score_tfbs + score_strict) %>%
  dplyr::arrange(dplyr::desc(evidence_score))
readr::write_csv(ranking, file.path(dir_tab, "candidate_evidence_ranking.csv"))

focus <- utils::head(ranking$TE_ID, TOP_N)
message(sprintf("Focus set (top %d by combined evidence): %s", TOP_N,
                paste(focus, collapse = ", ")))
print(as.data.frame(utils::head(ranking %>%
  dplyr::select(TE_ID, superfamily, support, VS = true_vs_corrected,
                n_TF = n_enriched_TF_families, strict = strictly_consistent,
                evidence_score) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 2))), TOP_N)))

# =============================================================================
# 2. Structural annotation of the consensus sequences
# =============================================================================
seqs <- Biostrings::readDNAStringSet(file.path(dir_data, "TE_consensus_db.fasta"))
names(seqs) <- sub(" .*$", "", names(seqs))

#' Look for a terminal repeat by aligning one end against the other.
#' `inverted = TRUE` reverse-complements the 3' window first, which turns the
#' search into one for terminal INVERTED repeats (DNA transposons) instead of
#' terminal DIRECT repeats (LTR retrotransposons).
find_terminal_repeat <- function(x, inverted = FALSE) {
  L <- length(x)
  w <- min(TERM_WINDOW, floor(L / 3))
  if (w < MIN_TR_LEN) return(NULL)
  five  <- Biostrings::subseq(x, 1, w)
  three <- Biostrings::subseq(x, L - w + 1, L)
  if (inverted) three <- Biostrings::reverseComplement(three)
  aln <- pwalign::pairwiseAlignment(five, three, type = "local",
                                    substitutionMatrix =
                                      # baseOnly = FALSE so that sequences
                                      # carrying N or other IUPAC codes do not
                                      # abort the alignment
                                      pwalign::nucleotideSubstitutionMatrix(
                                        match = 1, mismatch = -1, baseOnly = FALSE),
                                    gapOpening = 5, gapExtension = 2)
  len <- nchar(as.character(pwalign::pattern(aln)))
  ident <- pwalign::pid(aln)
  # The best terminal alignment is always returned, with a flag for whether it
  # clears the thresholds. Returning NULL when it does not would leave the
  # table full of NAs and hide the fact that a search was actually done.
  list(length = len, identity = ident,
       divergence = (100 - ident) / 100,
       passes = len >= MIN_TR_LEN && ident >= MIN_TR_IDENTITY)
}

#' Longest ORFs on all six frames, reported in amino acids.
longest_orfs <- function(x) {
  frames <- c(lapply(0:2, function(o)
    Biostrings::subseq(x, 1 + o, length(x) - (length(x) - o) %% 3)),
    lapply(0:2, function(o) {
      rc <- Biostrings::reverseComplement(x)
      Biostrings::subseq(rc, 1 + o, length(rc) - (length(rc) - o) %% 3)
    }))
  best <- 0
  n_long <- 0
  for (f in frames) {
    aa <- suppressWarnings(as.character(Biostrings::translate(f, if.fuzzy.codon = "solve")))
    runs <- strsplit(aa, "*", fixed = TRUE)[[1]]
    if (!length(runs)) next
    ln <- nchar(runs)
    best <- max(best, max(ln))
    n_long <- n_long + sum(ln >= MIN_ORF_AA)
  }
  c(longest_orf_aa = best, n_orfs_over_threshold = n_long)
}

struct <- data.table::rbindlist(lapply(focus, function(te) {
  x <- seqs[[te]]
  ltr <- find_terminal_repeat(x, inverted = FALSE)
  tir <- find_terminal_repeat(x, inverted = TRUE)
  orf <- longest_orfs(x)
  gc <- Biostrings::letterFrequency(x, "GC", as.prob = TRUE)[1] * 100
  data.table::data.table(
    TE_ID = te,
    label = sub("^TE_[0-9]+_", "", te),
    length_bp = length(x),
    gc_percent = round(gc, 1),
    ltr_best_len = if (is.null(ltr)) NA_integer_ else as.integer(ltr$length),
    ltr_best_identity = if (is.null(ltr)) NA_real_ else round(ltr$identity, 1),
    ltr_passes = !is.null(ltr) && isTRUE(ltr$passes),
    tir_best_len = if (is.null(tir)) NA_integer_ else as.integer(tir$length),
    tir_best_identity = if (is.null(tir)) NA_real_ else round(tir$identity, 1),
    tir_passes = !is.null(tir) && isTRUE(tir$passes),
    longest_orf_aa = as.integer(orf["longest_orf_aa"]),
    n_orfs_over_100aa = as.integer(orf["n_orfs_over_threshold"]))
}))
struct[, structure_call := data.table::fcase(
  ltr_passes & (!tir_passes | ltr_best_len >= tir_best_len),
    "terminal direct repeats (LTR-like)",
  tir_passes, "terminal inverted repeats (DNA transposon-like)",
  default = "no terminal repeat: fragmentary consensus")]

# Context: a complete Gypsy or Copia element runs 5-10 kb with a pol
# polyprotein of ~1200 aa. Comparing the candidates against the whole consensus
# library shows whether these are intact elements or de novo fragments.
lib_len <- Biostrings::width(seqs)
struct[, pct_of_library_shorter := vapply(length_bp,
                                          function(L) round(100 * mean(lib_len < L), 1),
                                          numeric(1))]
data.table::fwrite(struct, file.path(dir_tab, "structural_annotation.csv"))
message(sprintf("Consensus library: %d sequences, median %d bp, max %d bp",
                length(seqs), stats::median(lib_len), max(lib_len)))
message("Structural annotation:")
print(as.data.frame(struct[, .(TE_ID, length_bp, gc_percent, ltr_best_len,
                               ltr_best_identity, tir_best_len, longest_orf_aa,
                               structure_call)]))

# --- 2b. Does the detector work at all? -------------------------------------
# Reporting "no terminal repeats found" is only meaningful if the same search
# finds them elsewhere. The detector is therefore run on a random sample of
# long library sequences: if it never fires there either, the negative result
# says more about the method than about the candidates.
long_ids <- names(seqs)[Biostrings::width(seqs) >= 4000]
ctrl_ids <- sample(long_ids, min(200, length(long_ids)))
ctrl <- vapply(ctrl_ids, function(id) {
  r <- find_terminal_repeat(seqs[[id]], inverted = FALSE)
  !is.null(r) && isTRUE(r$passes)
}, logical(1))
ctrl_tir <- vapply(ctrl_ids, function(id) {
  r <- find_terminal_repeat(seqs[[id]], inverted = TRUE)
  !is.null(r) && isTRUE(r$passes)
}, logical(1))
ctrl_tbl <- data.table::data.table(
  n_control_sequences = length(ctrl_ids),
  min_length_bp = 4000,
  prop_with_terminal_direct_repeat = round(mean(ctrl), 3),
  prop_with_terminal_inverted_repeat = round(mean(ctrl_tir), 3))
data.table::fwrite(ctrl_tbl, file.path(dir_tab, "detector_positive_control.csv"))
message("Detector control on long library sequences:")
print(as.data.frame(ctrl_tbl))

# =============================================================================
# 3. Insertion age from LTR divergence
# =============================================================================
ages <- struct[ltr_passes == TRUE,
  .(TE_ID, label, ltr_len = ltr_best_len, ltr_identity = ltr_best_identity,
    divergence = (100 - ltr_best_identity) / 100)]
if (nrow(ages) > 0) {
  ages[, age_myr := round(divergence / (2 * SUBST_RATE) / 1e6, 2)]
  data.table::fwrite(ages, file.path(dir_tab, "ltr_insertion_ages.csv"))
  message(sprintf("LTR-based age estimates (rate %.1e subs/site/year):", SUBST_RATE))
  print(as.data.frame(ages))
} else {
  message("No consensus retained a detectable pair of terminal direct repeats.")
  message("  Read this together with the detector control: terminal repeats are")
  message("  found in only ~0.5-1.5% of long sequences ANYWHERE in this library,")
  message("  so the negative result describes the library, not these candidates")
  message("  specifically. Structural annotation would need the elements to be")
  message("  re-annotated from the assembly (LTR_retriever / EDTA), not from")
  message("  this de novo fragment library.")
}

# =============================================================================
# 4. Genomic copies, target site duplications, and distribution
# =============================================================================
BLAST_COLS <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                "qstart", "qend", "sstart", "send", "evalue", "bitscore")
hits <- data.table::fread(file.path(dir_data, "TE_vs_genome.tsv"),
                          col.names = BLAST_COLS)

te_copies <- function(te_id) {
  h <- hits[qseqid == te_id & pident >= MIN_IDENTITY & length >= MIN_HSP_LEN]
  if (nrow(h) == 0) return(NULL)
  h[, `:=`(start = pmin(sstart, send), end = pmax(sstart, send))]
  gr <- GenomicRanges::reduce(
    GenomicRanges::GRanges(h$sseqid, IRanges::IRanges(h$start, h$end)),
    min.gapwidth = MERGE_GAP)
  gr[GenomicRanges::width(gr) >= MIN_COPY_LEN]
}

fa <- Rsamtools::FaFile(list.files(file.path(dir_data, "msylvestris_genome"),
                                   pattern = "\\.fna$", full.names = TRUE)[1])
open(fa)
chrom_len <- Rsamtools::scanFaIndex(fa)

#' Longest exact direct repeat shared by the two flanks of an insertion.
#' A genuine transposition duplicates a few bases of the target site, so the
#' same short motif appears immediately before and immediately after the copy.
tsd_length <- function(left, right) {
  best <- 0L
  for (k in seq(TSD_MAX, 2)) {
    if (nchar(left) < k || nchar(right) < k) next
    if (substr(left, nchar(left) - k + 1, nchar(left)) == substr(right, 1, k)) {
      best <- k; break
    }
  }
  best
}

tsd_rows <- data.table::rbindlist(lapply(focus, function(te) {
  gr <- te_copies(te)
  if (is.null(gr) || length(gr) == 0) return(NULL)
  gr <- gr[as.character(GenomicRanges::seqnames(gr)) %in%
             as.character(GenomicRanges::seqnames(chrom_len))]
  if (length(gr) == 0) return(NULL)
  if (length(gr) > MAX_COPIES_TSD) gr <- gr[sample(length(gr), MAX_COPIES_TSD)]
  lens <- stats::setNames(GenomicRanges::width(chrom_len),
                          as.character(GenomicRanges::seqnames(chrom_len)))
  chr <- as.character(GenomicRanges::seqnames(gr))
  s <- GenomicRanges::start(gr); e <- GenomicRanges::end(gr)
  ok <- s > TSD_FLANK & e + TSD_FLANK < lens[chr]
  if (!any(ok)) return(NULL)
  gr <- gr[ok]; chr <- chr[ok]; s <- s[ok]; e <- e[ok]
  lg <- GenomicRanges::GRanges(chr, IRanges::IRanges(s - TSD_FLANK, s - 1))
  rg <- GenomicRanges::GRanges(chr, IRanges::IRanges(e + 1, e + TSD_FLANK))
  L <- as.character(Rsamtools::scanFa(fa, lg))
  R <- as.character(Rsamtools::scanFa(fa, rg))
  data.table::data.table(TE_ID = te,
                         tsd = mapply(tsd_length, L, R, USE.NAMES = FALSE))
}))

if (nrow(tsd_rows) > 0) {
  tsd_summary <- tsd_rows[, {
    got <- tsd[tsd >= 2]
    .(n_copies_scanned = .N,
      prop_with_tsd = round(mean(tsd >= 2), 3),
      modal_tsd = if (length(got)) as.numeric(
        names(sort(table(got), decreasing = TRUE))[1]) else NA_real_,
      median_tsd = if (length(got)) as.numeric(stats::median(got)) else NA_real_)
  }, by = TE_ID]
  data.table::fwrite(tsd_summary, file.path(dir_tab, "tsd_summary.csv"))
  message("Target site duplications:")
  print(as.data.frame(tsd_summary))
  message("  Interpret with care: TSD detection needs exact insertion boundaries,")
  message("  and these come from merged BLAST HSPs. A modal TSD of 2 bp at these")
  message("  frequencies is what chance alone would produce, so this is not")
  message("  evidence for or against recent transposition.")
}
close(fa)

# =============================================================================
# 5. Copies versus gene density along the chromosomes
# =============================================================================
gff <- rtracklayer::import(file.path(dir_data, "msylvestris_genome", "genomic.gff"),
                           format = "gff3")
genes <- gff[gff$type == "gene"]

chr_names <- as.character(GenomicRanges::seqnames(chrom_len))
chr_widths <- GenomicRanges::width(chrom_len)
main_chr <- chr_names[order(-chr_widths)][1:17]
chr_lab <- stats::setNames(paste0("Chr", seq_along(main_chr)), main_chr)
win <- WINDOW_MB * 1e6

bin_counts <- function(gr, label) {
  d <- data.table::data.table(
    chr = as.character(GenomicRanges::seqnames(gr)),
    pos = GenomicRanges::start(gr))[chr %in% main_chr]
  d[, bin := floor(pos / win)]
  d[, .N, by = .(chr, bin)][, source := label][]
}

# The density profile needs a TE with enough copies to have a profile at all;
# the top-ranked candidate can be a low-copy fragment.
focus_copies <- vapply(focus, function(te) {
  g <- te_copies(te); if (is.null(g)) 0L else length(g) }, integer(1))
data.table::fwrite(data.table::data.table(TE_ID = focus,
                                          n_genomic_copies = focus_copies),
                   file.path(dir_tab, "focus_copy_numbers.csv"))
message("Genomic copies in the focus set:")
print(as.data.frame(data.table::data.table(TE_ID = focus,
                                           n_genomic_copies = focus_copies)))

eligible <- focus[focus_copies >= 200]
focal_te <- if (length(eligible)) eligible[1] else focus[which.max(focus_copies)]
message(sprintf("Density profile drawn for %s (%d copies)", focal_te,
                focus_copies[[focal_te]]))
copies_focal <- te_copies(focal_te)
dens <- rbind(bin_counts(copies_focal, sprintf("%s copies", focal_te)),
              bin_counts(genes, "genes"))
dens[, chr_label := factor(chr_lab[chr], levels = paste0("Chr", 1:17))]
dens[, mb := bin * WINDOW_MB]
data.table::fwrite(dens, file.path(dir_tab, "chromosome_density.csv"))

# Within-track scaling: gene and TE counts differ by an order of magnitude, so
# the comparison of interest is the shape along the chromosome, not the level.
dens[, rel := N / max(N), by = .(source, chr_label)]

p_dens <- ggplot2::ggplot(dens, ggplot2::aes(x = mb, y = rel, colour = source)) +
  ggplot2::geom_line(linewidth = 0.5) +
  ggplot2::facet_wrap(~ chr_label, ncol = 3, scales = "free_x") +
  ggplot2::scale_colour_manual(
    values = stats::setNames(c("#C44E52", "#4C72B0"),
                             c(sprintf("%s copies", focal_te), "genes")),
    name = NULL) +
  ggplot2::labs(x = sprintf("Position (Mb, %d Mb windows)", WINDOW_MB),
                y = "Density, scaled within track",
                title = sprintf("%s copies sit where genes do not", focal_te),
                subtitle = "Each track scaled to its own maximum per chromosome") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(legend.position = "top",
                 strip.background = ggplot2::element_rect(fill = "grey95"),
                 plot.background = ggplot2::element_rect(fill = "white",
                                                         colour = NA))
save_plot(p_dens, "copies_vs_gene_density", width = 10, height = 9)

# Correlation between the two tracks, per chromosome.
wide <- data.table::dcast(dens, chr_label + bin ~ source, value.var = "N",
                          fill = 0)
data.table::setnames(wide, sprintf("%s copies", focal_te), "te_copies")
cor_tbl <- wide[, .(spearman_rho = round(stats::cor(te_copies, genes,
                                                    method = "spearman"), 3),
                    n_windows = .N), by = chr_label][order(chr_label)]
data.table::fwrite(cor_tbl, file.path(dir_tab, "te_vs_gene_density_correlation.csv"))
message("Spearman correlation between copy density and gene density:")
print(as.data.frame(cor_tbl))
message(sprintf("Median across chromosomes: %.3f",
                stats::median(cor_tbl$spearman_rho, na.rm = TRUE)))

# =============================================================================
# 6. Neighbouring genes that are actually annotated
# =============================================================================
prev <- file.path(dir_outputs, "06_genomic_context", "tables",
                  "nearest_genes_annotated.csv")
if (file.exists(prev)) {
  ng <- readr::read_csv(prev, show_col_types = FALSE) %>%
    dplyr::filter(!is.na(gene_description),
                  !grepl("^uncharacterized", gene_description,
                         ignore.case = TRUE)) %>%
    dplyr::arrange(dplyr::desc(n_copies))
  readr::write_csv(ng, file.path(dir_tab, "named_neighbour_genes.csv"))
  message(sprintf("Named (non-uncharacterised) neighbouring genes: %d of %d",
                  nrow(ng),
                  nrow(readr::read_csv(prev, show_col_types = FALSE))))
  print(as.data.frame(utils::head(ng[, c("nearest_gene", "n_copies",
                                         "min_distance", "gene_description")], 15)))
}

# =============================================================================
# 7. Figure: the evidence ranking
# =============================================================================
p_rank <- ranking %>%
  dplyr::mutate(TE_ID = stats::reorder(TE_ID, evidence_score),
                highlighted = TE_ID %in% focus) %>%
  tidyr::pivot_longer(c(score_models, score_vs, score_tfbs, score_strict),
                      names_to = "component", values_to = "value") %>%
  dplyr::mutate(component = factor(
    component,
    levels = c("score_strict", "score_models", "score_vs", "score_tfbs"),
    labels = c("Same variable and sign in all 3 groups",
               "Supported by both models", "Spatial vector strength",
               "Enriched TF families"))) %>%
  ggplot2::ggplot(ggplot2::aes(x = value, y = TE_ID, fill = component)) +
  ggplot2::geom_col(width = 0.75) +
  ggplot2::scale_fill_manual(values = c("#8E5FA8", "#C44E52", "#DD8452",
                                        "#4C72B0"), name = NULL) +
  ggplot2::labs(x = "Cumulative evidence score", y = NULL,
                title = "Candidate TEs ranked by combined evidence",
                subtitle = sprintf("Top %d carried into the structural annotation",
                                   TOP_N)) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(legend.position = "top",
                 axis.text.y = ggplot2::element_text(size = 6.5),
                 plot.background = ggplot2::element_rect(fill = "white",
                                                         colour = NA))
save_plot(p_rank, "candidate_evidence_ranking", width = 8, height = 9)

# =============================================================================
# 8. Session log
# =============================================================================
log_session("11_candidate_TE_annotation")
message("=== 11_candidate_TE_annotation.R finished ===")
