# =============================================================================
# Script:  06_genomic_context.R
# Purpose: Where does the focal candidate TE actually sit in the genome, and
#          is it close to genes?
#
#          The candidate carried forward from scripts 10 and 11 is a CONSENSUS
#          sequence, not a single locus, so "where does it map" means "where do
#          its genomic copies sit". This script therefore:
#            1. takes the existing BLAST alignment of the TE consensus against
#               the M. sylvestris reference (GCF_916048215.2, drMalSylv7.2),
#            2. merges the high-scoring segment pairs into distinct genomic
#               copies,
#            3. measures each copy's distance to the nearest annotated gene and
#               classifies it (inside a gene / in the 2 kb promoter / nearby /
#               distal),
#            4. tests that proximity against a matched random null, because a
#               gene-dense genome makes almost anything look "near a gene",
#            5. reports which genes those are, and does the same for the other
#               candidate TEs so the focal one has context.
#
#          THE NULL MATTERS. Roughly half the M. sylvestris assembly is within
#          a few kb of some gene, so "241 promoter overlaps" is meaningless on
#          its own. The null here resamples the same number of intervals, with
#          the same width distribution, on the same chromosomes in proportion
#          to their length, and recomputes the same statistics. Only the
#          comparison against that null supports any claim of enrichment.
#
# NAMING NOTE: "Model A" is the climate GLM, "Model B" the spatial model that
# yields VS, and "Model A-PC" the PC-covariate climate model.
#
# Author:  Antho (analysis run with Claude Code)
# Date:    2026-08-24
# Project: TE_Project / review (reviewer response, Genome Biology submission)
#
# Inputs:
#   data/TE_vs_genome.tsv                 BLAST outfmt 6, TE consensus vs genome
#   data/msylvestris_genome/genomic.gff   NCBI Gnomon annotation
#   data/msylvestris_genome/*.fna.fai     chromosome lengths
#   data/TE_consensus_db.fasta            consensus sequences (for lengths)
#   outputs/05_.../lists/TE_STRICT_same_variable_same_sign_p005.txt
#   outputs/05_.../tables/FINAL_union_candidates.csv
#
# Outputs (all under outputs/06_genomic_context/):
#   tables/focal_TE_copies.csv            every merged copy, with nearest gene
#   tables/focal_TE_summary.csv           counts by class, vs the null
#   tables/nearest_genes_annotated.csv    the genes involved, with descriptions
#   tables/all_candidates_gene_proximity.csv  the same stats for all candidates
#   figures/*.pdf + *.svg
#   ../../logs/06_genomic_context_*_sessionInfo.txt
# =============================================================================

# --- Setup -------------------------------------------------------------------
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(svglite)
})

set.seed(20260824)

# --- Parameters --------------------------------------------------------------
MIN_IDENTITY  <- 80      # per-HSP identity floor
MIN_HSP_LEN   <- 200     # per-HSP alignment length floor (bp)
MERGE_GAP     <- 1000    # HSPs closer than this are one copy
MIN_COPY_LEN  <- 400     # a merged copy must span at least this (bp)
PROMOTER_BP   <- 2000    # promoter window used by the original analysis
N_NULL        <- 200     # random resamples for the null

dir_out <- file.path(dir_outputs, "06_genomic_context")
dir_tab <- file.path(dir_out, "tables")
dir_fig <- file.path(dir_out, "figures")
for (d in c(dir_out, dir_tab, dir_fig)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

save_plot <- function(p, name, width = 7, height = 5) {
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".pdf")), p, width = width,
                  height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(dir_fig, paste0(name, ".svg")), p, width = width,
                  height = height, device = svglite::svglite)
}

message("=== 06_genomic_context.R ===")

# =============================================================================
# 1. Which TE, and the reference it was aligned to
# =============================================================================
strict <- readLines(file.path(dir_outputs, "02_candidate_TEs",
                              "lists", "TE_STRICT_same_variable_same_sign_p005.txt"))
FOCAL <- strict[nzchar(strict)][1]
message("Focal TE: ", FOCAL)

candidates <- readr::read_csv(
  file.path(dir_outputs, "02_candidate_TEs", "tables",
            "FINAL_union_candidates.csv"), show_col_types = FALSE)

consensus <- Biostrings::readDNAStringSet(file.path(dir_data, "TE_consensus_db.fasta"))
names(consensus) <- sub(" .*$", "", names(consensus))
focal_len <- as.integer(BiocGenerics::width(consensus[FOCAL]))
message(sprintf("Consensus length: %d bp", focal_len))

# Chromosome lengths, for the null.
fai <- data.table::fread(
  list.files(file.path(dir_data, "msylvestris_genome"), pattern = "\\.fai$",
             full.names = TRUE)[1],
  col.names = c("seqnames", "length", "offset", "linebases", "linewidth"))
chrom_len <- stats::setNames(fai$length, fai$seqnames)

# =============================================================================
# 2. BLAST hits -> distinct genomic copies
# =============================================================================
BLAST_COLS <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                "qstart", "qend", "sstart", "send", "evalue", "bitscore")
hits_all <- data.table::fread(file.path(dir_data, "TE_vs_genome.tsv"),
                              col.names = BLAST_COLS)
message(sprintf("BLAST table: %d HSPs over %d TE consensus sequences",
                nrow(hits_all), data.table::uniqueN(hits_all$qseqid)))

#' Merge the HSPs of one TE into distinct genomic copies.
#' BLAST reports minus-strand hits with sstart > send, so coordinates are
#' normalised before anything else.
te_copies <- function(te_id) {
  h <- hits_all[qseqid == te_id & pident >= MIN_IDENTITY & length >= MIN_HSP_LEN]
  if (nrow(h) == 0) return(NULL)
  h[, `:=`(start = pmin(sstart, send), end = pmax(sstart, send))]
  gr <- GenomicRanges::GRanges(h$sseqid,
                               IRanges::IRanges(h$start, h$end),
                               pident = h$pident, alnlen = h$length)
  red <- GenomicRanges::reduce(gr, min.gapwidth = MERGE_GAP,
                               with.revmap = TRUE)
  rv <- S4Vectors::mcols(red)$revmap
  S4Vectors::mcols(red)$n_hsp <- lengths(rv)
  S4Vectors::mcols(red)$mean_pident <- vapply(
    rv, function(i) mean(h$pident[i]), numeric(1))
  S4Vectors::mcols(red)$aligned_bp <- vapply(
    rv, function(i) sum(h$length[i]), numeric(1))
  S4Vectors::mcols(red)$revmap <- NULL
  red[GenomicRanges::width(red) >= MIN_COPY_LEN]
}

focal_gr <- te_copies(FOCAL)
message(sprintf("%s: %d HSPs pass the filters, merged into %d genomic copies",
                FOCAL,
                nrow(hits_all[qseqid == FOCAL & pident >= MIN_IDENTITY &
                                length >= MIN_HSP_LEN]),
                length(focal_gr)))

# =============================================================================
# 3. Gene annotation
# =============================================================================
gff <- rtracklayer::import(file.path(dir_data, "msylvestris_genome", "genomic.gff"),
                           format = "gff3")
genes <- gff[gff$type == "gene"]
genes <- genes[as.character(GenomicRanges::seqnames(genes)) %in% names(chrom_len)]
message(sprintf("Annotated genes: %d", length(genes)))

# Promoters: PROMOTER_BP upstream of the TSS, strand-aware, as in the original
# analysis.
proms <- GenomicRanges::promoters(genes, upstream = PROMOTER_BP, downstream = 0)

#' Distance of each interval to the nearest gene, plus a positional class.
classify <- function(gr) {
  if (length(gr) == 0) return(NULL)
  common <- intersect(GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(genes))
  gr <- gr[as.character(GenomicRanges::seqnames(gr)) %in% common]
  if (length(gr) == 0) return(NULL)
  d <- GenomicRanges::distanceToNearest(gr, genes, ignore.strand = TRUE)
  out <- data.table::data.table(
    idx = S4Vectors::queryHits(d),
    gene_idx = S4Vectors::subjectHits(d),
    distance = S4Vectors::mcols(d)$distance)
  out[, in_gene := distance == 0]
  ov_prom <- GenomicRanges::countOverlaps(gr, proms, ignore.strand = TRUE) > 0
  out[, in_promoter := ov_prom[idx]]
  out[, class := data.table::fcase(
    in_gene, "inside a gene",
    in_promoter, "in the 2 kb promoter",
    distance <= 5000, "within 5 kb",
    distance <= 20000, "within 20 kb",
    default = "distal (>20 kb)")]
  out[, gr_index := idx]
  out
}

focal_cls <- classify(focal_gr)

copies_tbl <- data.table::data.table(
  TE_ID = FOCAL,
  seqnames = as.character(GenomicRanges::seqnames(focal_gr))[focal_cls$gr_index],
  start = GenomicRanges::start(focal_gr)[focal_cls$gr_index],
  end = GenomicRanges::end(focal_gr)[focal_cls$gr_index],
  copy_width = GenomicRanges::width(focal_gr)[focal_cls$gr_index],
  pct_of_consensus = round(100 * GenomicRanges::width(focal_gr)[focal_cls$gr_index] /
                             focal_len, 1),
  mean_pident = round(S4Vectors::mcols(focal_gr)$mean_pident[focal_cls$gr_index], 2),
  n_hsp = S4Vectors::mcols(focal_gr)$n_hsp[focal_cls$gr_index],
  distance_to_gene = focal_cls$distance,
  class = focal_cls$class,
  nearest_gene = genes$Name[focal_cls$gene_idx],
  gene_id = genes$ID[focal_cls$gene_idx],
  gene_biotype = genes$gene_biotype[focal_cls$gene_idx],
  gene_description = genes$description[focal_cls$gene_idx])
data.table::fwrite(copies_tbl, file.path(dir_tab, "focal_TE_copies.csv"))

message("Positional classes of the copies:")
print(as.data.frame(copies_tbl[, .N, by = class][order(-N)]))
message(sprintf("Median distance to the nearest gene: %s bp",
                format(stats::median(copies_tbl$distance_to_gene), big.mark = " ")))

# =============================================================================
# 4. Matched random null
# =============================================================================
# Intervals of the same widths, on the same chromosomes in proportion to their
# length. Without this, "close to genes" says more about gene density than
# about the TE.
widths <- GenomicRanges::width(focal_gr)
usable <- chrom_len[chrom_len > max(widths) + 1]
probs <- usable / sum(usable)

null_stats <- data.table::rbindlist(lapply(seq_len(N_NULL), function(k) {
  chr <- sample(names(usable), length(widths), replace = TRUE, prob = probs)
  st <- floor(stats::runif(length(widths), 1, usable[chr] - widths))
  gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(st, width = widths))
  cl <- classify(gr)
  if (is.null(cl)) return(NULL)
  data.table::data.table(
    rep = k,
    median_distance = stats::median(cl$distance),
    prop_in_gene = mean(cl$in_gene),
    prop_in_promoter = mean(cl$in_promoter),
    prop_within_5kb = mean(cl$distance <= 5000))
}))

obs <- data.table::data.table(
  median_distance = stats::median(copies_tbl$distance_to_gene),
  prop_in_gene = mean(copies_tbl$class == "inside a gene"),
  prop_in_promoter = mean(copies_tbl$class == "in the 2 kb promoter" |
                            copies_tbl$class == "inside a gene" &
                            FALSE),
  prop_within_5kb = mean(copies_tbl$distance_to_gene <= 5000))
# Promoter overlap counted directly, independent of the class hierarchy.
obs$prop_in_promoter <- mean(focal_cls$in_promoter)

# Both tails are reported: a TE can be closer to genes than chance (enrichment)
# or further from them (depletion), and only looking for enrichment would hide
# the second case entirely.
emp_p <- function(nullv, o, lower) {
  if (lower) (sum(nullv <= o) + 1) / (length(nullv) + 1)
  else (sum(nullv >= o) + 1) / (length(nullv) + 1)
}

summary_tbl <- data.table::data.table(
  statistic = c("median distance to nearest gene (bp)",
                "proportion inside a gene",
                "proportion overlapping a 2 kb promoter",
                "proportion within 5 kb of a gene"),
  observed = c(obs$median_distance, obs$prop_in_gene, obs$prop_in_promoter,
               obs$prop_within_5kb),
  null_median = c(stats::median(null_stats$median_distance),
                  stats::median(null_stats$prop_in_gene),
                  stats::median(null_stats$prop_in_promoter),
                  stats::median(null_stats$prop_within_5kb)),
  null_q025 = c(stats::quantile(null_stats$median_distance, .025),
                stats::quantile(null_stats$prop_in_gene, .025),
                stats::quantile(null_stats$prop_in_promoter, .025),
                stats::quantile(null_stats$prop_within_5kb, .025)),
  null_q975 = c(stats::quantile(null_stats$median_distance, .975),
                stats::quantile(null_stats$prop_in_gene, .975),
                stats::quantile(null_stats$prop_in_promoter, .975),
                stats::quantile(null_stats$prop_within_5kb, .975)),
  # "closer to genes than chance" for the proximity statistics, and "further
  # from genes than chance" for the same statistics read the other way.
  p_closer_than_chance = c(
    emp_p(null_stats$median_distance, obs$median_distance, TRUE),
    emp_p(null_stats$prop_in_gene, obs$prop_in_gene, FALSE),
    emp_p(null_stats$prop_in_promoter, obs$prop_in_promoter, FALSE),
    emp_p(null_stats$prop_within_5kb, obs$prop_within_5kb, FALSE)),
  p_further_than_chance = c(
    emp_p(null_stats$median_distance, obs$median_distance, FALSE),
    emp_p(null_stats$prop_in_gene, obs$prop_in_gene, TRUE),
    emp_p(null_stats$prop_in_promoter, obs$prop_in_promoter, TRUE),
    emp_p(null_stats$prop_within_5kb, obs$prop_within_5kb, TRUE)),
  fold_vs_null = c(
    obs$median_distance / stats::median(null_stats$median_distance),
    obs$prop_in_gene / stats::median(null_stats$prop_in_gene),
    obs$prop_in_promoter / stats::median(null_stats$prop_in_promoter),
    obs$prop_within_5kb / stats::median(null_stats$prop_within_5kb)))
data.table::fwrite(summary_tbl, file.path(dir_tab, "focal_TE_summary.csv"))
message(sprintf("Observed vs %d matched random resamples:", N_NULL))
print(as.data.frame(summary_tbl))

# =============================================================================
# 5. Which genes?
# =============================================================================
near_genes <- copies_tbl[distance_to_gene <= 5000,
  .(n_copies = .N, min_distance = min(distance_to_gene)),
  by = .(gene_id, nearest_gene, gene_biotype, gene_description)][order(-n_copies)]
data.table::fwrite(near_genes, file.path(dir_tab, "nearest_genes_annotated.csv"))
message(sprintf("Distinct genes within 5 kb of a copy: %d", nrow(near_genes)))
print(as.data.frame(utils::head(near_genes[, .(nearest_gene, n_copies,
                                               min_distance,
                                               gene_description)], 12)))

# =============================================================================
# 6. The same statistics for every candidate TE, for context
# =============================================================================
all_te <- sort(unique(candidates$TE_ID))
ctx <- data.table::rbindlist(lapply(all_te, function(te) {
  gr <- te_copies(te)
  if (is.null(gr) || length(gr) == 0) return(NULL)
  cl <- classify(gr)
  if (is.null(cl)) return(NULL)
  data.table::data.table(
    TE_ID = te,
    n_copies = length(gr),
    median_copy_width = stats::median(GenomicRanges::width(gr)),
    median_distance_to_gene = stats::median(cl$distance),
    prop_in_gene = mean(cl$in_gene),
    prop_in_promoter = mean(cl$in_promoter),
    prop_within_5kb = mean(cl$distance <= 5000))
}), fill = TRUE)
ctx[, is_focal := TE_ID == FOCAL]
# The BLAST table predates the revision: it covers the original 34 candidates,
# so the TEs added by the union have no alignment and cannot be profiled here.
not_blasted <- setdiff(all_te, unique(hits_all$qseqid))
if (length(not_blasted) > 0) {
  message(sprintf("NOTE: %d candidate TE(s) absent from the BLAST table and not profiled: %s",
                  length(not_blasted), paste(not_blasted, collapse = ", ")))
  writeLines(not_blasted, file.path(dir_tab, "candidates_not_yet_blasted.txt"))
}
data.table::fwrite(ctx, file.path(dir_tab, "all_candidates_gene_proximity.csv"))
message("Candidate TEs profiled: ", nrow(ctx))
message(sprintf("Focal TE rank by copy number: %d / %d",
                which(order(-ctx$n_copies) == which(ctx$is_focal)), nrow(ctx)))

# =============================================================================
# 7. Figures
# =============================================================================
chr_main <- names(sort(chrom_len, decreasing = TRUE))[1:17]
chr_lab <- stats::setNames(paste0("Chr", seq_along(chr_main)), chr_main)

p_chr <- copies_tbl[seqnames %in% chr_main][
  , .N, by = seqnames][
  , chr := factor(chr_lab[seqnames], levels = paste0("Chr", 1:17))] |>
  ggplot2::ggplot(ggplot2::aes(x = chr, y = N)) +
  ggplot2::geom_col(fill = "#4C72B0", width = 0.75) +
  ggplot2::geom_text(ggplot2::aes(label = N), vjust = -0.3, size = 2.8) +
  ggplot2::labs(x = NULL, y = "Genomic copies",
                title = sprintf("%s copies per chromosome", FOCAL),
                subtitle = sprintf("%d copies passing >=%d%% identity and >=%d bp",
                                   nrow(copies_tbl), MIN_IDENTITY, MIN_COPY_LEN)) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
save_plot(p_chr, "copies_per_chromosome", width = 8, height = 4.5)

p_class <- copies_tbl[, .N, by = class][
  , class := factor(class, levels = c("inside a gene", "in the 2 kb promoter",
                                      "within 5 kb", "within 20 kb",
                                      "distal (>20 kb)"))] |>
  ggplot2::ggplot(ggplot2::aes(x = class, y = N, fill = class)) +
  ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
  ggplot2::geom_text(ggplot2::aes(label = N), vjust = -0.3, size = 3) +
  ggplot2::scale_fill_brewer(palette = "RdYlBu", direction = -1) +
  ggplot2::labs(x = NULL, y = "Genomic copies",
                title = "Position of the copies relative to annotated genes",
                subtitle = sprintf("M. sylvestris GCF_916048215.2, %d annotated genes",
                                   length(genes))) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
save_plot(p_class, "position_relative_to_genes", width = 7, height = 5)

p_null <- ggplot2::ggplot(null_stats, ggplot2::aes(x = prop_within_5kb)) +
  ggplot2::geom_histogram(bins = 25, fill = "grey72", colour = "white") +
  ggplot2::geom_vline(xintercept = obs$prop_within_5kb, colour = "#C44E52",
                      linewidth = 1.1) +
  ggplot2::labs(
    x = "Proportion of copies within 5 kb of a gene",
    y = sprintf("Random resamples (n = %d)", N_NULL),
    title = "Is the TE closer to genes than chance would put it?",
    subtitle = sprintf("Red = observed (%.1f%%); grey = matched random intervals (median %.1f%%)",
                       100 * obs$prop_within_5kb,
                       100 * stats::median(null_stats$prop_within_5kb))) +
  ggplot2::theme_classic(base_size = 10)
save_plot(p_null, "gene_proximity_vs_null", width = 7, height = 5)

p_ctx <- ggplot2::ggplot(ctx, ggplot2::aes(x = n_copies, y = prop_within_5kb)) +
  ggplot2::geom_point(ggplot2::aes(colour = is_focal, size = is_focal),
                      alpha = 0.85) +
  ggplot2::scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#C44E52"),
                               guide = "none") +
  ggplot2::scale_size_manual(values = c(`FALSE` = 1.8, `TRUE` = 3.6),
                             guide = "none") +
  ggplot2::scale_x_log10() +
  ggplot2::labs(x = "Genomic copies (log scale)",
                y = "Proportion within 5 kb of a gene",
                title = "The focal TE among the candidate set",
                subtitle = sprintf("Red = %s", FOCAL)) +
  ggplot2::theme_classic(base_size = 10)
save_plot(p_ctx, "focal_TE_among_candidates", width = 7, height = 5)

# =============================================================================
# 8. Session log
# =============================================================================
log_session("06_genomic_context")
message("=== 06_genomic_context.R finished ===")
