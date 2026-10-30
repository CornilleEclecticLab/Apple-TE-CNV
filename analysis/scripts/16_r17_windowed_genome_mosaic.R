# =============================================================================
# Script:  16_r17_windowed_genome_mosaic.R
# Purpose: Answer Reviewer 1, comment R1.7. The reviewer asked us to repeat,
#          on the full database, the windowed analysis they had run: map to
#          high-quality reference genomes, call SNPs, hard-filter, count
#          distinct SNPs per window, and plot the result as a heatmap.
#
#          Two complementary heatmaps are produced here.
#
#          Panel A -- local ancestry mosaic. Per domesticated accession and
#          per genomic window, the proportion of M. sieversii ancestry
#          averaged over both haplotypes, from the RFMix v2 run. This is the
#          patchy genome composition the reviewer wanted to visualise, at
#          the resolution RFMix supports.
#
#          Panel B -- windowed variant density. Per accession and per window,
#          the proportion of the available neutral SNPs at which the
#          accession carries a non-reference allele. This is the direct
#          analogue of the reviewer's per-window distinct-SNP count, built
#          on the 1,297 unlinked synonymous SNPs already used throughout the
#          paper.
#
#          LIMIT, stated plainly: the 1,297 SNPs are LD-pruned and genome-wide
#          sparse (about one SNP per 500 kb), so panel B is necessarily coarse
#          and cannot reach the resolution of a re-called, unpruned SNP set.
#          Reproducing the reviewer's analysis at their resolution requires
#          re-mapping the raw reads against the five reference genomes, which
#          is HPC work and is not done here.
#
# Author:  Antho
# Project: TE_Project / final_pipeline (Genome Biology revision)
# =============================================================================

source("scripts/00_setup.R")

library(ggplot2)
library(dplyr)
library(tidyr)

TAG     <- "16_r17_windowed_genome_mosaic"
out_dir <- file.path(dir_outputs, TAG)
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

SRC        <- "~/Documents/TE_Project/1-PCA"
WIN_ANC    <- 200e3   # window for the ancestry mosaic
WIN_SNP    <- 5e6     # window for the variant-density heatmap

# --- Metadata ----------------------------------------------------------------
meta <- read.csv(file.path(SRC, "250923_sample_list_TE.csv"), check.names = FALSE)
names(meta)[1] <- "Sample_ID"   # the file carries a UTF-8 BOM on the first header

# =============================================================================
# PANEL A -- local ancestry mosaic (RFMix)
# =============================================================================
rf <- read.delim(file.path(SRC, "RFMix_out.whole_genome_plot.tsv"))
names(rf)[names(rf) == "strat"] <- "start"

chr_len <- rf %>% group_by(chr) %>% summarise(len = max(end), .groups = "drop")

# Explode each ancestry block into the windows it overlaps, weighting each
# window by the number of base pairs the block contributes to it.
anc <- rf %>%
  mutate(w0 = floor(start / WIN_ANC), w1 = floor((end - 1) / WIN_ANC)) %>%
  rowwise() %>%
  do({
    b <- .
    ws <- b$w0:b$w1
    tibble(chr = b$chr, ind = b$ind, hap = b$hap, pop = b$pop, win = ws,
           bp = pmin(b$end, (ws + 1) * WIN_ANC) - pmax(b$start, ws * WIN_ANC))
  }) %>%
  ungroup()

mosaic <- anc %>%
  group_by(chr, win, ind) %>%
  summarise(frac_sie = sum(bp[pop == "Sie"]) / sum(bp), .groups = "drop") %>%
  mutate(pos_mb = (win * WIN_ANC + WIN_ANC / 2) / 1e6)

# Order accessions by mean sieversii ancestry so the mosaic is readable.
ord <- mosaic %>% group_by(ind) %>%
  summarise(m = mean(frac_sie, na.rm = TRUE), .groups = "drop") %>%
  arrange(m)
mosaic$ind <- factor(mosaic$ind, levels = ord$ind)

write.csv(mosaic, file.path(tab_dir, "ancestry_mosaic_200kb.csv"), row.names = FALSE)

pA <- ggplot(mosaic, aes(pos_mb, ind, fill = frac_sie)) +
  geom_raster() +
  facet_grid(. ~ chr, scales = "free_x", space = "free_x", switch = "x") +
  scale_fill_gradient2(low = "#2c7fb8", mid = "#f7f7f7", high = "#d95f0e",
                       midpoint = 0.5, limits = c(0, 1),
                       name = "proportion\nM. sieversii\nancestry") +
  labs(title = "Local ancestry mosaic across the 71 domesticated accessions",
       subtitle = sprintf("RFMix v2, %.0f kb windows, both haplotypes averaged",
                          WIN_ANC / 1e3),
       x = "Position along chromosome (Mb)", y = NULL) +
  theme_minimal(base_size = 8) +
  theme(axis.text.y = element_text(size = 3),
        axis.text.x = element_text(size = 4),
        panel.spacing = unit(0.06, "lines"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 5, angle = 90),
        strip.placement = "outside",
        legend.key.width = unit(0.3, "cm"))

ggsave(file.path(fig_dir, "FigR1.7A_ancestry_mosaic.pdf"), pA,
       width = 11, height = 6.5, device = cairo_pdf)
ggsave(file.path(fig_dir, "FigR1.7A_ancestry_mosaic.svg"), pA,
       width = 11, height = 6.5)

# =============================================================================
# PANEL B -- windowed variant density from the 1,297 neutral SNPs
# =============================================================================
e <- new.env(); load(file.path(SRC, "dosage_data.rda"), envir = e)
dos <- get("dosage_data", envir = e)
dos <- dos[rownames(dos) %in% meta$Sample_ID, , drop = FALSE]

sp <- do.call(rbind, strsplit(colnames(dos), "_"))
snp <- tibble(snp = colnames(dos), chr = sp[, 1], pos = as.numeric(sp[, 2])) %>%
  mutate(win = floor(pos / WIN_SNP))

long <- as.data.frame(dos) %>%
  tibble::rownames_to_column("Sample_ID") %>%
  pivot_longer(-Sample_ID, names_to = "snp", values_to = "dosage") %>%
  inner_join(snp, by = "snp") %>%
  filter(!is.na(dosage))

dens <- long %>%
  group_by(Sample_ID, chr, win) %>%
  summarise(n_snp = n(), n_var = sum(dosage > 0),
            prop_var = n_var / n_snp, .groups = "drop") %>%
  filter(n_snp >= 3) %>%                       # drop windows too sparse to read
  mutate(pos_mb = (win * WIN_SNP + WIN_SNP / 2) / 1e6) %>%
  left_join(meta[, c("Sample_ID", "Species", "group_structure_CORE_admix_K6")],
            by = "Sample_ID")

ordB <- dens %>% group_by(Sample_ID, Species) %>%
  summarise(m = mean(prop_var), .groups = "drop") %>%
  arrange(Species, m)
dens$Sample_ID <- factor(dens$Sample_ID, levels = ordB$Sample_ID)

write.csv(dens, file.path(tab_dir, "variant_density_5Mb.csv"), row.names = FALSE)

pB <- ggplot(dens, aes(pos_mb, Sample_ID, fill = prop_var)) +
  geom_raster() +
  facet_grid(Species ~ chr, scales = "free", space = "free", switch = "x") +
  scale_fill_viridis_c(option = "magma", name = "proportion of\nneutral SNPs\ncarrying a\nnon-reference\nallele") +
  labs(title = "Windowed variant density across the 189 Malus accessions",
       subtitle = sprintf("1,297 unlinked synonymous SNPs, %.0f Mb windows, windows with < 3 SNPs dropped",
                          WIN_SNP / 1e6),
       x = "Position along chromosome (Mb)", y = NULL) +
  theme_minimal(base_size = 8) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 4),
        panel.spacing = unit(0.06, "lines"),
        panel.grid = element_blank(),
        strip.text.x = element_text(size = 5, angle = 90),
        strip.text.y = element_text(size = 5, angle = 0),
        strip.placement = "outside",
        legend.key.width = unit(0.3, "cm"))

ggsave(file.path(fig_dir, "FigR1.7B_variant_density.pdf"), pB,
       width = 11, height = 7.5, device = cairo_pdf)
ggsave(file.path(fig_dir, "FigR1.7B_variant_density.svg"), pB,
       width = 11, height = 7.5)

# --- Numbers for the response letter ----------------------------------------
summ <- list(
  n_ind_rfmix     = length(unique(mosaic$ind)),
  n_win_anc       = nrow(distinct(mosaic, chr, win)),
  win_anc_kb      = WIN_ANC / 1e3,
  mean_frac_sie   = mean(mosaic$frac_sie, na.rm = TRUE),
  sd_across_win   = sd(tapply(mosaic$frac_sie, paste(mosaic$chr, mosaic$win), mean)),
  n_ind_snp       = length(unique(dens$Sample_ID)),
  n_win_snp       = nrow(distinct(dens, chr, win)),
  win_snp_mb      = WIN_SNP / 1e6,
  median_snp_per_window = median(dens$n_snp)
)
writeLines(paste(names(summ), unlist(summ), sep = " = "),
           file.path(tab_dir, "summary_R1.7.txt"))
print(unlist(summ))

log_session(TAG)
