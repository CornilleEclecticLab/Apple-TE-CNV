# =============================================================================
# Script:  31_figure3_cvd.R
# Purpose: Figure 3 in full (A map, B vector-strength barplot, C CpG
#          composition), redrawn with the colour-blind-safe palette of the
#          paper at the request of Genome Biology. Individual panels plus the
#          assembled figure.
#
#          NO model is re-run: the script reads the same tables as
#          26_figure3_final34.R (arrow_geometry.csv, VS_table_38_candidates.csv,
#          candidate_membership.csv). The plotted values are therefore exactly
#          those of the published figure; only the colours change.
#
#          Three legibility adjustments, all on the form, none on the data:
#            - the very pale blue of Gypsy (#CCE6FF) was indistinguishable from
#              the sea (Delta E 6.2). Every arrow is now backed by a dark
#              casing, and the sea moves from light blue to a neutral blue-grey.
#            - the sampling-group dots take a dark outline instead of a white
#              one, otherwise the light sand disappears on the beige land.
#            - the barplot keeps its reference line at VS 0.5.
#
#          Outputs in SVG + PDF, never PNG.
# Author:  Antho
# =============================================================================
source("scripts/00_setup.R")
suppressMessages({library(dplyr); library(readr); library(tibble); library(ggplot2)
                  library(sf); library(patchwork); library(Biostrings); library(svglite)})

# Output directory. Relative to the pipeline root so the script is portable;
# override with the FIG_OUT environment variable.
OUT <- Sys.getenv("FIG_OUT", unset = file.path("..", "revision_finale", "figures_cvd"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---------- palette of the paper ----------
## One colour per category, the same one in Figures 1, 2 and 3. No green: under
## deuteranopia and protanopia green folds onto yellow and merges with orange.
## Red is kept. Validated with CIEDE2000 on simulated colours; the Fig. 3A-B
## legend has a minimum Delta E of 23.8.
SF_COL <- c(LTR_Copia = "#A50F15",   # rouge fonce   (etait lightgreen)
            LTR_Gypsy = "#CCE6FF",   # bleu tres pale (etait lightblue)
            LTR       = "#E69F00",   # orange        (etait orange)
            noCat     = "#6E6E6E",   # gris moyen    (etait grey)
            TIR       = "#332288")   # indigo        (etait purple)
SF_LAB <- c(LTR_Copia = "Copia", LTR_Gypsy = "Gypsy", LTR = "LTR",
            noCat = "noCat", TIR = "TIR")

## Sampling groups: a legend of their own, chosen to stay apart from each
## other (Delta E 47.8) and from the five arrow colours (17.1).
GRP_COL <- c(SYLW = "#5599FF", SYLE = "#FFD6A5", SYL_Admix = "#000000")

## Panel C: the two TE sets. The pair from the reference swatch sheet,
## Delta E 71.8 under all three kinds of colour blindness, the most separated
## of everything tested. The orange is also the LTR orange of A-B; the two
## legends are separate and spelled out, but if that reuse is a problem,
## #4C72B0 / #C44E52 (Delta E 34.5) does the job.
GRP <- c(background = "#E69F00", climate = "#002255")

LAND   <- "#F4F1EA"
SEA    <- "#CBD5DB"   # etait #DCE9F2, trop proche du bleu pale de Gypsy
CASING <- "#3A3A3A"   # liiere sous les fleches

GROUP_OF <- function(x) dplyr::case_when(
  x == "Q1" ~ "SYLW", x == "Q5" ~ "SYLE",
  x %in% c("Q1_Admix", "Q5_Admix") ~ "SYL_Admix", TRUE ~ NA_character_)

## ---------- data, as they stand ----------
mem  <- read_csv(file.path(dir_outputs, "18_climate_pc_diagnostics", "tables",
                           "candidate_membership.csv"), show_col_types = FALSE)
CAND <- mem$TE_ID[mem$in_modelA]

arr <- read_csv(file.path(dir_outputs, "07_map_figure3A", "tables",
                          "arrow_geometry.csv"), show_col_types = FALSE) %>%
  filter(TE_ID %in% CAND)

sites <- read_csv(file.path(dir_data, "250923_sample_list_TE.csv"),
                  show_col_types = FALSE) %>%
  filter(Species == "M. sylvestris", !is.na(latitude), !is.na(longitude)) %>%
  mutate(group = GROUP_OF(group_structure_CORE_admix_K6)) %>%
  filter(!is.na(group))
sxy <- st_coordinates(st_transform(
  st_as_sf(sites, coords = c("longitude", "latitude"), crs = 4326), 3035))
sites$x <- sxy[, 1]; sites$y <- sxy[, 2]

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(3035)
pad  <- 210e3
xlim <- c(min(sites$x) - pad, max(sites$x) + pad)
ylim <- c(min(sites$y) - pad, max(sites$y) + pad)

## ---------- A: the map ----------
ARROW <- grid::arrow(length = unit(.14, "cm"), type = "closed")
pA <- ggplot() +
  geom_sf(data = world, fill = LAND, colour = "#BFB8AC", linewidth = .25) +
  geom_point(data = sites, aes(x, y, fill = group), shape = 21, size = 1.6,
             colour = "grey20", stroke = .3, alpha = .95) +
  # dark casing: keeps even the very pale Gypsy blue legible over the sea
  geom_segment(data = arr,
               aes(x, y, xend = xend, yend = yend,
                   linewidth = true_vs_corrected),
               colour = CASING, arrow = ARROW, lineend = "round", alpha = .95,
               show.legend = FALSE) +
  geom_segment(data = arr,
               aes(x, y, xend = xend, yend = yend, colour = superfamily,
                   linewidth = true_vs_corrected * .55),
               arrow = ARROW, lineend = "round", alpha = 1) +
  geom_point(data = arr[1, ], aes(x, y), shape = 21, size = 2.2,
             fill = "white", colour = "grey15", stroke = .6) +
  scale_colour_manual(values = SF_COL, labels = SF_LAB, name = "TE superfamily") +
  scale_fill_manual(values = GRP_COL, name = "Sampling group") +
  scale_linewidth_continuous(range = c(.4, 1.6), guide = "none") +
  guides(colour = guide_legend(order = 1, nrow = 2,
                               override.aes = list(linewidth = 1.2, alpha = 1)),
         fill = guide_legend(order = 2, nrow = 1,
                             override.aes = list(size = 2.2))) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = 3035) +
  labs(tag = "A") + theme_void(base_size = 8) +
  theme(panel.background = element_rect(fill = SEA, colour = "grey40"),
        legend.position = "bottom", legend.box = "vertical",
        legend.box.just = "left", legend.spacing.y = unit(.02, "cm"),
        legend.key.size = unit(.26, "cm"),
        legend.text = element_text(size = 6.2),
        legend.title = element_text(size = 6.6),
        plot.tag = element_text(size = 10.5, face = "bold"),
        plot.margin = margin(1, 2, 1, 2))

## ---------- B: spatial vector strength ----------
vs <- read_csv(file.path(dir_outputs, "03_spatial_vs", "tables",
                         "VS_table_38_candidates.csv"), show_col_types = FALSE) %>%
  filter(TE_ID %in% CAND)
## Two label sizes, and the geometry forces it: panel A is a map with a fixed
## aspect ratio (coord_sf), so in the assembled figure patchwork sizes the top
## row from the map, and the name column of panel B gets only 166 px whatever
## the figure height, i.e. 5.03 px per line for 34 names. Past about 4.5 pt the
## names overlap. On its own (3.6 x 4.6 in) the same column gets 8.8 px per
## line and 6 pt is comfortable.
##   - 6 pt   for Fig3B_VSbarplot_CVD, the panel meant for manual assembly
##   - 4.4 pt for the assembled figure
make_pB <- function(label_size) {
  ggplot(vs, aes(reorder(TE_ID, true_vs_corrected), true_vs_corrected,
                 fill = superfamily)) +
    geom_col(width = .78, colour = "grey25", linewidth = .12) +
    geom_hline(yintercept = .5, linetype = 2, colour = "grey40") +
    coord_flip() +
    scale_fill_manual(values = SF_COL, labels = SF_LAB, name = "TE superfamily") +
    labs(x = NULL, y = "Spatial vector strength (VS)", tag = "B") +
    theme_classic(base_size = 8) +
    theme(axis.text.y = element_text(size = label_size),
          axis.ticks.y = element_blank(), legend.position = "none",
          plot.tag = element_text(size = 10.5, face = "bold"),
          plot.margin = margin(1, 2, 1, 2))
}
pB      <- make_pB(6)      # panneau isole
pB_fig  <- make_pB(4.4)    # dans la figure assemblee

## ---------- C: CpG and its two deamination products ----------
stars <- function(p) ifelse(p < .001, "***",
                     ifelse(p < .01, "**", ifelse(p < .05, "*", "n.s.")))

db <- readDNAStringSet(file.path(dir_data, "TE_consensus_db.fasta"))
names(db) <- sub(" .*$", "", names(db))
base <- read_csv(file.path(dir_outputs, "15_sequence_composition", "tables",
                           "composition_per_sequence.csv"), show_col_types = FALSE)
ids  <- base$TE_ID[base$grp == "background" | base$TE_ID %in% CAND]
sq   <- db[intersect(ids, names(db))]
mono <- letterFrequency(sq, c("A", "C", "G", "T"))
di   <- dinucleotideFrequency(sq); Lw <- width(sq)
oe   <- function(d) di[, d] * Lw / (mono[, substr(d, 1, 1)] * mono[, substr(d, 2, 2)])
comp <- tibble(TE_ID = names(sq), CG = oe("CG"), TG = oe("TG"), CA = oe("CA")) %>%
  mutate(grp = factor(ifelse(TE_ID %in% CAND, "climate", "background"),
                      levels = c("background", "climate"))) %>%
  filter(is.finite(CG), is.finite(TG), is.finite(CA))
nb <- sum(comp$grp == "background"); nc <- sum(comp$grp == "climate")
labs2 <- c(sprintf("Non-climate-associated (n = %d)", nb),
           sprintf("Climate-associated (n = %d)", nc))

cg   <- comp %>% mutate(pct = 100 * CG)
medc <- cg %>% group_by(grp) %>% summarise(m = median(pct), .groups = "drop")
w_cg <- wilcox.test(pct ~ grp, cg, exact = FALSE)
dmax <- max(vapply(split(cg$pct, cg$grp),
                   function(x) max(density(x, adjust = 1.1)$y), numeric(1)))
bx <- sort(medc$m); by <- dmax * 1.10; bt <- dmax * .035
pC1 <- ggplot(cg, aes(pct, fill = grp, colour = grp)) +
  geom_density(alpha = .45, linewidth = .7, adjust = 1.1) +
  geom_vline(data = medc, aes(xintercept = m, colour = grp), linetype = 2,
             linewidth = .7, show.legend = FALSE) +
  annotate("segment", x = bx[1], xend = bx[2], y = by, yend = by,
           linewidth = .5, colour = "grey20") +
  annotate("segment", x = bx[1], xend = bx[1], y = by, yend = by - bt,
           linewidth = .5, colour = "grey20") +
  annotate("segment", x = bx[2], xend = bx[2], y = by, yend = by - bt,
           linewidth = .5, colour = "grey20") +
  annotate("text", x = mean(bx), y = by + bt * .35, label = stars(w_cg$p.value),
           size = 4.6, vjust = 0, fontface = "bold", colour = "grey20") +
  expand_limits(y = dmax * 1.24) +
  scale_fill_manual(values = GRP, labels = labs2, name = NULL) +
  scale_colour_manual(values = GRP, labels = labs2, name = NULL) +
  labs(x = "CpG content (% of the amount expected by chance)", y = "Density",
       tag = "C",
       caption = "Wilcoxon rank-sum test:  *** p < 0.001   ** p < 0.01   * p < 0.05   n.s. not significant") +
  theme_classic(base_size = 8) +
  theme(legend.position = "top", legend.key.size = unit(.3, "cm"),
        legend.text = element_text(size = 6.4),
        plot.caption = element_text(hjust = 0, size = 6.4, colour = "grey35"),
        plot.tag = element_text(size = 10.5, face = "bold"),
        plot.margin = margin(1, 2, 1, 2))

trio <- bind_rows(
  comp %>% transmute(grp, pct = 100 * TG, dinuc = "TpG  (deamination product)"),
  comp %>% transmute(grp, pct = 100 * CA, dinuc = "CpA  (deamination product)")) %>%
  mutate(dinuc = factor(dinuc, levels = c("TpG  (deamination product)",
                                          "CpA  (deamination product)")))
medt <- trio %>% group_by(dinuc, grp) %>% summarise(m = median(pct), .groups = "drop")
ann  <- trio %>% group_by(dinuc) %>%
  summarise(p = wilcox.test(pct ~ grp, exact = FALSE)$p.value,
            mc = median(pct[grp == "climate"]),
            mb = median(pct[grp == "background"]), .groups = "drop") %>%
  mutate(lab = sprintf("%.0f%% vs %.0f%%   %s", mc, mb, stars(p)))
pC2 <- ggplot(trio, aes(pct, fill = grp, colour = grp)) +
  geom_density(alpha = .45, linewidth = .6, adjust = 1.1) +
  geom_vline(data = medt, aes(xintercept = m, colour = grp), linetype = 2,
             linewidth = .6, show.legend = FALSE) +
  geom_vline(xintercept = 100, colour = "grey55", linewidth = .4) +
  geom_text(data = ann, aes(x = Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.6, size = 2.6, colour = "grey25") +
  facet_wrap(~ dinuc, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = GRP, guide = "none") +
  scale_colour_manual(values = GRP, guide = "none") +
  labs(x = "Content (% of expected; grey line = expected by chance)", y = "Density") +
  theme_classic(base_size = 8) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0, size = 7),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(1, 2, 1, 2))

save2 <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h,
         device = grDevices::cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".svg")), p, width = w, height = h,
         device = svglite::svglite)
  message("  ecrit : ", file.path(OUT, name), ".pdf / .svg")
}
save2(pA,  "Fig3A_map_CVD",       w = 5.2, h = 4.6)
save2(pB,  "Fig3B_VSbarplot_CVD", w = 3.6, h = 4.6)
save2(pC1, "Fig3C_CpG_CVD",       w = 4.2, h = 3.6)
save2(pC2, "Fig3C_deamination_CVD", w = 4.2, h = 3.6)

## ---------- assembled figure, same layout as the submitted version --------
fig <- (pA | pB_fig) / (pC1 | pC2) +
  plot_layout(heights = c(1.15, .85), widths = c(1.15, 1))
save2(fig, "FIG3_complete_CVD_v3", w = 180 / 25.4, h = 7.6)

message(sprintf("Fig. 3 redessinee sur %d candidats, %d fleches, %d/%d sequences",
                length(CAND), nrow(arr), nc, nb))
message(sprintf("  CpG p = %.3g | TpG p = %.3g | CpA p = %.3g",
                w_cg$p.value, ann$p[1], ann$p[2]))
print(as.data.frame(vs %>% count(superfamily, name = "n")))
