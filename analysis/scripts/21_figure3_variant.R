# =============================================================================
# Script:  21_figure3_variant.R
# Purpose: Rebuild Figure 3 for an arbitrary candidate set, so the effect of
#          changing the definition of "climate-associated" can be seen rather
#          than argued about. Set CAND below.
# Author:  Antho
# =============================================================================
source("scripts/00_setup.R")
library(dplyr); library(readr); library(ggplot2); library(patchwork); library(sf)

TAG <- "21_figure3_variant"
dir_out <- file.path(dir_outputs, TAG); dir_fig <- file.path(dir_out, "figures")
dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)

mem <- read_csv(file.path(dir_outputs, "18_climate_pc_diagnostics", "tables",
                          "candidate_membership.csv"), show_col_types = FALSE)
CAND  <- mem$TE_ID[mem$in_modelA]        # <- the 34 replicated candidates
LABEL <- sprintf("%d replicated candidates", length(CAND))

SF_COL <- c(LTR_Copia = "lightgreen", LTR_Gypsy = "lightblue", LTR = "orange",
            noCat = "grey", TIR = "purple")
SF_LAB <- c(LTR_Copia = "Copia", LTR_Gypsy = "Gypsy", LTR = "LTR",
            noCat = "noCat", TIR = "TIR")
GRP_COL <- c(background = "#4C72B0", climate_associated = "#C44E52")

# --- A : map -----------------------------------------------------------------
arr <- read_csv(file.path(dir_outputs, "07_map_figure3A", "tables",
                          "arrow_geometry.csv"), show_col_types = FALSE) %>%
  filter(TE_ID %in% CAND)
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(3035)
sites <- read_csv(file.path(dir_data, "250923_sample_list_TE.csv"), show_col_types = FALSE) %>%
  filter(Species == "M. sylvestris", !is.na(latitude), !is.na(longitude))
sxy <- st_coordinates(st_transform(st_as_sf(sites, coords = c("longitude","latitude"),
                                            crs = 4326), 3035))
sites$x <- sxy[,1]; sites$y <- sxy[,2]
pad <- 210e3
xlim <- c(min(sites$x)-pad, max(sites$x)+pad); ylim <- c(min(sites$y)-pad, max(sites$y)+pad)

pA <- ggplot() +
  geom_sf(data = world, fill = "#F4F1EA", colour = "#BFB8AC", linewidth = 0.25) +
  geom_point(data = sites, aes(x, y), colour = "grey45", size = 1.1, alpha = .7) +
  geom_segment(data = arr, aes(x = x, y = y, xend = xend, yend = yend,
                               colour = superfamily, linewidth = true_vs_corrected),
               arrow = grid::arrow(length = unit(0.16,"cm"), type = "closed"),
               lineend = "round", alpha = .95) +
  geom_point(data = arr[1,], aes(x, y), shape = 21, size = 2.4,
             fill = "white", colour = "grey15", stroke = .6) +
  scale_colour_manual(values = SF_COL, labels = SF_LAB, name = "TE superfamily") +
  scale_linewidth_continuous(range = c(.45,1.6), guide = "none") +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = 3035) +
  labs(tag = "A") +
  theme_void(base_size = 8) +
  theme(panel.background = element_rect(fill = "#DCE9F2", colour = "grey40"),
        legend.position = "bottom", legend.key.size = unit(.28,"cm"),
        plot.tag = element_text(size = 10.5, face = "bold"))

# --- B : vector strength ------------------------------------------------------
vs <- read_csv(file.path(dir_outputs, "03_spatial_vs", "tables",
                         "VS_table_38_candidates.csv"), show_col_types = FALSE) %>%
  filter(TE_ID %in% CAND)
pB <- ggplot(vs, aes(reorder(TE_ID, true_vs_corrected), true_vs_corrected,
                     fill = superfamily)) +
  geom_col(width = .7) +
  geom_hline(yintercept = .5, linetype = 2, colour = "grey40") +
  coord_flip() +
  scale_fill_manual(values = SF_COL, labels = SF_LAB, guide = "none") +
  labs(x = NULL, y = "Spatial vector strength (VS)", tag = "B") +
  theme_classic(base_size = 8) +
  theme(axis.text.y = element_text(size = 6),
        plot.tag = element_text(size = 10.5, face = "bold"))

# --- C : CpG ------------------------------------------------------------------
comp <- read_csv(file.path(dir_outputs, "15_sequence_composition", "tables",
                           "composition_per_sequence.csv"), show_col_types = FALSE) %>%
  filter(is.finite(cpg_obs_exp)) %>%
  mutate(grp = ifelse(TE_ID %in% CAND, "climate_associated", "background")) %>%
  filter(!(grp == "background" & TE_ID %in% mem$TE_ID)) %>%   # drop the other union members
  mutate(grp = factor(grp, levels = c("background","climate_associated")),
         pct = 100 * cpg_obs_exp)
med <- comp %>% group_by(grp) %>% summarise(m = median(pct), .groups = "drop")
w <- wilcox.test(pct ~ grp, data = comp, exact = FALSE)
stars <- ifelse(w$p.value < .001, "***", ifelse(w$p.value < .01, "**",
         ifelse(w$p.value < .05, "*", "n.s.")))
n_lab <- sprintf("Non-climate-associated (n = %d)", sum(comp$grp == "background"))
c_lab <- sprintf("Climate-associated (n = %d)", sum(comp$grp == "climate_associated"))
pC <- ggplot(comp, aes(pct, fill = grp, colour = grp)) +
  geom_density(alpha = .45, linewidth = .7, adjust = 1.1) +
  geom_vline(data = med, aes(xintercept = m, colour = grp), linetype = 2,
             linewidth = .7, show.legend = FALSE) +
  annotate("text", x = mean(med$m), y = Inf, label = stars, vjust = 1.8,
           size = 5, fontface = "bold", colour = "grey20") +
  scale_fill_manual(values = GRP_COL, labels = c(n_lab, c_lab), name = NULL) +
  scale_colour_manual(values = GRP_COL, labels = c(n_lab, c_lab), name = NULL) +
  labs(x = "CpG content (% of the amount expected by chance)", y = "Density", tag = "C",
       caption = sprintf("Wilcoxon p = %.3g   (median %.0f%% vs %.0f%%)", w$p.value,
                         med$m[2], med$m[1])) +
  theme_classic(base_size = 8) +
  theme(legend.position = "top", legend.key.size = unit(.3,"cm"),
        plot.caption = element_text(hjust = 0, size = 7, colour = "grey35"),
        plot.tag = element_text(size = 10.5, face = "bold"))

fig <- (pA | pB) / pC + plot_layout(heights = c(1.25, .8), widths = c(1.15, 1))
for (ext in c("pdf","svg")) ggsave(file.path(dir_fig, paste0("Fig3_", length(CAND), "TE.", ext)),
  fig, width = 180/25.4, height = 7.2, device = if (ext=="pdf") grDevices::cairo_pdf else svglite::svglite)
message(sprintf("Fig3 for %d candidates written. CpG: median %.1f%% vs %.1f%%, p = %.3g",
                length(CAND), med$m[2], med$m[1], w$p.value))
log_session(TAG)
