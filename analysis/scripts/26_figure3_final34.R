# =============================================================================
# Script:  26_figure3_final34.R
# Purpose: The published Figure 3, built on the 34 replicated candidates.
#          Panel C has two halves, as the legend states: the CpG distribution
#          on the left, and the two deamination products on the right. The
#          quick variant script only drew the left half.
# Author:  Antho
# =============================================================================
source("scripts/00_setup.R")
suppressMessages({library(dplyr); library(readr); library(tibble); library(ggplot2)
                  library(patchwork); library(sf); library(Biostrings)})
TAG <- "26_figure3_final34"
dir_fig <- file.path(dir_outputs, TAG, "figures")
dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)

mem  <- read_csv(file.path(dir_outputs,"18_climate_pc_diagnostics","tables",
                           "candidate_membership.csv"), show_col_types=FALSE)
CAND <- mem$TE_ID[mem$in_modelA]
SF_COL <- c(LTR_Copia="lightgreen", LTR_Gypsy="lightblue", LTR="orange", noCat="grey", TIR="purple")
SF_LAB <- c(LTR_Copia="Copia", LTR_Gypsy="Gypsy", LTR="LTR", noCat="noCat", TIR="TIR")
GRP <- c(background="#4C72B0", climate="#C44E52")
# Sampling-group palette, as in 07_map_figure3A.R. The three groups are the
# subsets the candidate rule replicates across, so they belong on the map.
GRP_COL <- c(SYLW="#2E6DA4", SYLE="#B8A317", SYL_Admix="#C4571F")
GROUP_OF <- function(x) dplyr::case_when(
  x == "Q1" ~ "SYLW", x == "Q5" ~ "SYLE",
  x %in% c("Q1_Admix","Q5_Admix") ~ "SYL_Admix", TRUE ~ NA_character_)
stars <- function(p) ifelse(p<.001,"***", ifelse(p<.01,"**", ifelse(p<.05,"*","n.s.")))

# --- A : map -----------------------------------------------------------------
arr <- read_csv(file.path(dir_outputs,"07_map_figure3A","tables","arrow_geometry.csv"),
                show_col_types=FALSE) %>% filter(TE_ID %in% CAND)
sites <- read_csv(file.path(dir_data,"250923_sample_list_TE.csv"), show_col_types=FALSE) %>%
  filter(Species=="M. sylvestris", !is.na(latitude), !is.na(longitude)) %>%
  mutate(group = GROUP_OF(group_structure_CORE_admix_K6)) %>%
  filter(!is.na(group))
sxy <- st_coordinates(st_transform(st_as_sf(sites, coords=c("longitude","latitude"), crs=4326), 3035))
sites$x <- sxy[,1]; sites$y <- sxy[,2]
world <- rnaturalearth::ne_countries(scale="medium", returnclass="sf") %>% st_transform(3035)
pad <- 210e3
xlim <- c(min(sites$x)-pad, max(sites$x)+pad); ylim <- c(min(sites$y)-pad, max(sites$y)+pad)
pA <- ggplot() +
  geom_sf(data=world, fill="#F4F1EA", colour="#BFB8AC", linewidth=.25) +
  geom_point(data=sites, aes(x,y,fill=group), shape=21, size=1.5,
             colour="white", stroke=.25, alpha=.95) +
  geom_segment(data=arr, aes(x,y,xend=xend,yend=yend,colour=superfamily,linewidth=true_vs_corrected),
    arrow=grid::arrow(length=unit(.14,"cm"),type="closed"), lineend="round", alpha=.9) +
  geom_point(data=arr[1,], aes(x,y), shape=21, size=2.2, fill="white", colour="grey15", stroke=.6) +
  scale_colour_manual(values=SF_COL, labels=SF_LAB, name="TE superfamily") +
  scale_fill_manual(values=GRP_COL, name="Sampling group") +
  scale_linewidth_continuous(range=c(.4,1.5), guide="none") +
  guides(colour=guide_legend(order=1, nrow=2,
           override.aes=list(linewidth=1.2, alpha=1)),
         fill=guide_legend(order=2, nrow=1, override.aes=list(size=2.2))) +
  coord_sf(xlim=xlim, ylim=ylim, expand=FALSE, crs=3035) +
  labs(tag="A") + theme_void(base_size=8) +
  theme(panel.background=element_rect(fill="#DCE9F2", colour="grey40"),
        legend.position="bottom", legend.box="vertical", legend.box.just="left",
        legend.spacing.y=unit(.02,"cm"), legend.key.size=unit(.26,"cm"),
        legend.text=element_text(size=6.2), legend.title=element_text(size=6.6),
        plot.tag=element_text(size=10.5, face="bold"), plot.margin=margin(1,2,1,2))

# --- B : vector strength ------------------------------------------------------
vs <- read_csv(file.path(dir_outputs,"03_spatial_vs","tables","VS_table_38_candidates.csv"),
               show_col_types=FALSE) %>% filter(TE_ID %in% CAND)
pB <- ggplot(vs, aes(reorder(TE_ID,true_vs_corrected), true_vs_corrected, fill=superfamily)) +
  geom_col(width=.78) + geom_hline(yintercept=.5, linetype=2, colour="grey40") +
  coord_flip() + scale_fill_manual(values=SF_COL, guide="none") +
  labs(x=NULL, y="Spatial vector strength (VS)", tag="B") +
  theme_classic(base_size=8) +
  theme(axis.text.y=element_text(size=4.4), axis.ticks.y=element_blank(),
        plot.tag=element_text(size=10.5, face="bold"), plot.margin=margin(1,2,1,2))

# --- C : CpG and its two deamination products --------------------------------
db <- readDNAStringSet(file.path(dir_data,"TE_consensus_db.fasta"))
names(db) <- sub(" .*$","",names(db))
base <- read_csv(file.path(dir_outputs,"15_sequence_composition","tables",
                           "composition_per_sequence.csv"), show_col_types=FALSE)
ids <- base$TE_ID[base$grp=="background" | base$TE_ID %in% CAND]
s <- db[intersect(ids, names(db))]
mono <- letterFrequency(s, c("A","C","G","T")); di <- dinucleotideFrequency(s); L <- width(s)
oe <- function(d) di[,d]*L/(mono[,substr(d,1,1)]*mono[,substr(d,2,2)])
comp <- tibble(TE_ID=names(s), CG=oe("CG"), TG=oe("TG"), CA=oe("CA")) %>%
  mutate(grp=factor(ifelse(TE_ID %in% CAND,"climate","background"),
                    levels=c("background","climate"))) %>%
  filter(is.finite(CG), is.finite(TG), is.finite(CA))
nb <- sum(comp$grp=="background"); nc <- sum(comp$grp=="climate")
labs2 <- c(sprintf("Non-climate-associated (n = %d)", nb), sprintf("Climate-associated (n = %d)", nc))

cg <- comp %>% mutate(pct=100*CG)
medc <- cg %>% group_by(grp) %>% summarise(m=median(pct), .groups="drop")
w_cg <- wilcox.test(pct ~ grp, cg, exact=FALSE)
dmax <- max(vapply(split(cg$pct, cg$grp), function(x) max(density(x, adjust=1.1)$y), numeric(1)))
bx <- sort(medc$m); by <- dmax*1.10; bt <- dmax*.035
pC1 <- ggplot(cg, aes(pct, fill=grp, colour=grp)) +
  geom_density(alpha=.45, linewidth=.7, adjust=1.1) +
  geom_vline(data=medc, aes(xintercept=m, colour=grp), linetype=2, linewidth=.7, show.legend=FALSE) +
  annotate("segment", x=bx[1], xend=bx[2], y=by, yend=by, linewidth=.5, colour="grey20") +
  annotate("segment", x=bx[1], xend=bx[1], y=by, yend=by-bt, linewidth=.5, colour="grey20") +
  annotate("segment", x=bx[2], xend=bx[2], y=by, yend=by-bt, linewidth=.5, colour="grey20") +
  annotate("text", x=mean(bx), y=by+bt*.35, label=stars(w_cg$p.value), size=4.6,
           vjust=0, fontface="bold", colour="grey20") +
  expand_limits(y=dmax*1.24) +
  scale_fill_manual(values=GRP, labels=labs2, name=NULL) +
  scale_colour_manual(values=GRP, labels=labs2, name=NULL) +
  labs(x="CpG content (% of the amount expected by chance)", y="Density", tag="C",
       caption="Wilcoxon rank-sum test:  *** p < 0.001   ** p < 0.01   * p < 0.05   n.s. not significant") +
  theme_classic(base_size=8) +
  theme(legend.position="top", legend.key.size=unit(.3,"cm"), legend.text=element_text(size=6.4),
        plot.caption=element_text(hjust=0, size=6.4, colour="grey35"),
        plot.tag=element_text(size=10.5, face="bold"), plot.margin=margin(1,2,1,2))

trio <- bind_rows(
  comp %>% transmute(grp, pct=100*TG, dinuc="TpG  (deamination product)"),
  comp %>% transmute(grp, pct=100*CA, dinuc="CpA  (deamination product)")) %>%
  mutate(dinuc=factor(dinuc, levels=c("TpG  (deamination product)","CpA  (deamination product)")))
medt <- trio %>% group_by(dinuc, grp) %>% summarise(m=median(pct), .groups="drop")
ann <- trio %>% group_by(dinuc) %>%
  summarise(p=wilcox.test(pct ~ grp, exact=FALSE)$p.value,
            mc=median(pct[grp=="climate"]), mb=median(pct[grp=="background"]), .groups="drop") %>%
  mutate(lab=sprintf("%.0f%% vs %.0f%%   %s", mc, mb, stars(p)))
pC2 <- ggplot(trio, aes(pct, fill=grp, colour=grp)) +
  geom_density(alpha=.45, linewidth=.6, adjust=1.1) +
  geom_vline(data=medt, aes(xintercept=m, colour=grp), linetype=2, linewidth=.6, show.legend=FALSE) +
  geom_vline(xintercept=100, colour="grey55", linewidth=.4) +
  geom_text(data=ann, aes(x=Inf, y=Inf, label=lab), inherit.aes=FALSE,
            hjust=1.05, vjust=1.6, size=2.6, colour="grey25") +
  facet_wrap(~dinuc, ncol=1, scales="free_y") +
  scale_fill_manual(values=GRP, guide="none") + scale_colour_manual(values=GRP, guide="none") +
  # No caption here: side by side with the left panel the two captions collided.
  # The manuscript legend states that CpG itself is shown on the left.
  labs(x="Content (% of expected; grey line = expected by chance)", y="Density") +
  theme_classic(base_size=8) +
  theme(strip.background=element_blank(), strip.text=element_text(face="bold", hjust=0, size=7),
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        plot.margin=margin(1,2,1,2))

fig <- (pA | pB) / (pC1 | pC2) + plot_layout(heights=c(1.15,.85), widths=c(1.15,1))
for (e in c("pdf","svg"))
  ggsave(file.path(dir_fig, paste0("FIG3.",e)), fig, width=180/25.4, height=7.6,
         device=if(e=="pdf") grDevices::cairo_pdf else svglite::svglite)
message(sprintf("FIG3 (%d candidates) | CpG p = %.3g | TpG p = %.3g | CpA p = %.3g",
                length(CAND), w_cg$p.value, ann$p[1], ann$p[2]))
log_session(TAG)
