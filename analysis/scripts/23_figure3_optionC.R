# Figure 3 as it would look under option C (82 TEs, climate + lat/long, FDR).
source("scripts/00_setup.R")
suppressMessages({library(dplyr); library(readr); library(ggplot2); library(patchwork); library(sf)})
TAG <- "23_figure3_optionC"; dir_fig <- file.path(dir_outputs, TAG, "figures")
dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)
SF_COL <- c(LTR_Copia="lightgreen", LTR_Gypsy="lightblue", LTR="orange", noCat="grey", TIR="purple")
SF_LAB <- c(LTR_Copia="Copia", LTR_Gypsy="Gypsy", LTR="LTR", noCat="noCat", TIR="TIR")
GRP <- c(background="#4C72B0", climate="#C44E52")
vs <- read_csv(file.path(dir_outputs,"22_optionC_check","tables","VS_optionC.csv"), show_col_types=FALSE)
comp <- read_csv(file.path(dir_outputs,"22_optionC_check","tables","composition_optionC.csv"), show_col_types=FALSE) %>%
  mutate(status=factor(status, levels=c("background","climate")), pct=100*cpg_oe)

# A : map, arrows rebuilt from the fitted angle and VS
sites <- read_csv(file.path(dir_data,"250923_sample_list_TE.csv"), show_col_types=FALSE) %>%
  filter(Species=="M. sylvestris", !is.na(latitude), !is.na(longitude))
o_lon <- mean(sites$longitude); o_lat <- mean(sites$latitude)
KM <- 900
arr <- vs %>% mutate(len=true_vs*KM*1000)
pr <- function(lon,lat) st_coordinates(st_transform(st_as_sf(data.frame(lon,lat),
        coords=c("lon","lat"), crs=4326), 3035))
o <- pr(o_lon,o_lat); arr$x <- o[1,1]; arr$y <- o[1,2]
arr$xend <- arr$x + cos(arr$angle_rad)*arr$len; arr$yend <- arr$y + sin(arr$angle_rad)*arr$len
sxy <- pr(sites$longitude, sites$latitude); sites$x <- sxy[,1]; sites$y <- sxy[,2]
world <- rnaturalearth::ne_countries(scale="medium", returnclass="sf") %>% st_transform(3035)
pad <- 210e3
pA <- ggplot() +
  geom_sf(data=world, fill="#F4F1EA", colour="#BFB8AC", linewidth=.25) +
  geom_point(data=sites, aes(x,y), colour="grey45", size=1, alpha=.7) +
  geom_segment(data=arr, aes(x,y,xend=xend,yend=yend,colour=superfamily,linewidth=true_vs),
    arrow=grid::arrow(length=unit(.12,"cm"),type="closed"), lineend="round", alpha=.8) +
  scale_colour_manual(values=SF_COL, labels=SF_LAB, name="TE superfamily") +
  scale_linewidth_continuous(range=c(.3,1.3), guide="none") +
  coord_sf(xlim=c(min(sites$x)-pad,max(sites$x)+pad),
           ylim=c(min(sites$y)-pad,max(sites$y)+pad), expand=FALSE, crs=3035) +
  labs(tag="A") + theme_void(base_size=8) +
  theme(panel.background=element_rect(fill="#DCE9F2", colour="grey40"),
        legend.position="bottom", legend.key.size=unit(.26,"cm"),
        plot.tag=element_text(size=10.5, face="bold"))

pB <- ggplot(vs, aes(reorder(TE_ID,true_vs), true_vs, fill=superfamily)) +
  geom_col(width=.8) + geom_hline(yintercept=.5, linetype=2, colour="grey40") +
  coord_flip() + scale_fill_manual(values=SF_COL, guide="none") +
  labs(x=NULL, y="Spatial vector strength (VS)", tag="B") +
  theme_classic(base_size=8) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        plot.tag=element_text(size=10.5, face="bold"))

med <- comp %>% group_by(status) %>% summarise(m=median(pct), .groups="drop")
w <- wilcox.test(pct ~ status, data=comp, exact=FALSE)
pC <- ggplot(comp, aes(pct, fill=status, colour=status)) +
  geom_density(alpha=.45, linewidth=.7, adjust=1.1) +
  geom_vline(data=med, aes(xintercept=m, colour=status), linetype=2, linewidth=.7, show.legend=FALSE) +
  annotate("text", x=mean(med$m), y=Inf, label="n.s.", vjust=1.8, size=4.2, fontface="bold", colour="grey20") +
  scale_fill_manual(values=GRP, labels=c(sprintf("Non-climate-associated (n = %d)", sum(comp$status=="background")),
                                          sprintf("Climate-associated (n = %d)", sum(comp$status=="climate"))), name=NULL) +
  scale_colour_manual(values=GRP, labels=c(sprintf("Non-climate-associated (n = %d)", sum(comp$status=="background")),
                                            sprintf("Climate-associated (n = %d)", sum(comp$status=="climate"))), name=NULL) +
  labs(x="CpG content (% of the amount expected by chance)", y="Density", tag="C",
       caption=sprintf("Wilcoxon p = %.3g   (median %.0f%% vs %.0f%%)   effect reversed", w$p.value, med$m[2], med$m[1])) +
  theme_classic(base_size=8) +
  theme(legend.position="top", legend.key.size=unit(.3,"cm"),
        plot.caption=element_text(hjust=0, size=7, colour="grey35"),
        plot.tag=element_text(size=10.5, face="bold"))

fig <- (pA | pB) / pC + plot_layout(heights=c(1.25,.8), widths=c(1.15,1))
for (e in c("pdf","svg")) ggsave(file.path(dir_fig, paste0("Fig3_82TE.",e)), fig,
  width=180/25.4, height=7.2, device=if(e=="pdf") grDevices::cairo_pdf else svglite::svglite)
message(sprintf("VS > 0.5 : %d / %d | CpG p = %.3g", sum(vs$true_vs>.5), nrow(vs), w$p.value))
