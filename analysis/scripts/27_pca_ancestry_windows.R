# =============================================================================
# Script:  27_pca_ancestry_windows.R
# Purpose: The PCA asked for on R1.7: instead of using all 1,297 neutral SNPs,
#          restrict them to the genomic windows that RFMix assigns to one
#          ancestry, and see where the domesticated accessions fall.
#
#          If the domesticated genome is a mosaic rather than a uniform blend,
#          the same accessions should sit next to M. sieversii when only
#          sieversii-ancestry windows are used, and next to M. sylvestris when
#          only sylvestris-ancestry windows are used. A uniform blend would put
#          them in the same intermediate position in both.
# Author:  Antho
# =============================================================================
source("scripts/00_setup.R")
suppressMessages({library(dplyr); library(readr); library(tibble); library(ggplot2); library(patchwork)})
TAG <- "27_pca_ancestry_windows"
dir_out <- file.path(dir_outputs, TAG); dir_tab <- file.path(dir_out,"tables"); dir_fig <- file.path(dir_out,"figures")
for (d in c(dir_out,dir_tab,dir_fig)) dir.create(d, showWarnings=FALSE, recursive=TRUE)
SRC <- "~/Documents/TE_Project/1-PCA"

e <- new.env(); load(file.path(SRC,"dosage_data.rda"), envir=e); dos <- get("dosage_data", e)
meta <- read_csv(file.path(dir_data,"250923_sample_list_TE.csv"), show_col_types=FALSE)
names(meta)[1] <- "Sample_ID"
meta <- meta %>% filter(Sample_ID %in% rownames(dos))
dos  <- dos[meta$Sample_ID, , drop=FALSE]
message(sprintf("%d accessions, %d SNPs", nrow(dos), ncol(dos)))

sp  <- do.call(rbind, strsplit(colnames(dos), "_"))
snp <- tibble(id=colnames(dos), chr=sp[,1], pos=as.numeric(sp[,2]))
rf  <- read.delim(file.path(SRC,"RFMix_out.whole_genome_plot.tsv")); names(rf)[2] <- "start"
snp$frac_sie <- vapply(seq_len(nrow(snp)), function(i) {
  b <- rf[rf$chr==snp$chr[i] & rf$start<=snp$pos[i] & rf$end>=snp$pos[i], ]
  if (nrow(b)==0) NA_real_ else mean(b$pop=="Sie")
}, numeric(1))
write_csv(snp, file.path(dir_tab,"snp_ancestry_fraction.csv"))

sets <- list(
  "All neutral SNPs"                = snp$id,
  "sieversii-ancestry windows"      = snp$id[which(snp$frac_sie > .70)],
  "sylvestris-ancestry windows"     = snp$id[which(snp$frac_sie < .30)])
message(paste(sprintf("  %s : %d SNPs", names(sets), lengths(sets)), collapse="\n"))

grp <- meta$Code_group_geo_residual
sp_of <- function(g) case_when(grepl("^DOM", g) ~ "M. domestica", grepl("^SIE", g) ~ "M. sieversii",
                               grepl("^SYL", g) ~ "M. sylvestris", grepl("^OR",  g) ~ "M. orientalis",
                               grepl("^BAC", g) ~ "M. baccata", TRUE ~ "other")
meta$species <- sp_of(grp)
COL <- c("M. domestica"="#C4571F","M. sieversii"="#7B3294","M. sylvestris"="#2E6DA4",
         "M. orientalis"="#B8A317","M. baccata"="#4DAF4A","other"="grey70")

run_pca <- function(ids, lab) {
  X <- dos[, ids, drop=FALSE]
  X <- X[, apply(X, 2, function(c) var(c, na.rm=TRUE) > 0), drop=FALSE]
  X <- apply(X, 2, function(c) { c[is.na(c)] <- mean(c, na.rm=TRUE); c })
  p <- prcomp(X, scale.=TRUE)
  ve <- 100*p$sdev^2/sum(p$sdev^2)
  tibble(Sample_ID=meta$Sample_ID, species=meta$species,
         PC1=p$x[,1], PC2=p$x[,2], set=lab,
         lab=sprintf("%s  (n = %d)\nPC1 %.1f%%, PC2 %.1f%%", lab, ncol(X), ve[1], ve[2]))
}
pca <- bind_rows(Map(run_pca, sets, names(sets)))
write_csv(pca, file.path(dir_tab,"pca_coordinates.csv"))

# Where do the domesticated sit, between the two wild progenitors, on PC1?
pos <- pca %>% group_by(set) %>% group_modify(~{
  cen <- .x %>% group_by(species) %>% summarise(m=mean(PC1), .groups="drop")
  g <- function(s) cen$m[cen$species==s]
  tibble(sieversii=g("M. sieversii"), sylvestris=g("M. sylvestris"),
         domestica=g("M. domestica"),
         # 0 = on top of sylvestris, 1 = on top of sieversii
         relative_position=(g("M. domestica")-g("M. sylvestris"))/(g("M. sieversii")-g("M. sylvestris")))
}) %>% ungroup()
write_csv(pos, file.path(dir_tab,"domesticated_position.csv"))
message("\nPosition of the domesticated centroid on PC1 (0 = sylvestris, 1 = sieversii):")
print(as.data.frame(pos %>% select(set, relative_position)))

pca$lab <- factor(pca$lab, levels=unique(pca$lab[order(match(pca$set, names(sets)))]))
p <- ggplot(pca, aes(PC1, PC2, fill=species)) +
  geom_point(shape=21, size=1.9, colour="white", stroke=.25, alpha=.9) +
  facet_wrap(~lab, nrow=1, scales="free") +
  scale_fill_manual(values=COL, name=NULL) +
  labs(x="PC1", y="PC2") +
  theme_bw(base_size=8) +
  theme(legend.position="bottom", legend.key.size=unit(.3,"cm"),
        strip.background=element_blank(), strip.text=element_text(size=7, face="bold"),
        panel.grid.minor=element_blank())
for (ext in c("pdf","svg"))
  ggsave(file.path(dir_fig, paste0("FigS23_PCA_ancestry_windows.",ext)), p,
         width=180/25.4, height=3.3,
         device=if(ext=="pdf") grDevices::cairo_pdf else svglite::svglite)
log_session(TAG)
