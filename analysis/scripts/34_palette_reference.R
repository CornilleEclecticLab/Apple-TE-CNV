# =============================================================================
# Script:  34_palette_reference.R
# Purpose: Reference swatch sheet for the palette used across Figures 1, 2
#          and 3. One colour per category, the same one everywhere. NO GREEN
#          (green merges with red under deuteranopia and protanopia); red is
#          kept. Validated legend by legend with CIEDE2000 on simulated
#          colours (normal + deuteranopia + protanopia + tritanopia).
# Author:  Antho
# Outputs: PALETTE_reference.pdf and .svg
# =============================================================================
OUT <- Sys.getenv("FIG_OUT", unset = file.path("..", "revision_finale", "figures_cvd"))

blocks <- list(
  list(t="TE class   -- Fig. 1E-F   (colours kept as published)       dE 23",
       k=c("Class I        red","Class II       blue","n.a.           mid grey"),
       v=c("#FF0000","#0000FF","#6E6E6E")),
  list(t="TE order / family  -- Fig. 1E-F, Fig. 2E-J, Fig. 3A-B     dE 30",
       k=c("LTR            orange","LINE           sky blue","TIR            indigo",
           "noCat          mid grey","other          black"),
       v=c("#E69F00","#56B4E9","#332288","#6E6E6E","#000000")),
  list(t="TE superfamily     -- Fig. 1E-F, Fig. 3A-B                dE 24",
       k=c("Gypsy          very pale blue","Copia          dark red","n.a.           mid grey"),
       v=c("#CCE6FF","#A50F15","#6E6E6E")),
  list(t="Cell gradient       -- Fig. 1E-F   blue / white / red       dE 53",
       k=c("low  (negative)","middle","high (positive)"), v=c("#0000FF","#FFFFFF","#FF0000")),
  list(t="Genome intersections -- Fig. 2A-D                         dE 24",
       k=c("Domesticated","Wild","Both","Shared by all"),
       v=c("#E69F00","#0072B2","#000000","#994455")),
  list(t="Family composition  -- Fig. 2E-F   (10 categories)        dE 15",
       k=c("LTR_Other      orange","LINE           sky blue","TIR            indigo",
           "noCat          mid grey","LTR_Gypsy      very pale blue","LTR_Copia      dark red",
           "Helitron       strong blue","MITE           near black","SSR            light sand",
           "pararetrovirus light grey"),
       v=c("#E69F00","#56B4E9","#332288","#6E6E6E","#CCE6FF","#A50F15",
           "#0072B2","#1A1A1A","#FFD6A5","#AAAAAA")),
  list(t="Nearby genes, GO    -- Fig. 2M-N   ramp + accent          dE 29",
       k=c("No genes at 2,000 bp","GO unavailable","No response to stimulus","Response to stimulus"),
       v=c("#F7F7F7","#888888","#3A3A3A","#E69F00")),
  list(t="Fig. 2G-L  violins / boxes, single fill",
       k=c("fill","stroke"), v=c("#BBBBBB","#000000")),
  list(t="Fig. 3C             -- the two TE sets                    dE 72",
       k=c("34 candidates","523 background"), v=c("#002255","#E69F00")),
  list(t="Fig. 1E-F only: \"other\" pools SINE, SSR, Helitron,",
       k=character(0), v=character(0)),
  list(t="pararetrovirus, MITE and mixture -- 12 TEs of 2420 (0.5%), 2 of 1360,",
       k=character(0), v=character(0)),
  list(t="i.e. one-pixel columns. In Fig. 2E-F they are named individually.",
       k=character(0), v=character(0))
)

draw <- function() {
  par(mar=c(0,0,0,0)); plot.new()
  n <- sum(sapply(blocks, function(b) length(b$k))) + 2*length(blocks)
  plot.window(c(0,1), c(n,0))
  y <- 0.4
  for (b in blocks) {
    text(0.02, y, b$t, adj=0, cex=0.72, font=2); y <- y + 1.1
    for (i in seq_along(b$k)) {
      rect(0.05, y-0.34, 0.14, y+0.34, col=b$v[i], border="grey40", lwd=0.4)
      text(0.17, y, b$k[i], adj=0, cex=0.72)
      text(0.66, y, b$v[i], adj=0, cex=0.72, family="mono")
      y <- y + 0.85
    }
    y <- y + 0.7
  }
}
h <- 15.5
cairo_pdf(file.path(OUT,"PALETTE_reference.pdf"), width=6.2, height=h); draw(); dev.off()
svg(file.path(OUT,"PALETTE_reference.svg"), width=6.2, height=h); draw(); dev.off()
cat("written: PALETTE_reference.pdf / .svg\n")
