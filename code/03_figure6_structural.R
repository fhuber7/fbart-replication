## =====================================================================
## 03_figure6_structural.R -- MAIN STRUCTURAL RESULT (paper Figure 6).
## GIRFs of the 8 focus variables to large/small financial and monetary
## shocks: positive (blue) vs sign-flipped negative (red), 16/84 bands,
## with the tightening-vs-easing asymmetry probability
##   P = Pr(|positive response| > |negative response|)
## printed above each panel at selected horizons.
##
##   Rscript code/03_figure6_structural.R
## Reads cache/{irf_quantiles_allsizes,irf_asymmetry_measures}.rds
## Writes output/figures/irf_focus8_*_2x4.{pdf,png} + irf_focus8_legend.pdf
## =====================================================================
local({
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  root <- if (length(f)) normalizePath(file.path(dirname(f), "..")) else normalizePath(getwd())
  assign("ROOT", root, envir = .GlobalEnv)
})
source(file.path(ROOT, "code", "config.R"))
suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) })

S <- readRDS(file.path(CACHE_DIR, "irf_quantiles_allsizes.rds"))
A <- readRDS(file.path(CACHE_DIR, "irf_asymmetry_measures.rds"))
Q <- S$Q; nhor <- S$nhor; BAND <- c("q16","q84")
xs <- seq_len(nhor)
PHOR  <- c(2, 5, 8, 11, 14, 17, 20, 23)  # horizons annotated with P
PGRN  <- "black"      # P numbers
PGRNL <- "#bdbdbd"    # light grey dashed connectors
BL <- "#1f4e79"; RD <- "#9e1b1b"
focus <- FOCUS; flab <- FLAB

base_thm <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "none",
        strip.text = element_text(face = "bold", size = 10))

mk24 <- function(shock, sz){
  sp <- as.character(sz); sn <- as.character(-sz)
  pdat <- bind_rows(lapply(focus, function(v){
    data.frame(series = v, label = unname(flab[v]), horizon = xs,
               pm =  Q["q50",v,,sp,shock], plo =  Q[BAND[1],v,,sp,shock], phi =  Q[BAND[2],v,,sp,shock],
               nm = -Q["q50",v,,sn,shock], nlo = -Q[BAND[2],v,,sn,shock], nhi = -Q[BAND[1],v,,sn,shock])
  }))
  pdat$label <- factor(pdat$label, levels = unname(flab[focus]))

  pair <- if (abs(sz) >= 10) "large" else "small"
  Pdat <- bind_rows(lapply(focus, function(v){
    phi <-  Q[BAND[2],v,,sp,shock]; nhi <- -Q[BAND[1],v,,sn,shock]
    plo <-  Q[BAND[1],v,,sp,shock]; nlo <- -Q[BAND[2],v,,sn,shock]
    ub  <- pmax(phi, nhi)
    dmax <- max(c(phi,nhi,0)); dmin <- min(c(plo,nlo,0)); rng <- dmax - dmin
    data.frame(label = unname(flab[v]), horizon = PHOR,
               Plab = sub("^0","", formatC(A$P[v,PHOR,pair,shock], format="f", digits=2)),
               ynum = dmax + 0.11*rng, yend = ub[PHOR])
  }))
  Pdat$label <- factor(Pdat$label, levels = unname(flab[focus]))

  p <- ggplot(pdat, aes(horizon)) +
    geom_hline(yintercept = 0, colour = "grey55", linetype = 3) +
    geom_ribbon(aes(ymin = plo, ymax = phi), fill = BL, alpha = .25) +
    geom_ribbon(aes(ymin = nlo, ymax = nhi), fill = RD, alpha = .25) +
    geom_line(aes(y = pm), colour = BL, linewidth = 1.0) +
    geom_line(aes(y = nm), colour = RD, linewidth = 1.0) +
    geom_segment(data = Pdat, aes(x = horizon, xend = horizon, y = yend, yend = ynum),
                 colour = PGRNL, linewidth = 0.25, linetype = "22", inherit.aes = FALSE) +
    geom_text(data = Pdat, aes(x = horizon, y = ynum, label = Plab),
              vjust = -0.15, size = 3.1, colour = PGRN, inherit.aes = FALSE) +
    facet_wrap(~ label, nrow = 2, ncol = 4, scales = "free_y") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    labs(x = "Horizon (quarters)", y = NULL) + base_thm
  fn <- file.path(FIG_DIR, sprintf("irf_focus8_%s_size_%d_2x4.pdf", shock, sz))
  ggsave(fn, p, width = 12, height = 6.0)
  ggsave(sub("\\.pdf$", ".png", fn), p, width = 12, height = 6.0, dpi = 150)
  cat("wrote", basename(fn), "\n")
}
for (shk in c("Financial","Monetary")) for (sz in c(5, 14)) mk24(shk, sz)

## shared horizontal legend for the two Figure-6 panels
pdf(file.path(FIG_DIR, "irf_focus8_legend.pdf"), width = 6.6, height = 0.34)
par(mar = c(0,0,0,0)); plot.new()
legend("center", legend = c("Positive", "Negative (mirrored)"),
       col = c(BL, RD), lwd = 3, horiz = TRUE, bty = "n", seg.len = 2.4, cex = 1.05)
invisible(dev.off()); cat("wrote irf_focus8_legend.pdf\n")
cat("Figure 6 written to:", FIG_DIR, "\n")
