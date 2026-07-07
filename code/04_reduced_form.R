## =====================================================================
## 04_reduced_form.R -- REDUCED-FORM findings of the US application.
## (1) In-sample conditional-mean decomposition for the 8 focus series:
##     nonlinear adjustment Lambda_mu mu(x_t) (blue) vs linear effect
##     A* x_t (red), 16/84 bands, NBER recessions shaded.
## (2) Share of each variable's variance carried by the (orthogonalised)
##     nonlinear factor block.
##
##   Rscript code/04_reduced_form.R
## Reads cache/{nl_adjustment_summary,lin_effect_summary,varexpl_orth_draws}.rds
## Writes output/figures/{insample_decomp_focus8.pdf, variation_expl.pdf}
## =====================================================================
local({
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  root <- if (length(f)) normalizePath(file.path(dirname(f), "..")) else normalizePath(getwd())
  assign("ROOT", root, envir = .GlobalEnv)
})
source(file.path(ROOT, "code", "config.R"))

## ---- (1) in-sample decomposition ------------------------------------
NL <- readRDS(file.path(CACHE_DIR, "nl_adjustment_summary.rds"))
LN <- readRDS(file.path(CACHE_DIR, "lin_effect_summary.rds"))
dts <- NL$dates
focus <- FOCUS; flab <- FLAB
rec <- list(c("1990-07-01","1991-03-01"), c("2001-03-01","2001-11-01"),
            c("2007-12-01","2009-06-01"))
BL <- "dodgerblue4"; RD <- "red3"

pdf(file.path(FIG_DIR, "insample_decomp_focus8.pdf"), width = 10, height = 8.0)
par(mfrow = c(3,3), mar = c(2.2,3.0,2.0,3.0), oma = c(1.2,1.6,2.4,1.6),
    mgp = c(2,0.5,0), tcl = -0.3, cex.axis = 0.78)
for (v in focus) {
  nm <- NL$NL.med[,v]; nlo <- NL$NL.lo[,v]; nhi <- NL$NL.hi[,v]
  lm <- LN$LIN.med[,v]; llo <- LN$LIN.lo[,v]; lhi <- LN$LIN.hi[,v]
  lyl <- range(c(lm,llo,lhi,0))
  plot(dts, lm, type="n", ylim=lyl, axes=FALSE, xlab="", ylab="", main=flab[v])
  for (r in rec) rect(as.Date(r[1]), lyl[1], as.Date(r[2]), lyl[2], col=grey(0.9), border=NA)
  polygon(c(dts,rev(dts)), c(llo,rev(lhi)), col=adjustcolor("orange",0.40), border=NA)
  lines(dts, lm, col=RD, lwd=1.5); axis(4, col=RD, col.axis=RD, col.ticks=RD)
  axis.Date(1, at = seq(as.Date("1985-01-01"), as.Date("2020-01-01"), by="5 years"), format="%Y"); box()
  par(new=TRUE); nyl <- range(c(nm,nlo,nhi,0))
  plot(dts, nm, type="n", ylim=nyl, axes=FALSE, xlab="", ylab="")
  polygon(c(dts,rev(dts)), c(nlo,rev(nhi)), col=adjustcolor("dodgerblue",0.30), border=NA)
  lines(dts, nm, col=BL, lwd=1.8); abline(h=0, col=BL, lty=3)
  axis(2, col=BL, col.axis=BL, col.ticks=BL)
  legend("topright", legend = sprintf("rel. st. dev. = %.2f",
         sd(nm,na.rm=TRUE)/sd(lm,na.rm=TRUE)), bty="n", cex=0.78, inset=c(0.01,0.01))
}
plot.new()
legend("center", c("Nonlinear adj. (left axis)","Linear effect (right axis)",
                   "16/84 band","NBER recession"),
       col=c(BL,RD,NA,NA), lwd=c(2,2,NA,NA), fill=c(NA,NA,NA,grey(0.9)),
       border=NA, pch=c(NA,NA,15,NA), bty="n", cex=0.95)
mtext(expression("In-sample conditional-mean decomposition:  " * Lambda[mu]*mu*"("*x[t]*")  (blue, left)  and  " * A^"*"*x[t] * "  (red, right)"),
      outer=TRUE, side=3, line=0.6, font=2, cex=0.85)
invisible(dev.off()); cat("wrote insample_decomp_focus8.pdf\n")

## ---- (2) variance share of the nonlinear block ----------------------
V <- readRDS(file.path(CACHE_DIR, "varexpl_orth_draws.rds"))
shares <- V$vexpl; vn <- V$var.names
ord <- order(apply(shares, 2, median), decreasing = TRUE)
shares <- shares[, ord]; vn <- vn[ord]
pdf(file.path(FIG_DIR, "variation_expl.pdf"), width = 7.5, height = 4.2)
par(mar = c(6.5, 4.2, 1.4, 1.2), mgp = c(2.6, 0.6, 0))
boxplot(shares, las = 2, outline = FALSE, names = vn, border = "gray30",
        whisklty = 0, staplelty = 0, ylab = "% variance from nonlinear block",
        col = adjustcolor("forestgreen", 0.25), cex.axis = 0.95, ylim = c(0, max(shares)))
abline(h = mean(apply(shares, 2, mean)), col = "red", lwd = 2)
invisible(dev.off()); cat("wrote variation_expl.pdf\n")
cat("Reduced-form figures written to:", FIG_DIR, "\n")
