## =====================================================================
## 02_build_caches.R -- pool the posterior-draw files from 01_estimate.R
## and build the small summary caches consumed by the figure scripts
## (03_figure6_structural.R, 04_reduced_form.R).
##
##   Rscript code/02_build_caches.R
##
## Rebuilds, in cache/: irf_quantiles_allsizes.rds, irf_asymmetry_measures.rds,
## nl_adjustment_summary.rds, lin_effect_summary.rds, varexpl_orth_draws.rds.
## (The package already ships these caches, so 03/04 can run without 01-02.)
## =====================================================================
local({
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  root <- if (length(f)) normalizePath(file.path(dirname(f), "..")) else normalizePath(getwd())
  assign("ROOT", root, envir = .GlobalEnv)
})
source(file.path(ROOT, "code", "config.R"))
suppressMessages(library(abind))

reps <- sort(list.files(DRAWS_DIR, pattern = "^US_baseline_rep[0-9]+\\.rds$", full.names = TRUE))
if (!length(reps)) stop("No draw files in ", DRAWS_DIR, " -- run 01_estimate.R first.")
message("Pooling ", length(reps), " replicate file(s) ...")

IRF <- A <- Lm <- mu <- vexpl <- NULL
X <- vn <- NULL
for (f in reps) {
  r <- readRDS(f)
  IRF   <- if (is.null(IRF))   r$IRF      else abind(IRF,   r$IRF,      along = 1)
  A     <- if (is.null(A))     r$A        else abind(A,     r$A,        along = 1)
  Lm    <- if (is.null(Lm))    r$Lambdamu else abind(Lm,    r$Lambdamu, along = 1)
  mu    <- if (is.null(mu))    r$mu       else abind(mu,    r$mu,       along = 1)
  vexpl <- if (is.null(vexpl)) r$var.expl else abind(vexpl, r$var.expl, along = 1)
  if (is.null(X)) { X <- r$X; vn <- r$var.names }
  rm(r); gc()
}
ndraw <- dim(IRF)[1]; nhor <- dim(IRF)[3]
size.grid <- as.numeric(dimnames(IRF)[[4]]); if (any(is.na(size.grid))) size.grid <- c(-25,-14,-5,5,14,25)
shock.labels <- c("Demand","Monetary","Supply","Financial")
message(sprintf("pooled draws = %d | vars = %d | horizons = %d", ndraw, length(vn), nhor))

## ---- (1) IRF quantile cache -----------------------------------------
qs <- c(0.16, 0.25, 0.50, 0.75, 0.84)
Q  <- apply(IRF, c(2,3,4,5), quantile, qs, na.rm = TRUE)
dimnames(Q) <- list(c("q16","q25","q50","q75","q84"), vn, NULL,
                    as.character(size.grid), shock.labels)
saveRDS(list(Q = Q, var.names = vn, size.grid = size.grid,
             shock.labels = shock.labels, nhor = nhor),
        file.path(CACHE_DIR, "irf_quantiles_allsizes.rds"))
message("  wrote irf_quantiles_allsizes.rds")

## ---- (2) asymmetry measures  P = Pr(|tightening| > |easing|) ---------
hlab   <- dimnames(IRF)[[3]]
shocks <- c(Monetary = "2", Financial = "4")
pairs  <- list(small = c(pos = "5",  neg = "-5"),
               large = c(pos = "14", neg = "-14"))
dn <- list(var = vn, hor = hlab, pair = names(pairs), shock = names(shocks))
mk <- function() array(NA_real_, dim = lengths(dn), dimnames = dn)
P<-mk(); Dmean<-mk(); Dmed<-mk(); Dq05<-mk(); Dq16<-mk(); Dq84<-mk(); Dq95<-mk()
Tmed<-mk(); Tq16<-mk(); Tq84<-mk(); Emed<-mk(); Eq16<-mk(); Eq84<-mk()
for (sh in names(shocks)) { ki <- shocks[sh]
  for (pn in names(pairs)) { pp <- pairs[[pn]]
    for (v in vn) {
      Tdr <- IRF[, v, , pp["pos"], ki]; Edr <- IRF[, v, , pp["neg"], ki]
      absT <- abs(Tdr); absE <- abs(Edr); D <- absT - absE
      ## tie-split so structural mirror at impact maps to the symmetry null 0.5
      P   [v,,pn,sh] <- colMeans(absT > absE) + 0.5 * colMeans(absT == absE)
      Dmean[v,,pn,sh] <- colMeans(D)
      qd <- apply(D, 2, quantile, c(.05,.16,.50,.84,.95), na.rm = TRUE)
      Dq05[v,,pn,sh]<-qd[1,]; Dq16[v,,pn,sh]<-qd[2,]; Dmed[v,,pn,sh]<-qd[3,]
      Dq84[v,,pn,sh]<-qd[4,]; Dq95[v,,pn,sh]<-qd[5,]
      Tq<-apply(Tdr,2,quantile,c(.16,.50,.84),na.rm=TRUE)
      Eq<-apply(Edr,2,quantile,c(.16,.50,.84),na.rm=TRUE)
      Tq16[v,,pn,sh]<-Tq[1,]; Tmed[v,,pn,sh]<-Tq[2,]; Tq84[v,,pn,sh]<-Tq[3,]
      Eq16[v,,pn,sh]<-Eq[1,]; Emed[v,,pn,sh]<-Eq[2,]; Eq84[v,,pn,sh]<-Eq[3,]
    }
  }
}
saveRDS(list(P=P,Dmean=Dmean,Dmed=Dmed,Dq05=Dq05,Dq16=Dq16,Dq84=Dq84,Dq95=Dq95,
             Tmed=Tmed,Tq16=Tq16,Tq84=Tq84,Emed=Emed,Eq16=Eq16,Eq84=Eq84,
             ndraw=ndraw,nhor=nhor,hlab=hlab,shocks=shocks,pairs=pairs,var.names=vn),
        file.path(CACHE_DIR, "irf_asymmetry_measures.rds"))
message("  wrote irf_asymmetry_measures.rds")

## ---- (3) in-sample conditional-mean decomposition -------------------
T_ <- dim(mu)[2]
NL <- array(NA_real_, c(ndraw, T_, length(vn)))
for (d in seq_len(ndraw)) NL[d,,] <- mu[d,,] %*% t(Lm[d,,])
NLm<-apply(NL,c(2,3),median); NLlo<-apply(NL,c(2,3),quantile,0.16,na.rm=TRUE)
NLhi<-apply(NL,c(2,3),quantile,0.84,na.rm=TRUE)
dimnames(NLm)<-dimnames(NLlo)<-dimnames(NLhi)<-list(NULL,vn)
dates <- seq(as.Date("1985-07-01"), by = "quarter", length.out = T_)
saveRDS(list(NL.med=NLm,NL.lo=NLlo,NL.hi=NLhi,dates=dates,var.names=vn),
        file.path(CACHE_DIR, "nl_adjustment_summary.rds"))

LIN <- array(NA_real_, c(ndraw, nrow(X), length(vn)))
for (d in seq_len(ndraw)) LIN[d,,] <- X %*% A[d,,]
LMm<-apply(LIN,c(2,3),median); LMlo<-apply(LIN,c(2,3),quantile,0.16,na.rm=TRUE)
LMhi<-apply(LIN,c(2,3),quantile,0.84,na.rm=TRUE)
dimnames(LMm)<-dimnames(LMlo)<-dimnames(LMhi)<-list(NULL,vn)
saveRDS(list(LIN.med=LMm,LIN.lo=LMlo,LIN.hi=LMhi,var.names=vn),
        file.path(CACHE_DIR, "lin_effect_summary.rds"))
message("  wrote nl_adjustment_summary.rds, lin_effect_summary.rds")

## ---- (4) variance share of the nonlinear block ---------------------
colnames(vexpl) <- vn
if (max(vexpl, na.rm = TRUE) <= 1.5) vexpl <- 100 * vexpl   # -> percent
saveRDS(list(vexpl = vexpl, var.names = vn),
        file.path(CACHE_DIR, "varexpl_orth_draws.rds"))
message("  wrote varexpl_orth_draws.rds")

message("\nAll caches rebuilt in ", CACHE_DIR,
        "\nNext: Rscript code/03_figure6_structural.R  and  code/04_reduced_form.R")
