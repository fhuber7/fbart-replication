## =====================================================================
## 01_estimate.R -- estimate the factor-BART VAR on the US quarterly data
## and identify the factors ex post (Hilbert projection).  Produces the
## posterior-draw files consumed by 02_build_caches.R.
##
##   Rscript code/01_estimate.R              # FULL run (paper settings; hours)
##   REPL_MODE=smoke Rscript code/01_estimate.R   # quick smoke test (minutes)
##
## Individual knobs can also be overridden, e.g.
##   NSAVE=2000 NBURN=2000 NREPS=4 NCORES=4 Rscript code/01_estimate.R
## =====================================================================
local({
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  root <- if (length(f)) normalizePath(file.path(dirname(f), "..")) else normalizePath(getwd())
  assign("ROOT", root, envir = .GlobalEnv)
})
source(file.path(ROOT, "code", "config.R"))

Sys.setenv(OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)
suppressMessages({ library(parallel); library(factorBART); library(readxl) })

## ---- data ------------------------------------------------------------
source(file.path(HELP_DIR, "aux_vec_irfs.R"))   # provides fred.trans() etc.
source(file.path(HELP_DIR, "plot.aux.R"))
trans    <- "stationary"
mod.size <- "large"
source(file.path(DATA_DIR, "prepare_data.R"))   # -> Data, transforms, sign.mat

## ---- MCMC / model settings ------------------------------------------
SMOKE <- identical(tolower(Sys.getenv("REPL_MODE", "full")), "smoke")
geti  <- function(v, d) as.integer(Sys.getenv(v, d))
if (SMOKE) {
  NSAVE<-geti("NSAVE",100); NBURN<-geti("NBURN",100); NOUT<-geti("NOUT",20)
  NREPS<-geti("NREPS",2);   N_PATHS<-geti("NPATHS",5)
  message(">>> SMOKE TEST settings (results are NOT paper-accurate) <<<")
} else {
  NSAVE<-geti("NSAVE",5000); NBURN<-geti("NBURN",5000); NOUT<-geti("NOUT",500)
  NREPS<-geti("NREPS",8);    N_PATHS<-geti("NPATHS",50)
}
NCORES <- geti("NCORES", min(NREPS, parallel::detectCores() - 1))
P_LAG <- 2; Q_MU <- 8; Q_Q <- ncol(sign.mat); NHOR <- 24
SL_SHOCK  <- c(1, 2, 3, 4)                       # Demand, Monetary, Supply, Financial
SIZE_GRID <- c(-25, -14, -5, 5, 14, 25)

cat(sprintf("[01_estimate] mode=%s NREPS=%d NSAVE=%d NBURN=%d NOUT=%d NPATHS=%d cores=%d\n",
            ifelse(SMOKE,"smoke","full"), NREPS, NSAVE, NBURN, NOUT, N_PATHS, NCORES))

## ---- one replicate ---------------------------------------------------
run_one <- function(rep_idx) {
  fn <- file.path(DRAWS_DIR, sprintf("US_baseline_rep%d.rds", rep_idx))
  if (file.exists(fn)) return(invisible(fn))
  set.seed(2024L + rep_idx)

  Yraw <- Data
  Y.sd <- apply(Yraw, 2, sd); Y.m <- apply(Yraw, 2, mean)
  Yraw <- scale(Yraw)
  Yraw <- ts(Yraw, start = c(1977, 1), frequency = 4)
  Yraw <- window(Yraw, start = c(1985, 1), end = c(2019, 4))
  colnames(Yraw) <- colnames(Data)
  sums.1 <- sort(c(which(transforms == 3), which(transforms == 5), which(transforms == 6)))

  t0 <- Sys.time()
  fit <- factorBART::bartfm(
    Yraw = Yraw, nsave = NSAVE, nburn = NBURN, nout = NOUT,
    p = P_LAG, Q.mu = Q_MU, Q.q = Q_Q,
    prior = "HS", sign = TRUE, sign.mat = sign.mat,
    trans.1 = sums.1, trans.2 = NULL, Y.sd = Y.sd, Y.m = Y.m,
    VAR = TRUE, IRF = TRUE, scens = TRUE,
    nhor = NHOR, sl.shock = SL_SHOCK, size.grid = SIZE_GRID,
    shrink.load = "rowwise", ident = "unident", column.shrink = TRUE,
    n.paths = N_PATHS)
  ## ex-post identification (Hilbert projection); IRF is invariant under it
  id <- factorBART::identify_factors(fit, rotate = TRUE, rescale = TRUE, signident = TRUE)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("[rep %d] %.1f min | IRF dim %s\n", rep_idx, wall,
              paste(dim(fit$IRF), collapse = "x")))

  out <- list(IRF = fit$IRF, A = id$A, Lambdamu = id$Lambdamu, mu = id$mu,
              var.expl = fit$var.expl, X = fit$X,
              Y.sd = Y.sd, Y.m = Y.m, var.names = colnames(Yraw))
  saveRDS(out, fn); invisible(fn)
}

cat("[01_estimate] running", NREPS, "replicate(s) on", NCORES, "core(s) ...\n")
t_all <- Sys.time()
invisible(mclapply(seq_len(NREPS), run_one, mc.cores = NCORES))
cat(sprintf("[01_estimate] done in %.1f min. Draws in: %s\n",
            as.numeric(difftime(Sys.time(), t_all, units = "mins")), DRAWS_DIR))
cat("Next: Rscript code/02_build_caches.R\n")
