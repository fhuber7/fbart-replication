## =====================================================================
## 00_install.R -- install CRAN dependencies and the factorBART package.
## Run once:  Rscript code/00_install.R   (from the Replication_JBES folder)
## =====================================================================
local({
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  root <- if (length(f)) normalizePath(file.path(dirname(f), "..")) else normalizePath(getwd())
  assign("ROOT", root, envir = .GlobalEnv)
})
source(file.path(ROOT, "code", "config.R"))

repos <- "https://cloud.r-project.org"

## CRAN packages used by the packaged replication pipeline
cran <- c("abind", "coda", "dbarts", "dplyr", "ggplot2", "invgamma",
          "MASS", "Matrix", "patchwork", "readxl", "reshape",
          "stochvol", "TruncatedNormal", "tidyr")
missing <- setdiff(cran, rownames(installed.packages()))
if (length(missing)) {
  message("Installing CRAN dependencies: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = repos)
} else message("All CRAN dependencies already installed.")

## factorBART (pure R -- no compiler / Rtools required).
## The source tarball installs identically on Windows and macOS.
tgz <- file.path(PKG_DIR, "factorBART_1.0.0.tar.gz")
message("Installing factorBART from: ", tgz)
install.packages(tgz, repos = NULL, type = "source")

## sanity check
suppressMessages(library(factorBART))
stopifnot(exists("bartfm"), exists("identify_factors"))
message("\nfactorBART ", as.character(packageVersion("factorBART")),
        " installed successfully. Setup complete.")
