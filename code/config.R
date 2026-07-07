## =====================================================================
## config.R -- shared paths and settings for the factorBART JBES
## replication package.  Sourced by every 0x_*.R driver.
## =====================================================================
if (!exists("ROOT")) ROOT <- normalizePath(getwd())

DATA_DIR  <- file.path(ROOT, "data")
HELP_DIR  <- file.path(ROOT, "code", "helpers")
PKG_DIR   <- file.path(ROOT, "package")
## the next three may be redirected via environment variables (handy for testing)
CACHE_DIR <- Sys.getenv("REPL_CACHE", file.path(ROOT, "cache"))
DRAWS_DIR <- Sys.getenv("REPL_DRAWS", file.path(ROOT, "output", "draws"))
FIG_DIR   <- Sys.getenv("REPL_FIGS",  file.path(ROOT, "output", "figures"))
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DRAWS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,   recursive = TRUE, showWarnings = FALSE)

XLSX <- file.path(DATA_DIR, "CHKquarterlydata.xlsx")

## eight focus response variables shown in Figure 6 / the reduced-form exhibits
FOCUS  <- c("GDPC1","PCECC96","PNFIx","UNRATE","FGRECPTx","CPIAUCSL","FEDFUNDS","EBP")
FLAB   <- c(GDPC1="Real GDP", PCECC96="Real consumption",
            PNFIx="Bus. fixed investment", UNRATE="Unemployment rate",
            FGRECPTx="Federal tax receipts", CPIAUCSL="CPI inflation",
            FEDFUNDS="Fed funds rate", EBP="Excess bond premium")

message(sprintf("[factorBART replication] ROOT = %s", ROOT))
