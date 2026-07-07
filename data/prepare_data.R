## =====================================================================
## prepare_data.R -- build the US quarterly dataset, transform codes and
## sign-restriction matrix for the factorBART application.
##
## Expects the following to be defined BEFORE sourcing (see 01_estimate.R):
##   XLSX      : full path to CHKquarterlydata.xlsx  (set in config.R)
##   trans     : which transform-code column to use  ("stationary")
##   mod.size  : variable set to use                 ("large" -> 22 series)
##   fred.trans: transform helper (from code/helpers/aux_vec_irfs.R)
## Produces: Data (T x 22), transforms (length 22), sign.mat (22 x 4 x 2).
## =====================================================================
library(readxl)

data_raw   <- as.matrix(read_excel(XLSX, sheet = "N35_Tshorter", skip = 0))
data       <- data_raw[, 2:ncol(data_raw)]                     # drop the date column
transforms <- as.numeric(as.matrix(read_excel(XLSX, sheet = "N35_Tshorter_trans",
                                              skip = 0))[, trans])
class(data) <- "numeric"
data        <- data[-1, ]                                      # drop the trans-code row
Data        <- fred.trans(data, transforms)                   # transform to stationarity
names(transforms) <- colnames(Data)

if (mod.size == "large") {
  set.vars <- c("GDPC1","PCECC96","EBP","PRFIx","PNFIx","UNRATE","GCEC1",
                "FGRECPTx","PCEPILFE","CPIAUCSL","HRLYCOMP","OPHNFB","dtfp_util",
                "PAYEMS","UMCSENTx","INDPRO","HOUST","FEDFUNDS","TB3MS","GS5",
                "GS10","BAA")
} else {
  set.vars <- 1:ncol(Data)
}
Data       <- Data[, set.vars]
transforms <- transforms[set.vars]

S <- as.matrix(read_excel(XLSX, sheet = "Sign_N35_Tshorter", skip = 0))
rowS <- S[, 1]
S <- S[, 2:(ncol(S) - 1)]; class(S) <- "numeric"; rownames(S) <- rowS
S <- S[set.vars, ]

sign.mat <- array(NA, c(nrow(S), ncol(S), 2))
dimnames(sign.mat) <- list(rownames(S), colnames(S), c("lb", "ub"))
S[is.na(S)] <- 10^4
for (j in 1:nrow(S)) for (i in 1:ncol(S)) {
  if      (S[j, i] ==  1)    sign.mat[j, i, ] <- c(0,     Inf)
  else if (S[j, i] == -1)    sign.mat[j, i, ] <- c(-Inf,  0)
  else if (S[j, i] == 10^4)  sign.mat[j, i, ] <- c(-Inf,  Inf)
  else if (S[j, i] ==  0)    sign.mat[j, i, ] <- c(0,     0)
}
sign.mat <- sign.mat[, c("Demand","Monetary","Supply","Financial"), ]
