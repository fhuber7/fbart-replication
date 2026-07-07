# Replication package — factorBART US IRF application (JBES revision)

This package reproduces the **main results of the US empirical application** in
*"Factor-augmented BART VAR with ex-post identification"* (factorBART):

* **Structural results — paper Figure 6.** Generalized impulse responses (GIRFs)
  of eight focus variables to large/small **financial-conditions** and
  **monetary-policy** shocks, positive vs. sign-flipped negative, with the
  **posterior asymmetry probability** `P = Pr(|positive response| > |negative
  response|)` printed above each panel.
* **Reduced-form findings.** The in-sample conditional-mean decomposition
  (nonlinear factor block `Λ_μ μ(x_t)` vs. linear block `A* x_t`) and the share
  of each variable's variance carried by the nonlinear block.

It ships the model as an installable R package (source + Windows + macOS
binaries) and small precomputed caches so every figure can be regenerated in
**seconds**, or from scratch by re-estimating the model.

---

## 1. Contents

```
Replication_JBES/
├── README.md
├── package/                      factorBART R package (v1.0.0)
│   ├── factorBART_1.0.0.tar.gz     source — installs on Windows AND macOS (pure R, no compiler)
│   ├── factorBART_1.0.0.zip        Windows binary
│   └── factorBART_1.0.0.tgz        macOS binary
├── data/
│   ├── CHKquarterlydata.xlsx       US quarterly dataset (22 series + sign restrictions)
│   └── prepare_data.R              builds the data matrix, transforms and sign matrix
├── code/
│   ├── config.R                    shared paths/settings (sourced by every script)
│   ├── helpers/                    data-transform helpers used during estimation
│   ├── 00_install.R                install CRAN deps + factorBART
│   ├── 01_estimate.R               estimate the model  → output/draws/
│   ├── 02_build_caches.R           draws → small summary caches → cache/
│   ├── 03_figure6_structural.R     cache → Figure 6 (GIRFs + asymmetry P)
│   └── 04_reduced_form.R           cache → reduced-form exhibits
├── cache/                          precomputed summaries (ship-with; ~2 MB)
└── output/
    ├── draws/                      posterior draws written by 01 (empty until you run it)
    └── figures/                    figures land here (ships with the reference output)
```

## 2. Requirements

* **R ≥ 4.1** (developed and tested under R 4.5.2).
* CRAN packages: `abind, coda, dbarts, dplyr, ggplot2, invgamma, MASS, Matrix,
  patchwork, readxl, reshape, stochvol, TruncatedNormal, tidyr`
  (installed automatically by `00_install.R`).
* **No compiler / Rtools needed.** `factorBART` is a *pure-R* package, so the
  source tarball installs identically on Windows and macOS.

## 3. Install

From a shell, `cd` into this folder, then:

```bash
Rscript code/00_install.R
```

This installs the CRAN dependencies and `factorBART` (from the bundled source
tarball). Alternatively, install the platform binary by hand inside R:

```r
# Windows:
install.packages("package/factorBART_1.0.0.zip",     repos = NULL)
# macOS:
install.packages("package/factorBART_1.0.0.tgz",     repos = NULL)
# any platform (recommended, no build tools required):
install.packages("package/factorBART_1.0.0.tar.gz",  repos = NULL, type = "source")
```

> The `.zip`/`.tgz` binaries were built on macOS; because the package contains
> no compiled code they install on either OS, but if a platform binary is ever
> refused simply use the **source tarball**, which always works.

## 4. Reproduce the results

Run everything from **this folder** (the scripts locate themselves; on Windows
use `Rscript.exe`).

### 4a. Quick — figures from the shipped caches (seconds)

```bash
Rscript code/03_figure6_structural.R   # Figure 6 (structural + asymmetry P)
Rscript code/04_reduced_form.R         # reduced-form exhibits
```

Outputs are written to `output/figures/` (which already contains the reference
copies for comparison):

| File | Paper exhibit |
|------|---------------|
| `irf_focus8_Financial_size_14_2x4.pdf` | Figure 6(a), large financial shock |
| `irf_focus8_Monetary_size_14_2x4.pdf`  | Figure 6(b), large monetary shock |
| `irf_focus8_*_size_5_2x4.pdf`          | small-shock (appendix) variants |
| `irf_focus8_legend.pdf`                | shared legend for Figure 6 |
| `insample_decomp_focus8.pdf`           | in-sample nonlinear vs. linear decomposition |
| `variation_expl.pdf`                   | variance share of the nonlinear block |

### 4b. Full — re-estimate the model, then rebuild everything

```bash
Rscript code/01_estimate.R          # posterior draws  → output/draws/
Rscript code/02_build_caches.R      # draws → cache/   (overwrites shipped caches)
Rscript code/03_figure6_structural.R
Rscript code/04_reduced_form.R
```

**Runtime.** `01_estimate.R` re-runs the Gibbs/BART sampler (default: 8
replicates, 5000 burn-in + 5000 saved draws each, GIRFs over 6 shock sizes × 4
shocks × 24 horizons with 50 Monte-Carlo paths). This is heavy — on the order of
a few hours of wall-clock on a multi-core machine (the 8 replicates run in
parallel via `parallel::mclapply`). Memory ≈ 2 GB per replicate.

**Quick smoke test** (minutes, *not* paper-accurate — just checks the pipeline):

```bash
REPL_MODE=smoke Rscript code/01_estimate.R
Rscript code/02_build_caches.R
Rscript code/03_figure6_structural.R
```

Individual knobs are environment-overridable, e.g.
`NSAVE=2000 NBURN=2000 NREPS=4 NCORES=4 Rscript code/01_estimate.R`.

## 5. The asymmetry measure (Figure 6)

For every posterior draw the response to a **positive** (contractionary) shock
and to an equal-sized **negative** (expansionary) shock are computed *jointly*.
Above each panel we report, at selected horizons,

```
P_h = Pr( | response to positive shock |  >  | response to negative shock | )
```

with **0.5 = symmetry**. Because the two responses come from the same draws they
are strongly positively correlated, so this probability is estimated far more
sharply than the overlap of the two marginal 16/84 bands suggests. The measure
is built in `02_build_caches.R` (`irf_asymmetry_measures.rds`) and drawn by
`03_figure6_structural.R`.

## 6. Data

`data/CHKquarterlydata.xlsx` (US quarterly, 1985Q1–2019Q4) provides 22 series,
their FRED-style transform codes, and the sign-restriction matrix for the four
structural shocks (Demand, Monetary, Supply, Financial). `data/prepare_data.R`
assembles the estimation inputs; sign restrictions follow Korobilis (2022).

## 7. Notes

* Posterior draws depend on RNG seeds (`set.seed(2024 + rep)`); tiny numerical
  differences across machines/BLAS are expected and do not affect the reported
  figures.
* The **full run in 4b overwrites the shipped `cache/`**. To keep the originals,
  copy `cache/` aside first, or redirect via `REPL_CACHE=/path Rscript ...`.
* `identify_factors()` performs the ex-post Hilbert-projection identification;
  the GIRFs are invariant to it (only the factor labelling changes).

## 8. Citation / contact

Florian Huber — <florian.huber@plus.ac.at>. Please cite the paper when using
this code or the `factorBART` package.
