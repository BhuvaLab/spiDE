# spiDE

<!-- badges: start -->
<!-- badges: end -->

**spiDE** identifies context-specific, neighbourhood-dependent differential
expression in spatial transcriptomics data. Within an *index* cell type, it
tests how gene expression changes with an experimental *condition* as a function
of the local density (the *niche*) of surrounding cell types.

The method:

1. builds per-cell **niche** covariates from Gaussian kernel density estimates of
   each cell type at multiple spatial bandwidths (`buildNiches()`);
2. fits a per-gene **negative binomial GLM** over a design containing the
   three-way `cell type : condition : niche` interactions, using the
   [SpaNorm](https://bioconductor.org/packages/SpaNorm) `fitNB()` engine
   (`fitSpiDE()`);
3. **tests** the neighbourhood interactions with Wald statistics combined across
   correlated covariates (Brown's method) and bandwidths (Cauchy combination),
   under a hierarchical (gene → index cell type → niche cell type) FDR
   (`testSpiDE()`).

## Installation

spiDE depends on SpaNorm (>= 1.7.4), which exposes the negative binomial fitting
engine (`fitNB`) and the `calculateMu` / `invert_mat` helpers.

```r
# install.packages("BiocManager")
BiocManager::install("bhuvad/spiDE")
```

## Quick start

```r
library(spiDE)
data(toySpiDE)

res <- spiDE(toySpiDE, condition = "condition", covariates = "Age")
head(results(res))
```

See the vignette (`vignette("spiDE")`) for a full walk-through.
