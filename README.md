# spiDE

<!-- badges: start -->
<!-- badges: end -->

**spiDE** finds neighbourhood-dependent differential expression in spatial
transcriptomics. Within an *index* cell type, it asks whether a gene's
expression changes with the local density (the *niche*) of another cell type,
within patients, and whether that dependence differs between two conditions.
Patients are the unit of replication throughout: a triplet's evidence is how
its niche slope varies between patients, never the number of cells.

1. **Niches** (`buildNiches()`): per sample, a Gaussian kernel density of every
   cell type at every cell.
2. **Fit** (`fitSpiDE()`): for every index cell type, one negative binomial GLM
   per gene on the niche densities with an intercept per patient, fitted with
   [SpaNorm](https://github.com/bhuvad/SpaNorm)'s Newton solver. Each niche
   coefficient is a within-patient slope.
3. **Test** (`testSpiDE()`): the **pooled** test (is the slope non-zero,
   consistently across patients?) always, and the **condition-specific** test
   (does it differ between conditions?) when given a condition.

Two engines estimate the between-patient error. The **slopes** engine (the
default) estimates each patient's own slopes and combines them with weighted
`limma`; the **sandwich** engine fits the condition-specific model and uses a
patient-clustered CR2 sandwich. `testNicheAbundance()` asks the different,
between-patient question of whether patients with more of a niche type express
genes differently.

## Installation

spiDE needs SpaNorm (>= 1.7.14).

```r
# install.packages("BiocManager")
BiocManager::install("BhuvaLab/spiDE")
```

## Quick start

```r
library(spiDE)
data(toySpiDE)

res <- spiDE(toySpiDE, condition = "condition", sigma = 30)
results(res, test = "pooled")      # niche-dependent expression, across patients
results(res, test = "condition")   # its difference between the conditions
```

The vignettes cover the walk-through (`vignette("spiDE")`), the model
(`vignette("spiDE-model")`) and what its calibration rests on
(`vignette("spiDE-calibration")`).

## Earlier versions

spiDE 0.99.22 and earlier fitted a joint mixed-effects model, which was not
calibrated. It is archived, with its tests and documentation, as the research
package `spiDEmixed` (`research/mixed` in this repository); objects saved by it
are read with `spiDEmixed::readSpiDE()`.
