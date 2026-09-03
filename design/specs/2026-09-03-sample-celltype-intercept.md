# A (sample × cell type) random intercept in the niche design

**Status:** proposed; the research equivalent is validated on all twelve raw shuffle grids (2026-09-03). Implementation next. Evidence in
`research/fdr-ordering/FINDINGS.md` (2026-09-03 entry, §7) and `REPORT.md` §5d.

## The defect

`.buildNicheDesign()` with `random = "intercept"` adds one ridge-penalised
intercept per sample, shared across cell types (`.buildRandomEffects()`,
`re_group = "SampleInt"`). There is no (sample × cell type) intercept. The
tested `CellType:condition:niche` slopes are therefore estimated from the
*total* covariance of niche density and expression among the cells of that
type — within-sample **and** between-sample — and the between-sample part is a
composition effect: samples whose type-`k` cells sit, on average, in a denser
niche of type `n` also have a different mean expression of the gene in type
`k`. That is a patient-level association with S = 55 units, not
neighbourhood-dependent DE, and it is reported as the latter with a cell-level
standard error.

Measured on the real YTMA cohort at the converged fit
(`research/fdr-ordering/R/perm_decomposition.R`, `perm_centred.R`): with the
fit held fixed and the niche rows permuted within (sample, cell type), the
permutation distribution of each tested `t` has the HC0 spread (0.8–1.0) but a
per-column **mean** of up to ±3 that equals the real-data `t` column by column
(CDV3 Tumor: r = 0.66 over 132 columns). The free shuffle preserves each
group's mean niche density, so it preserves the confound; every shuffle grid
carried it, which is why real data and "null" were indistinguishable and why
twelve variance-side candidates failed — the excess is a bias of the estimand,
not a variance. It scales with cells per index type (Tumor 1.36, types under
3,000 cells ~1.0) and is confined to bright genes.

## The fix

Add a second penalised block to `.buildRandomEffects()`:

```r
Zct <- stats::model.matrix(~ 0 + smp:ct)            # drop all-zero columns
colnames(Zct) <- paste0("Sample", ..., ":CellType", ...)
re_group <- c(re_group, rep("SampleCellTypeInt", ncol(Zct)))
```

- Tagged `Random` in `covtype`, so `.testedCols()` and the inference are
  untouched.
- Its own variance component in `.fitNBmixed()`: `tau2` is already estimated
  per `re_group` (`SampleInt`, `SampleSlope`); add `SampleCellTypeInt`.
- The Satterthwaite machinery reads `re_group`, so the `ResponseNiche` df
  stays a within-sample contrast.
- With the block present, every `CellType:niche` and `CellType:condition:niche`
  slope is a within-(sample, cell type) slope (Frisch–Waugh). The research
  equivalent — centring the niche-dependent columns within (sample, cell type),
  `CONVNULL_CENTRE=1` in `converged_null.R` — is exact for the slopes but does
  not model the group means; the package should add the intercepts.
- `random = "slope"` does **not** address this: per-sample slopes on the
  `CellType:niche` bases leave the group means untouched.
- Cost: up to S × K extra columns (55 × 12 = 660 here, fewer where a type is
  absent from a sample), so the per-gene gram grows from 345² to ~1000²
  (≈ 2.5 s per gene on 4 threads). The per-gene Newton polish that the fit
  defect needs (`converged_null.R`) is the natural place to pay it.

`checkSample()` currently rejects sample-constant covariates under
`random != "none"` because the per-sample intercept absorbs them; a (sample ×
cell type) intercept additionally absorbs (sample × cell type)-constant
covariates. Document it.

## What to expose

The between-sample association is worth having, on its own terms: a
patient-level regression of mean expression in type `k` on mean niche density
of type `n` around type `k`, S units, `limma` or plain `lm`, with the
patient-level covariates `fitSpiDE()` rejects. It answers "do patients with
B-cell-rich tumour compartments express X differently?", which is a different
question from spiDE's, and must not be reported as niche-dependent DE.

## Validation (done for the research equivalent)

`convnull_centred.sbatch` (job 27963912): the 12-grid converged null with every
niche-dependent column centred within (sample, cell type), scored with
`CONVNULL_SUFFIX=_centred R/score_convnull.R`. Median per-gene `sd(null t)`,
ML `psi` + Pearson:

| band | production | converged | centred, block | centred, free |
|---|---|---|---|---|
| below the top 5% | 0.929 | 1.002 | 0.988 | 0.984 |
| top 5%, aveLogCPM < 8.2 | 1.078 | 1.086 | 1.019 | 0.984 |
| 8.2–9 | 0.878 | 1.252 | 1.095 | 0.977 |
| 9–9.9 | 0.760 | 1.395 | 1.163 | 0.988 |
| ≥ 9.9 | 0.221 | 1.423 | 1.188 | 0.963 |
| gradient | 1.149 / 1.198 | 1.230 / 1.132 | 1.037 | 1.000 |
| `|t| > 4.89` exceedances | 555 / 95 | 1,091 / 133 | 138 | 0 |

Under `free` shuffles: flat, no gradient, no extreme tail, per-gene spread
0.92–1.05. Under `block`: a residual 1.10–1.19 in the ~100 brightest genes
(the spatially-smooth-covariate component). Real-data calls beyond
`|t| > 4.89` on the 769-gene panel: 182 → 87.

## Validation still owed (package)

1. The package implementation above, checked against the centred research
   fit gene by gene (the two should agree on the tested slopes up to the
   shrinkage of the new intercepts).
2. The simulation study (`research/reports/benchmarks/spiDE-simulation.Rmd`)
   re-run for type-I and power — the simulation has no between-sample
   composition effect, so power should be unchanged and the fix costs nothing
   there.
3. Real-data discoveries on the fixed design against `block` shuffles, scored
   with `.claude/skills/calibration-check/`; `block` is the conservative null
   because it keeps the spatial component the centred null still shows.
