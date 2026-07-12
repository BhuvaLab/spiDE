# Mixed-effects (PQL) speed-ups design

Date: 2026-07-11
Branch: `mixed-effects-ridge`

## Problem

The mixed-effects fit (`fitSpiDE(..., random = "intercept" | "slope")`) is ~8×
slower than the fixed-effects fit. Measured on the toy data (480 cells, 20
genes, 6 samples, σ=20): fixed 7.8 s, intercept 61.3 s, slope 57.2 s.

The entire penalty is the Schall/PQL outer loop in `.fitNBmixed()`
(`R/fitSpiDE.R`), which calls `SpaNorm::fitNB(Y, W, ...)` over the **whole gene
set** once per iteration, for up to `re.maxit` (default 10) iterations. Two
compounding wastes:

1. Every iteration re-fits all genes on all cells, but the variance components
   `tau2` are a *single shared scalar per random-effect group* — a population
   quantity that does not need every cell to estimate.
2. Every iteration re-runs the full NB dispersion moderation
   (`edgeR::estimateDisp`, robust/tagwise), although only the ridge penalty
   changes between iterations and the dispersion barely moves after iteration 1.

Observed secondary issue: with `random = "slope"` the loop never converges
(`SampleSlope` tau2 drifts 0.63 → 0.15, still ~8 %/iter at the cap), so all 10
full fits are always paid for.

Separately, `.fitOneBandwidth()` eagerly densifies the entire counts matrix
(`as.matrix(Y)` at `R/fitSpiDE.R:132`) to compute a per-gene loglik that
`.blockedInference()` unconditionally overwrites — dead work that also violates
the "never densify `Y`" invariant.

## Goals

Trade a little `tau2` accuracy for speed (the latitude the user granted, and the
same strategy `fitNB` already uses to subsample cells for its dispersion step),
bringing the mixed penalty from ~8× back toward ~1×, and quantify the trade-off.

## Design

### 1. Stratified cell subsampling inside the PQL loop (`re.prop`)

New helper `.stratifiedCellIdx(cell_type, sample, prop, min.cells)` returns a
logical vector over cells. For each `(cell_type, sample)` stratum of size `n`:

```
n_sampled = min(n, max(ceil(prop * n), min.cells))
```

so `n ≤ min.cells` → all cells; `n = 500, prop = 0.1` → 100; `n = 2000, prop =
0.1` → 200. The index is drawn **once** before the loop and reused every
iteration (stable `tau2` trajectory; only the penalty changes between steps).
No internal seed is set — reproducibility is the user's responsibility via an
external `set.seed()`, per R convention.

The sampled `idx` is passed to `SpaNorm::fitNB(Y, W, idx = idx, ...)` (fitNB's
existing cell-subsampling argument). The Schall update (representative weights,
information matrix, BLUP pool) is computed on the **same** sampled cells
(`W[idx, ]`, `alpha`) so the update is internally consistent with the fit.

### 2. One dispersion iteration inside the loop, full dispersion at the end (`re.maxit.psi`)

Inside the loop, `fitNB` is called with `maxit.psi = re.maxit.psi` (default
`1L`). After the loop converges (or hits `re.maxit`), a single **final fit** is
run on **all cells** with **full dispersion** (fitNB's default `maxit.psi`, or a
user value forwarded via `...`). The final fit supplies the `alpha`, `psi`,
`gmean` used for inference; the loop only supplies the converged `tau2` / penalty
/ df. `maxit.psi` is stripped from `...` for the internal call so `re.maxit.psi`
wins there, and forwarded to the final fit.

### 3. New / changed API

On `fitSpiDE()` (and via `...` on `spiDE()`; documented through
`@inheritParams`):

| param          | default | meaning |
|----------------|---------|---------|
| `re.prop`      | `0.1`   | cell proportion sampled per cell type × sample for the PQL loop. `1` = all cells (reproducible, old behaviour). |
| `re.maxit.psi` | `1L`    | dispersion iterations inside each PQL step. |
| `re.min.cells` | `100L`  | per-stratum floor. |

Defaults change the mixed-fit numbers (this is the opt-in-by-default speed-up the
user approved). `random = "none"` is untouched. If the benchmark shows the
`re.prop = 0.1` RNG variation is large, the default reverts to `re.prop = 1`.

### 4. Convergence metric (light touch)

Switch the `.fitNBmixed()` convergence test to a relative change on `log(tau2)`
(natural for a scale parameter). No Aitken acceleration — each iteration is now
cheap, and perturbing the fixed point risks the calibration tests. `re.maxit` /
`re.tol` stay user-tunable.

### 5. Finding 4 — stop densifying `Y`

Replace the `as.matrix(Y)` loglik at `R/fitSpiDE.R:132` with a gene-chunked
computation (reuse `.chunkGenes`), so the fit stage never realises the whole
counts matrix. Benefits fixed and mixed paths.

## Benchmark vignette

New `vignettes/spiDE-mixed-benchmark.Rmd` (BiocStyle, matching
`spiDE-cauchy-vs-brown.Rmd`), on a scaled-up `.toyClustered()` dataset (more
cells/genes so timings separate):

- **Time vs accuracy sweep:** `re.prop ∈ {1, 0.5, 0.25, 0.1, 0.05}` ×
  `re.maxit.psi ∈ {1, full}`, each compared against the `re.prop = 1`,
  full-dispersion reference. Accuracy metrics: ResponseNiche t-stat correlation,
  combined-p concordance (−log10), `tau2` values, planted-signal recovery.
- **RNG variation:** `re.prop = 0.1` over several seeds → spread of `tau2` and
  response p-values (the evidence for keeping sampling on by default).
- Figures rendered inline; a self-contained HTML Artifact summarises time-vs-prop,
  accuracy-vs-prop, and the RNG spread for sharing.

## Testing

- `random = "none"` results unchanged (existing test).
- New: `re.prop = 1, re.maxit.psi = <full>` reproduces the pre-change mixed fit
  (within tolerance) — guards the refactor.
- New: `.stratifiedCellIdx()` honours the min-cells floor, the proportion, and
  the cap at stratum size; every stratum is represented.
- Existing calibration tests (`test-mixedEffects.R`) still pass with the new
  defaults (this is the accuracy gate).

## Out of scope

- Aitken/step acceleration of the Schall loop.
- Gene subsampling (cell subsampling chosen instead; keeps every sample/celltype
  represented).
- Any change to `SpaNorm`.
```

## Non-goals

Warm-starting `fitNB` across PQL iterations (needs a SpaNorm change).
