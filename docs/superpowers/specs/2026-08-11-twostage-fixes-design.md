# Two-stage pipeline fixes — design

- **Date:** 2026-08-11
- **Status:** draft for review
- **Branch strategy:** all work on a new branch `twostage-fixes`, created from
  `origin/twostage-estimator` (f041f89). Nothing lands on `main` until the
  validation plan passes. `twoStageSpiDE()` is unpublished/experimental, so
  breaking changes are made directly, without deprecation shims.

## Background

`twoStageSpiDE()` (branch `twostage-estimator`) estimates per-patient niche
slopes (stage 1) and contrasts them between conditions at patient level
(stage 2). The architecture is sound — stage-1 estimation error propagates
implicitly through the between-patient empirical variance — but a critical
review identified six problems. This design fixes all six without giving up
the property that makes the estimator attractive: validity comes from the
patient-level contrast, not from stage 1 being right.

The problems, in one line each:

1. The default stage-1 input (`"nbresid"`) is an affine transform of raw
   counts (intercept-only NB mean is constant per gene), so slopes are on the
   identity scale: condition main effects leak into the interaction contrast
   (slope ∝ patient baseline × β), and depth–density correlation inside a
   patient is unadjusted.
2. Slopes are marginal simple regressions per niche — correlated niches
   cross-contaminate, so triplet attribution is unreliable.
3. Stage 2 is an unmoderated Welch t-test (no shrinkage across genes, no
   precision weighting), and supplying `patient.covariates` silently swaps
   Welch for a homoscedastic `lm`.
4. `min.cells` dropout is informative: patients are dropped per index type,
   and cell-type abundance can correlate with condition — silently.
5. `patient == sample_id` is assumed; multiple cores per patient would
   re-enter stage 2 as fake patients.
6. The published benchmark cannot detect flaws 1–2 (its simulator never
   plants those failure modes), so fixes need new adversarial scenarios.

## Design

### Stage 1: SpaNorm-anchored one-step estimator

Input contract: a SpaNorm-normalised SPE whose stored fit (design split into
biology and library-size blocks, per-gene coefficients and dispersion ψ) is
retrievable from metadata. Per cell c and gene g, from the stored fit:

- linear predictor split: η̂ = η̂_bio + η̂_LS
- full fitted mean: μ̂ = exp(η̂) (depth included)
- working residual: z = (y − μ̂)/μ̂
- NB working weight: w = μ̂/(1 + ψμ̂)
- stage-1 response: **ε = η̂_bio + z** (≈ log y − η̂_LS, linearised at μ̂)

Rationale (established in review, locked):

- The residual and weights must use the **full** mean: it is the correct
  linearisation point (count → log scale conversion), it handles zeros
  (ε = η̂_bio − 1 at y = 0), and the weights down-weight cells where the
  linearisation is weakest. Residualising against a biology-only or LS-only
  mean reintroduces depth confounding or baseline leakage respectively.
- Adding η̂_bio back means, net, only the LS effect is removed — matching
  SpaNorm's design philosophy (the biology splines are a rough anchor, not
  truth) — and it exactly undoes the (1 − R²) attenuation from the smooth
  biology basis absorbing part of the niche signal.
- Slopes are on the log-mean scale, so condition main effects do not leak
  into the slope contrast at any expression level, and the LS fix inherits
  SpaNorm's anchoring (no naive-depth-adjustment signal loss).

Estimation: within each (sample, index cell type) subset, one **joint**
weighted least-squares fit of ε on all centred niche columns (self-niche
excluded), weights w. Joint, not marginal: the add-back re-admits smooth
spatial trends, and joint estimation forces such a trend to compete across
correlated niche columns instead of loading fully on every one of them —
the add-back and joint estimation ship together. Outputs per subset: the
slope vector β̂ and its covariance (from the inverted weighted Gram matrix
times the working dispersion φ̂, estimated per gene per subset from the
weighted squared residuals), giving a per-coefficient variance v.
Degenerate/collinear niche columns within a subset are dropped by a rank
guard (generalising the current `ss > 1e-8` check); dropped columns yield NA
for that subset.

Joint estimation and the (slope, v) extraction apply to **every** stage-1
path, because stage 2's weights require v: `"ols"` becomes an unweighted
joint OLS of log-CPM on the centred niches, and `"nb"` reads the slopes and
v from the GLM's coefficients and covariance. Only the response and weights
differ per path.

Diagnostic (`"spanorm"` path only): per subset, R² of each niche column on
the biology basis — reported so users can see, per patient, how much
smooth-trend overlap the add-back is exposed to (and how much attenuation
the conservative variant would have suffered).

Option surface: `stage1 = c("spanorm", "ols")`, default `"spanorm"` (clear
error naming the SpaNorm call when no stored fit is present); `"ols"`
(log-CPM) kept as the no-dependency fallback; `"nb"` (per-subset NB GLM)
kept as the slow reference implementation for validating the one-step
estimator. **`"nbresid"`, `"analytic"`, and `"auto"` are removed outright**
— no deprecation messages (unpublished code). `.looksLikeCounts()` remains
only as an input-validation guard for the paths that need counts.

The response is selected by `epsilon = c("addback", "residual")`: the
default `"addback"` is ε = η̂_bio + z; `"residual"` (ε = z, attenuated but
smooth-trend-adjusted) is kept for sensitivity analysis.

### Between stages: patient pooling

New `patient` argument (colData column, default = `sample_id`). Stage 1
runs per sample (spatial fields are per core). Per-patient slopes are the
precision-weighted (1/v) average of that patient's per-core slopes — cores
share the patient's β, so fixed-effect pooling is correct here — with pooled
variance updated accordingly. Condition constancy is checked per patient.
This step precedes all stage-2 computation, including τ² estimation.

### Stage 2: weighted, moderated limma

Per (index, niche) pair: assemble the genes × patients slope matrix and the
per-slope variances v; fit `limma::lmFit` with design
`~ condition + patient.covariates` and observation weights **1/(v + τ̂²)**;
`eBayes(robust = TRUE)`; report the condition coefficient's t and p into the
existing tidy results table (same columns), flat BH per triplet as now.

τ̂² per (index, niche): DerSimonian–Laird moment estimate per gene, then a
robust pool (median across genes), floored at 0. Chosen over 1/v-only
weights deliberately: when between-patient variation dominates, 1/(v + τ²)
degrades toward equal weights (Welch-like validity preserved — avoiding the
pseudobulk pitfall of over-weighting cell-rich patients); when v dominates
(small/imbalanced subsets) the weights do real work. This replaces both the
Welch path and the covariate `lm` path with one machinery.

New package dependency: `limma` (Bioconductor).

### Dropout diagnostics

`min.cells` default stays 30 (slope estimability needs within-subset niche
variance regardless). Added: a per-index patient inclusion table in the
diagnostics, plus a warning when inclusion is condition-associated (Fisher
test, warn at p < 0.05). Near-floor subsets now enter stage 2 with honestly
large v, so the cliff is also less consequential.

### Results object

`SpiDEResults` gains a `diagnostics` slot (list, default empty, tolerated as
NULL for old objects): the R² table, the inclusion table, and the τ̂² table.
No other slot changes; the results table schema is unchanged.

## Data flow

```
SpaNorm-fitted SPE
  └─ stage 1 (per sample × index): joint WLS of ε on centred niches, weights w
       → β̂, v, R² diagnostic
  └─ patient pooling: precision-weighted average of per-core slopes
  └─ τ̂² per (index, niche): DL moment estimate, median-pooled, floored
  └─ stage 2: limma lmFit(~ condition + covariates, weights 1/(v + τ̂²)),
       eBayes(robust) → results table (unchanged schema), flat BH
  └─ diagnostics: R², inclusion, τ̂² tables + condition-imbalance warning
```

## Validation plan (acceptance criteria)

Reuses the local sweep infrastructure in
`/Users/uqdbhuva/R_packages/spiDE-twostage-bench` (runner/driver/aggregate;
worktree repointed at `twostage-fixes`). Four scenario families:

1. **Removable-interaction null** (targets fix 1): condition main effect on
   baseline + identical nonzero niche slope in every patient; no true
   interaction. Old default fails by construction; new stage 1 must hold
   type-I ≈ 0.05.
2. **Depth-coupling null** (targets fix 1): within-patient depth–density
   correlation whose strength differs by condition; no expression effect.
   Must hold type-I ≈ 0.05.
3. **Attribution scenario** (targets fix 2): two strongly correlated niche
   densities, true effect through one only. The true niche's discoveries
   must dominate the decoy's (marginal slopes light up both).
4. **Regression grid**: the original null (2 layouts × S) and power sweeps —
   type-I calibration preserved; TPR at S ≥ 24 not worse than the current
   estimator (moderation + weighting should improve it below the crossover).

## Out of scope

- Multiplicity beyond flat BH (the pre-specified-hypothesis guidance
  stands).
- Cross-bandwidth combination for the two-stage path.
- Any change to the one-stage `fitSpiDE()` pipeline.

## Decision log

- Stage-1 construction: full-mean working residual + η̂_bio add-back
  (user-approved; `ε = z` kept as conservative option).
- Joint niche regression coupled to the add-back (design requirement).
- Stage-2 weights: 1/(v + τ̂²) (user-approved over unweighted and 1/v).
- Separate branch `twostage-fixes`; no deprecation shims for removed
  `stage1` options (user-directed).
