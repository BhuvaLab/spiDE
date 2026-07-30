# `random="nested"`: a deterministic fixed-effects calibrated path

**Branch:** `satterthwaite-df` (or a fresh `nested-fixed`) · **Date:** 2026-07-30 · **Status:** design (awaiting review)

## Context

`fitSpiDE(random=...)` currently offers `"none"` (fixed, cell-level Wald —
documented anti-conservative baseline), `"intercept"`, and `"slope"` (mixed, ridge
= random effects, estimated via a Schall/PQL loop). The mixed path is calibrated
but its dominant cost is the iterative variance-component loop (re-fits the whole
gene set ~4–11×), and — per the statistician review — the MoM variance estimate is
statistically inefficient.

Two results from `analysis/statistician-review/` motivate a new option:

1. A **nested contr.sum fixed model, tested against the between-sample stratum**,
   is calibrated (null type-I ≈ 0.05) and **robust to a 1×–10× cell imbalance**,
   whereas the same nesting tested at the *cell* SE stays anti-conservative
   (≈ 0.8, regardless of df — the df is a red herring; the SE is the issue).
2. The only statistical advantage of the mixed model over this fixed-nested test
   — shrinkage/pooling of the between-sample variance across genes — **vanishes
   with sample size**: the power gap is ~10% at S=12, ~5% at S=24, and **~0.5% at
   S=55** (the real cohort scale), at every sparsity tier, with both calibrated.

So at realistic cohort sizes a **deterministic, single-fit** fixed-nested engine
matches the mixed model on power and calibration while avoiding the PQL loop. This
spec adds it as `random="nested"`, leaving `random="none"` (the documented
back-compat baseline) and the mixed paths untouched. Recommendation for users:
`"nested"` for large S (efficient), `"intercept"` for small S (shrinkage helps).

## Key insight (what makes it calibrated)

The response **main effect** η (condition) is a pure between-sample contrast; its
honest error is between-sample (S−2 df). The **niche interactions** β
(condition×celltype×niche) multiply a *within*-sample covariate, and with
τ²_slope ≈ 0 they are within-sample contrasts whose honest error is
within-sample (large df). A single test can't use one df/stratum for both — this
is the same per-coefficient logic the Satterthwaite work established, realised
here with **fixed** strata instead of variance components.

Crucially, adding the contr.sum columns and reading η's SE off the joint model's
(cell-level) Wald covariance does **not** calibrate it (empirically 0.835). The
calibrated between-sample SE is the **samples-as-clusters aggregated score
variance** — mathematically a cluster-robust (CR) variance for the response
block, which the contr.sum framing and the cluster-robust framing both reduce to.

## Design

### 1. Design matrix (`R/design.R`, `.buildNicheDesign`)
- New `random = "nested"`: append **contr.sum-nested sample-intercept columns**
  — one contr.sum basis of the samples *within each condition group* (S−2 columns
  total), so `resp`/η stays estimable without collinearity and without a penalty.
- Tag them `covtype = "Random"` with `re_group = "SampleNested"` and **penalty 0**
  (they are fixed, not ridge-penalised). Reuse the existing `.buildRandomEffects`
  plumbing shape but with contr.sum coding and zero penalty.

### 2. Fit (`R/fitSpiDE.R`, `.fitOneBandwidth`)
- `random = "nested"` → a **single** `SpaNorm::fitNB` with these columns included
  and `lambda.a = 0` on them. **No Schall loop, no τ² estimation, no `re.prop`
  subsampling** — deterministic. Populate `@re_group` (`"SampleNested"`),
  `@penalty` (0 on those columns), and the between-sample `@df` (below); `@tau2`
  stays `NULL`.

### 3. Inference (`R/inference.R`, `.blockedInference`) — the crux
Per gene, per tested response column, a **between-sample (samples-as-clusters)
variance** for the response block, with **per-coefficient df**:
- Reuse the per-gene working weights `w_{c,g}` and residuals already formed in the
  block loop. For the tested sub-design `Wsub` (Response + ResponseNiche columns),
  form the per-cell score `s_{c} = Wsub_c · w_{c,g} · e_{c,g}` and **aggregate to
  the sample level** (`Σ_{c∈s}`), giving a `S × k` cluster-score matrix `G_s`.
- The cluster-robust covariance of the response block is
  `V_CR = A⁻¹ (Σ_s G_s' G_s) A⁻¹`, `A = Wsub' diag(w) Wsub`, with a **CR2**
  small-sample bias correction (per-cluster leverage adjustment). `se = sqrt(diag(V_CR))`.
- **Per-coefficient df** (stored in `@df`, the named vector the Satterthwaite work
  already introduced): `covtype == "Response"` (η) → **S − 2** (between-sample);
  `covtype == "ResponseNiche"` (β) → the within-sample residual df (large; a
  Satterthwaite-style df on the cluster variance, which for β returns ≈ the cell
  residual because β carries no between-sample structure). Both consumed through
  the existing `.ptByCol()`.
- This is O(cells) aggregation per gene on top of the existing block loop —
  cheap, and it slots into the same GPU/CPU block machinery.

**Reuse:** the `@df` named-vector slot + `.ptByCol()` (from the Satterthwaite
feature) already carry per-coefficient df through `.blockedInference` and
`.nicheRecords`; `"nested"` populates the same slot, so downstream (FDR,
`results()`) needs no change.

## Validation (the plan's oracle)
- **Reproduce calibration:** on the `analysis/statistician-review` null, `"nested"`
  gives η null type-I ≈ 0.05 (not 0.835) and β null type-I ≈ 0.05, robust across
  the 1×–10× imbalance sweep. This is the load-bearing test — the 0.835-vs-0.06
  lesson means the SE computation must be verified, not assumed.
- **RE-equivalence at large S:** `"nested"` matches `random="intercept"` +
  `df.method="satterthwaite"` on power/calibration to within a few % at S≈24–55
  (the shrinkage study), and is cheaper (no PQL loop — assert far fewer `fitNB`
  calls / wall-clock).
- **Cross-check η's between-sample SE** against the standalone sample-mean/ANOVA-F
  computation validated in `imbalance.R` (they must agree).
- **`random="none"`/`"intercept"`/`"slope"` untouched** (byte-identical).

## Scope / risks
- **New mode only**; does not change `random="none"` (the documented baseline) or
  the mixed paths. `df.method` is ignored for `"nested"` (it computes its own
  per-coefficient df).
- **Small-S caveat documented:** recommend `"intercept"` for S ≲ 20 (shrinkage),
  `"nested"` for larger cohorts (efficiency); both calibrated.
- **Main risk — the CR2 small-sample correction.** A naive (CR0) cluster variance
  is anti-conservative at small S; CR2 + the S−2 (η) / Satterthwaite (β) df is what
  calibrates it. The plan validates against the 0.06 target and, if the hand-rolled
  CR2 proves finicky, falls back to the `clubSandwich` reference (test-only) to pin
  correctness — the analogue of the Satterthwaite work's `lmerTest` oracle.
- Vignette: add `"nested"` to the model vignette's mode guidance (efficient large-S
  option; the nesting aside already sets it up) + NEWS.

## Out of scope
- Cluster-robust as the *default* for other modes; changing `random="none"`;
  shrinkage/pooling for `"nested"` (its defining trade-off is *no* pooling).
