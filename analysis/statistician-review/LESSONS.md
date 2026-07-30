# Statistician review — lessons learnt (and the shelved `random="nested"` design)

This folder holds a self-contained investigation prompted by a statistician's
review of the spiDE mixed-effects model (the `spiDE-model` vignette). It settled
several questions about *why* the mixed-effects correction is needed and whether a
deterministic fixed-effects alternative is worth building. **The fixed-effects
`random="nested"` mode was designed but deliberately NOT built** — this file
records everything so it can be picked up later.

## TL;DR

- The recommended default is **`random="intercept"` + `df.method="satterthwaite"`**.
- **Random slopes are inert on real data** (τ²_slope ≈ 0) — do not use `"slope"` by
  default; the niche interactions need the *per-coefficient df*, not a random slope.
- A **fixed nested (contr.sum) + between-sample-stratum** test is calibrated and
  robust to cell imbalance, and is ~equivalent to `random="intercept"` — but its
  only advantage (avoiding the PQL loop) buys nothing statistically once S is
  moderately large, so it was **shelved** rather than built.

## The scripts (reproduce everything)

| script | what it shows |
|---|---|
| `re-vs-fixed.R` | the reviewer's exact requests: how per-sample effects are generated; zero-mean-within-group vs unequal group means; the nested contr.sum model (cell-SE vs between-sample stratum vs RE); the slope stress test; MoM vs REML `τ²` |
| `imbalance.R` | the nested-fixed vs RE comparison across a **1×–10× cell-imbalance sweep** |
| `shrinkage.R` | the pooled (RE) vs per-gene (fixed-nested) between-sample-variance power gap **as a function of S** |
| `exp7-slope-sweep.R` | the random-slope calibration at S = 12 vs 24 |

Results are the sibling `*_results.rds` / `reply.md`.

## Findings

### 1. The page-5 "bias" is genuine pseudo-replication, not a simulation artifact
`sim_clustered()` draws `u[g,s] ~ N(0, τ²)` i.i.d. and **independent of condition**.
Each gene's responder-vs-non-responder group-mean of `u` differs by chance
(~`N(0, 4τ²/S)`), which the fixed cell-level Wald test picks up at ~C df → **null
type-I ≈ 0.835** (grossly anti-conservative). The reviewer was right that this
mechanism exists; but real patient effects are *not* constrained to zero-mean
within group, so it is a real between-patient signal an honest test must handle.

### 2. Equalising the group means does NOT remove it
Centering `u` to zero mean within each group still left the fixed test at
**≈ 0.61–0.64** (mean|Δ| = 0). Two reasons, the second decisive:
- the log link means equal *log-mean* group means ≠ equal *count*-scale means;
- **the real issue is the error stratum, not the group mean.**

### 3. Nesting alone doesn't fix it — the error stratum does
For the reviewer's `contr.sum`-nested fixed model, on the null:

| test of the response contrast | null type-I |
|---|---|
| cell SE + cell df | 0.835 |
| cell SE + **S−2 df** | 0.835 |
| **between-sample-variance test** | 0.06 |
| random intercept (RE) | 0.033–0.04 |

**The df is a red herring** — the SE itself is cell-level, so S−2 df changes
nothing. Calibration comes only from the **between-sample variance**, which is
mathematically the random-effect treatment. So a *properly-tested* nested fixed
model is calibrated and equivalent to RE for the main effect; the nesting alone
(cell SE) is not.

### 4. Robust to cell imbalance (1×–10×)
Both `nested-between` and RE stay ~0.04 across the whole 1×–10× imbalance sweep;
only the plain cell-level fixed test degrades (0.73 → 0.91). Imbalance is **not** a
differentiator.

### 5. Slopes are inert; the niche interactions need df, not random slopes
On real data τ²_slope ≈ 0 (the niche→expression slope does not vary between
patients), so a random slope is shrunk away and adds nothing (confirmed: prior
studies show no benefit). The niche interactions β multiply a *within*-patient
covariate, so with τ²_slope ≈ 0 they are within-patient contrasts carrying large
effective df (≈ 775 in the toy fit vs S−2 = 8 for the main effect). Testing them at
S−2 is badly **over-conservative** — the real-data symptom — and the fix is the
per-coefficient **Satterthwaite df**, *not* a random slope. (A stress test that
plants τ²_slope > 0 leaves even the random-slope model partly anti-conservative —
a covariance-level effect — but that regime is not where real data sits.)

### 6. Kenward–Roger points the wrong way here
KR *inflates* the covariance (more conservative). Real data is *over*-conservative
(which Satterthwaite relieves), so KR was **held** — it only helps the
anti-conservative regime (τ²_slope > 0 / very small S) that real data isn't in.

### 7. MoM vs MLE for τ²
Pooled Schall/MoM `τ̂² = 0.454` (true 0.49); per-gene REML (glmmPQL, θ=5) mean
0.395 / sd 0.21 — noisier and biased low, at ~100× the per-gene cost. Pooling
across genes compensates for MoM's inefficiency.

### 8. The shrinkage advantage of RE vanishes with S — why `nested` was shelved
The one statistical edge a mixed model has over a calibrated fixed-nested test is
pooling the between-sample variance across genes (shared τ²). Its power advantage:

| S | sparse | mid | dense |
|---|---|---|---|
| 12 | +0.10 | +0.09 | +0.08 |
| 24 | +0.06 | +0.05 | +0.04 |
| **55** | **+0.006** | **+0.005** | **+0.006** |

At realistic cohort sizes (S ≈ 55) the gap is ~0.5% — because each gene's own
between-sample variance then has ~53 df and pooling adds essentially nothing. So a
deterministic fixed-nested engine would match RE with no statistical cost at large
S; RE's shrinkage only pays off at small S (≲ 20).

## The shelved `random="nested"` design (for future pickup)

A deterministic, single-fit calibrated fixed-effects path:

- **Design:** append `contr.sum`-nested sample-intercept columns (S−2, unpenalised)
  so the response contrast is estimable without a penalty.
- **Fit:** one `SpaNorm::fitNB` with `lambda.a = 0` on those columns — no Schall
  loop, no τ², deterministic.
- **Inference (the crux):** a **samples-as-clusters between-sample variance** for
  the response block — mathematically a cluster-robust (CR2) sandwich restricted to
  the response columns, *not* the joint model's cell-level Wald covariance (which
  gives 0.835). Per-coefficient df: main effect → **S−2** (between-sample); niche
  interactions → within-sample residual df (large). Reuse the `@df` named-vector +
  `.ptByCol()` that the Satterthwaite feature already added.
- **Validation oracle:** must reproduce the 0.06 calibration (not 0.835) + the RE
  equivalence at S = 55; the **CR2 small-sample correction is the main risk** (naive
  CR0 is anti-conservative at small S), with `clubSandwich` as a test-only reference
  (the analogue of the `lmerTest` oracle that caught the Satterthwaite `M⁻²Λ` bug).

### Why it was shelved
- At the cohort sizes spiDE targets (S ≈ 55) it is **statistically equivalent** to
  the already-shipped `random="intercept"` + `df.method="satterthwaite"` (§8), so it
  is largely redundant.
- Its only benefit (skip the PQL loop) is a compute saving, not a capability;
  `random="intercept"` is already deterministic at `re.prop = 1`.
- The between-sample-SE computation is subtle and error-prone (the 0.835-vs-0.06
  lesson), so the build carries real risk for a redundant feature.

If picked up: the full design is in this repo's history at
`specs/2026-07-30-nested-fixed-design.md` (commit `446ec7b`), and the calibrated
between-sample test is already prototyped in `imbalance.R`
(`arms_type1`, the `nested-between` arm).
