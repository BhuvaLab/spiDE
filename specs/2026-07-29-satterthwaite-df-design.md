# Per-coefficient Satterthwaite degrees of freedom for the mixed-effects fit

**Branch:** `satterthwaite-df` · **Date:** 2026-07-29 · **Status:** design (awaiting review)

## Problem

The mixed-effects fit (`fitSpiDE(..., random != "none")`) references every
Response / ResponseNiche Wald statistic against a **single global**
`df = max(S - 2, 1)` (`SpiDEFit@df`, set in `.fitNBmixed()`,
[R/fitSpiDE.R](../R/fitSpiDE.R)). On the 55-sample YTMA cohort that is
`df = 53`, applied uniformly even though the tested coefficients span an
effective sample size of roughly `n_eff = 4 → 2414`.

A single df cannot be right for both ends:

- The **Response main effect** is a pure between-sample contrast (the condition
  is a patient-level label, enforced by `checkSample()`). Its genuine effective
  df *is* ≈ `S - 2`, which is why it calibrates almost perfectly — observed
  type-I ≈ 0.99× nominal.
- The **ResponseNiche** terms mix in *within-sample* niche-density variation, so
  their effective df is far larger. Forcing them to `S - 2 = 53` over-fattens the
  t reference → over-conservative tests → the SE/power inflation seen on realistic
  data. A rare index cell type (`n_eff = 4`) has the opposite failure.

This is a property of spiDE's inference machinery, not of the data. The fix is a
**per-coefficient** reference df — Satterthwaite — standard in `lmerTest` /
`pbkrtest` for exactly this small-sample mixed-model regime.

## Scope decision

Chosen (agreed with the author):

- **Satterthwaite first.** Implement per-coefficient Satterthwaite df. Assess it on
  the simulation harness. Add the heavier **Kenward–Roger** covariance
  bias-correction *only if* the assessment shows residual inflation the df alone
  does not remove (see "KR decision gate"). KR is otherwise out of scope for this
  branch.
- **Assess on the simulation harness** — the research-branch `null` scenario
  (known ground-truth type-I error), not (for this branch) the real YTMA cohort.

## Statistical core

### The Satterthwaite df

Write the penalised covariance of the tested (fixed) coefficients as
`Ĉ = φ · M⁻¹`, with

```
M   = C' W̄ C + Λ            full penalised information (C = [fixed | random])
Λ   = diag(0 on fixed cols, 1/τ²_m on the random cols of group m)
W̄   = diag(w̄_c)             representative (gene-averaged) working weights
```

For a tested coefficient `j`, Satterthwaite df is `2·(se²_j)² / Var(se²_j)`, with
`Var(se²_j)` a delta-method sum over the variance parameters
`θ = (φ, τ²_int, τ²_slope)`:

```
df_j = 2·v_jj²  /  [ 2·v_jj²/df_resid  +  Σ_{m,m'} g_jm · g_jm' · Â[m,m'] ]

  v_jj = [M⁻¹]_jj                                    fixed-block diagonal, col j
  g_jm = (1/τ²_m²) · Σ_{k ∈ group m} ([M⁻¹]_jk)²     sensitivity of v_jj to τ²_m
  Â    = Cov(τ̂²)                                     variance-component covariance
  df_resid                                           residual df for φ̂
```

Derivation of `g_jm`: `Λ` depends on `τ²_m` only through the group-`m` diagonal
`1/τ²_m`, so `∂M⁻¹/∂τ²_m = (1/τ²_m²) M⁻¹ E_m M⁻¹` (`E_m` = diagonal indicator of
group `m`), whose `(j,j)` entry is `(1/τ²_m²) Σ_{k∈m}([M⁻¹]_jk)²`. The φ term uses
`Var(φ̂) ≈ 2φ²/df_resid` (scaled-chi-square). φ and τ² are treated as
asymptotically orthogonal (no cross term) — the standard Satterthwaite
simplification, to be checked numerically.

### Key insight — `df_j` is invariant to the per-gene dispersion `φ_g`

Both `se²_j = φ·v_jj` and every gradient term scale homogeneously in `φ`, so `φ`
cancels exactly in the ratio (verify: numerator `2φ²v_jj²`; φ-part of denom
`v_jj²·2φ²/df_resid`; τ²-part `φ²·g·Â·g` — every term carries `φ²`). Therefore:

- The **only** per-gene dependence left is the working-weight *pattern* `w̄`
  (through `M⁻¹`).
- The `n_eff = 4 → 2414` spread is a **per-coefficient (structural)** spread — it
  is about which error stratum a coefficient lives in, not about which gene.

**Consequence (the efficiency win):** a single **shared per-coefficient df
vector** (length `k` = number of tested columns), computed once per bandwidth on
the representative weights, captures essentially all of the adaptation. This is
the same "single shared population quantity" logic already used for `τ²`. No
per-gene df matrix; nothing added to the batched per-gene inference loop.

*Deliberately deferred:* genuine per-gene df (a `genes × k` object built from the
per-gene `M_g⁻¹` already inverted in `.waldCauchyBlock()`). The φ-invariance
argument says it buys little; if the assessment ever shows it matters, the slot
contract below (`@df` may be scalar / vector / matrix) already admits it without
reshaping the consumers.

### Response ≈ S − 2 falls out — the regression guarantee

For the between-sample Response contrast, `v_jj` and the sensitivity `g` are
dominated by the `τ²_int` direction, so the `g·Â·g` term dominates the
denominator and `df_j` computes back to ≈ `S - 2`. The 0.99× calibration is
reproduced *by the same formula*, not preserved by a special case. This is the
core regression test: **the Response column's Satterthwaite df must stay ≈ S − 2
and its null type-I ≈ 0.05.**

### The one genuinely new quantity — `Â = Cov(τ̂²)`

The asymptotic covariance of the variance-component estimator, via the
**reduced-form REML expected information** at representative weights:
`I_{m,m'} = ½·tr(P V_m P V_{m'})` with `V_m = Z_m Z_m'` the group-`m`
random-effect outer product and `P` the REML projection — reduced to traces over
`p × p` `M⁻¹` sub-blocks and `Z'W̄Z`, `Z'W̄C` cross terms, **never forming any
`n × n` matrix**. Then `Â = I⁻¹` (the τ² block). This is the analog of what
`pbkrtest` computes, adapted to the PQL working model.

- **Primary validation** — a Gaussian-response analog (identity link, known
  weights ⇒ the PQL working model *is* an exact LMM): assert our df matches
  `lmerTest`'s Satterthwaite df on a small dataset (guarded by
  `requireNamespace("lmerTest")`; add to `Suggests`). This pins the df algebra
  independently of the simulation.
- **Fallback** — if the 2-component `Â` (intercept + slope) is unstable for the
  slope term at very small `S`, fall back to intercept-only τ² sensitivity for
  the slope model, or a moment approximation `Cov(τ̂²_m) ≈ 2·τ⁴_m / df_m`. The
  chosen fallback is recorded and surfaced in `verbose` output.

## Where it lives & data flow

**Computed once per bandwidth, in the fit stage.** `.fitNBmixed()` already builds
the representative penalised inverse (`wbar` → `info` → `minv`,
[R/fitSpiDE.R](../R/fitSpiDE.R) lines ~117–120). A new helper
`.satterthwaiteDF()` reuses that machinery — recomputed once on the final
all-cell fit at the converged penalty — to produce the length-`k` df vector.
Cost: **one `p × p` inverse per bandwidth, never batched over genes** → trivially
within any GPU/CPU memory budget. The expensive per-gene inference path is
untouched.

**Representation.** `SpiDEFit@df` generalises from a scalar to a **named numeric
vector aligned to the tested (Response + ResponseNiche) columns** — i.e. the same
columns, in the same order, as `@t_stat` / `@se`. Contract:

- `NULL` → fixed-effects fit, normal reference (unchanged).
- length-1 (scalar) → `df.method = "between"`, the current global `S - 2`,
  recycled across all tested columns (back-compatible).
- length-`k` vector → `df.method = "satterthwaite"`.
- (`genes × k` matrix admissible by the same consumers, reserved for a future
  per-gene df; not produced by this branch.)

**Consumers** — both apply a `t` reference from `@df` to a `genes × (tested cols)`
statistic matrix. Route both through one small helper:

```r
.ptByCol(tmat, df, lower.tail)   # df: NULL | scalar | length-ncol(tmat) | matrix
```

which recycles a scalar, broadcasts a per-column vector down the columns, or uses
a matrix elementwise, and returns `pnorm` when `df` is `NULL`.

1. `.blockedInference()` / `.waldCauchyBlock()` / `.waldBrownGene()`
   ([R/inference.R](../R/inference.R)) — the `ptail()` closures become
   `.ptByCol()` calls. `df` there covers all `k` tested columns (Response +
   ResponseNiche), matching `@t_stat`'s column order.
2. `.nicheRecords()` ([R/fdr.R](../R/fdr.R) lines ~84–88) — subsets the
   ResponseNiche columns (`rn`) of both `t_stat` **and** `@df` identically before
   calling `.ptByCol()`.

**The per-gene inference loop and all GPU/tensor code stay unchanged** — `@df` is
precomputed and only indexed, so no new per-gene tensor work, no new memory in
the hot path.

**API.** New argument `df.method = c("between", "satterthwaite")` on `fitSpiDE()`
and `spiDE()` (forwarded through `.fitOneBandwidth()` → `.fitNBmixed()`). Ignored
for `random = "none"`. **Default `"between"`** for this branch, so the behaviour
change is opt-in during assessment; flipping the default to `"satterthwaite"` for
mixed fits is a follow-up the assessment justifies.

## Vignette / report updates

The statistical model report is the home for the df statistics; the mixed
benchmark report gains the calibration evidence.

### `vignettes/spiDE-model.Rmd` (the statistical model report)

- **§ "Corrected inference", point 3 (reference distribution).** Currently claims
  the reference is `t` with `S - 2` df "regardless of how many niche coefficients
  are simultaneously tested." Rewrite: `S - 2` is the correct df **for the
  between-sample Response contrast**, but the ResponseNiche contrasts are informed
  by within-sample niche variation and warrant a **larger, per-coefficient** df.
  Add a new subsection **"Per-coefficient degrees of freedom (Satterthwaite)"**
  presenting: the df formula, the φ-invariance argument (why one shared
  per-coefficient vector suffices — same logic as the shared `τ²`), the
  `Â = Cov(τ̂²)` REML-information quantity, and the derivation that the Response
  column returns ≈ `S - 2` (so the between-patient df is the *special case for the
  main effect*, not a global constant). A small runnable chunk fits the toy
  clustered data both ways and shows: Response df ≈ `S - 2` under both; the
  ResponseNiche columns get larger df under `"satterthwaite"`.
- **New subsection "The testing framework".** Document how calibration is judged:
  the `null` simulation scenario (no planted effect, real between-patient
  clustering), the type-I metric at α = 0.05, the Response 0.99× target as a
  regression anchor, and the requirement that ResponseNiche type-I rises *toward*
  0.05 from below without overshooting. This satisfies "detailing the testing
  framework."
- **§ "Interpreting a mixed `SpiDEFit`", `@df` row.** Update from "the
  between-patient reference df, `S - 2`" to: a per-coefficient vector under
  `"satterthwaite"` (Response ≈ `S - 2`, niche terms larger), or the scalar
  `S - 2` under `"between"`; how to read a large vs small entry.
- **§ "Usage and caveats", "Power needs samples and effect size".** Note that
  per-coefficient df restores legitimate power for niche terms that carry
  within-sample information, without touching the honest between-sample df of the
  main effect.

### `vignettes/spiDE-mixed-benchmark.Rmd`  (calibration figure — in scope for this branch)

Add a **calibration** section (figure + table, read from a shipped benchmark CSV
via the existing `bench()` helper, exactly like the timing/accuracy figures):
null type-I error for `fixed` / `intercept` / `slope` × `df.method` across the
`S` grid, for the Response effect and the ResponseNiche effect separately — the
visual evidence that Satterthwaite relieves ResponseNiche over-conservatism while
leaving Response at ≈ 0.99×.

**Artifact.** The figure reads a shipped CSV
`inst/extdata/benchmark/null_calibration.csv` (columns roughly `S`, `layout`,
`method`, `df.method`, `effect` ∈ {`response`, `responseniche`}, `type1`). This
branch **produces and ships that CSV** from the local small-grid calibration run
(assessment step 1–2 below) so the vignette renders self-contained, matching how
`pql_timing.csv` / `timing.csv` are already shipped. Regenerating it at full grid
on HPC (the research-branch `null` scenario) is the same copy-into-`inst` step the
other benchmark CSVs use.

## Assessment plan (simulation harness)

Primary, runnable in-repo at small scale (a few `S` values, `sd_patient > 0`,
`null` scenario) — the deliverable that answers "does it work":

1. **Response regression** — computed df ≈ `S - 2` *and* type-I ≈ 0.05 under both
   `df.method` values (the 0.99× anchor is not disturbed).
2. **ResponseNiche fix** — type-I moves *up toward* 0.05 from below under
   `"satterthwaite"` (over-conservatism relieved) **without overshooting** into
   anti-conservatism (type-I ≤ ~0.05).
3. **HPC-ready path** — reuse the research-branch `null` scenario
   (`research/R/scenarios.R`, already wired to `fit@df`) to run the full grid
   (`S ∈ {4,10,16,24,30}`, both layouts) for the publication numbers. The harness
   metric extraction (`.ab_p`, and the `df <- if (is.null(fit@df)) Inf else
   fit@df` line in `scenarios.R`) must be updated to read the **per-coefficient**
   df for the tested column rather than a scalar.

### KR decision gate

If, after Satterthwaite, the ResponseNiche null type-I is **> ~0.05**
(anti-conservative) — i.e. the residual miscalibration lives in the covariance,
not the df — then add the Kenward–Roger covariance bias-correction. The design
leaves a clean seam for this: `df.method` gains `"kenward-roger"`, and a
covariance-adjustment hook sits at the point where `.blockedInference()` forms
`(C'W̄C + Λ)⁻¹`. No rework of the Satterthwaite plumbing is required.

## Tests

- **`lmerTest` cross-check** (new) — Satterthwaite df on a Gaussian-response
  working analog matches `lmerTest` within tolerance (guarded).
- **Response df ≈ S − 2 regression** (new) — on `.toyClustered()`, the Response
  column's Satterthwaite df is within tolerance of `S - 2`.
- **`"between"` path unchanged** (new) — `df.method = "between"` produces
  byte-identical `@df`, `p.combined.*`, and `results()` to the current code.
- **Per-column threading** (new) — `@df` vector length/names align with `@t_stat`
  columns; `.nicheRecords()` uses the matching per-niche df; `results()` p-values
  reflect the per-coefficient df.
- **Back-compat of the slot** (extend `test-mixedEffects.R`) — `NULL` for
  `random = "none"`; correct type/shape for mixed fits.
- **`random = "none"` untouched** — fixed-effects path bit-identical.

## Out of scope

- Kenward–Roger (unless the decision gate trips).
- Per-gene (`genes × k`) df.
- Changing the default `df.method`.
- Any change to the `random = "none"` fixed-effects path.

## Risks

- **`Â = Cov(τ̂²)` accuracy** — the main technical risk; mitigated by the
  `lmerTest` cross-check and the moment fallback.
- **Slope-component instability at small `S`** — `τ²_slope` is collinear with the
  tested fixed effect; its variance is noisy. Fallback to intercept-only
  sensitivity for the slope model.
- **Overshoot into anti-conservatism** — larger df could reintroduce the very
  anti-conservatism the mixed correction fixed. Directly measured by the
  assessment; the KR gate is the remedy.
```
