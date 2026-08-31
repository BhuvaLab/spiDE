# Quasi-likelihood dispersion (edgeR v4 style) — design

- **Date:** 2026-08-31
- **Status:** **Gate 0 passed 2026-09-01 — build Stage A.** H1 (the winsorisation
  clamp) and H3 (the sandwich) are refuted; the inflation sits on the
  dispersion/abundance axis. Two results tighten the design: the v4 adjusted
  deviance and effective df are **hard requirements**, because the raw deviance
  against `n - p` has median 0.267 against Pearson's 0.837 and would inflate `t`
  by 1.94x instead of the current 1.09x; and `zeroish` is 0 for every gene, so
  the df adjustment cannot come from edgeR's legacy zero-count rule. Full result
  in `research/fdr-ordering/FINDINGS.md`.
- **Branch strategy:** `quasi-likelihood-dispersion`, created from `main`. It
  also carries two independent correctness fixes (`.isSelfNiche`, the one-sided
  ACAT input) that are unrelated to QL and can be reviewed separately.

## Background: the defect this is meant to fix

`research/fdr-ordering/FINDINGS.md` (addendum, 2026-08-31) localises spiDE's
loss of FDR control to the **shape of the null p-value tail**, and then
localises the tail to a **per-gene scale error**:

| standardisation of the null `t` | excess at `p < 1e-6` |
|---|---|
| raw | **59.9x** |
| per column (index x niche) | 37.3x |
| **per gene** (across its 136 columns) | **2.6x** |

Per-gene `sd(null t)` spans 0.837–2.292; per-column only 0.908–1.245. 90.8% of
`|t| > 4.89` exceedances come from the top 5% of genes by `sd(t)`. The affected
genes are the **highest-expressed** ones (`cor(sd(t), log median SE) = -0.60`;
worst are HLA-DPA1 2.29, CDV3 2.16, AEBP1 2.07, CST3 2.05, LAPTM5 2.04). It
survives `free` shuffles, so it is not spurious spatial regression.

**The awkward fact this design has to answer.** `.blockedInference()` already
applies a per-gene scalar variance correction — the working Pearson dispersion
at `R/inference.R:539`:

```r
dispb <- rowSums((Yb - mub)^2 / (mub + psib * mub^2)) / disp_df
```

A per-gene scalar is therefore already in the SEs, and the null `t` is still
per-gene over-dispersed by up to 2.3x. So "add a per-gene dispersion" is not
by itself the fix. Whatever is adopted must differ from the Pearson dispersion
in a way that matters, or must stop being a per-gene scalar altogether.

## What edgeR v4 actually does

Read from the installed **edgeR 4.8.0**, `edgeR:::glmQLFit.default` with
`legacy = FALSE` (not from the paper's prose — v4 differs from most summaries
of it):

1. The NB dispersion becomes a **single common value**, estimated by
   `estimateGLMCommonDisp` on the **top-abundance genes only**
   (`top.proportion` from `chooseLowessSpan(ngenes * sqrt(df.residual))`).
   Tagwise NB dispersion is abandoned.
2. First `glmFit` at that common dispersion; `.cxx_compute_ave_qd` returns a
   per-gene **average QL dispersion**; the model is **refit** at
   `dispersion / ave.ql.disp`.
3. `.cxx_compute_adj_vec` returns the **adjusted deviance** and an **adjusted
   residual df** — a per-observation effective df, not `n - p`. This is the
   small-count correction. `s2 = deviance.adj / df.adj`.
4. `squeezeVar(s2, df = df.adj, covariate = AveLogCPM, robust = TRUE)` — EB
   moderation against an abundance trend, giving `s2.post` and `df.prior`.
5. Testing is **F(df1, df.adj + df.prior)** scaled by `s2.post`.

The architectural idea, and the reason it is a candidate here: v4 deliberately
moves gene-specific variability **out** of the NB dispersion — where it enters
the IRLS weights non-linearly — and **into** a multiplicative QL dispersion with
an EB prior and an F reference. The measured defect is a per-gene multiplicative
scale error. The shapes match.

### Where the gain would have to come from

Since a per-gene scalar is already applied, the improvement must come from the
parts of v4 that are *not* a rescale:

- **deviance rather than Pearson residuals.** Pearson `chi^2/df` is badly biased
  for small counts, and CosMx counts are small.
- **the adjusted per-observation df.** spiDE uses `disp_df = n_cells - n_fixed`,
  which counts a full df for cells whose fitted `mu` is near zero. That inflates
  `disp_df`, deflates the dispersion, and shrinks the SE — a concrete candidate
  for the observed under-correction, and it is *not* fixed by any rescale.
- **robust EB moderation**, which stops outlier genes dragging the trend.
- **the F reference**, which accounts for uncertainty in the dispersion. At
  N = 77,454 this is a minor term and is not expected to carry the result.

### The competing hypothesis this design must not assume away

If the residual variance is heterogeneous **within** a gene — larger in the
cells that carry leverage on the niche coefficient — then no per-gene scalar
can fix it, and the answer is a **cluster-robust (sandwich) covariance**
clustered by patient, not QL. Gate 0 decides between the two.

## Gate 0 — the discriminating experiment (do this first)

One fit on one existing shuffle grid (~1 GPU-hour). Record per gene:
`psi_hat`, mean `mu`, Pearson dispersion, deviance-based `s2`, adjusted
`df`, the fraction of counts winsorised, and the spread of per-patient residual
variance. Regress the per-gene `sd(t)` inflation on those.

| outcome | reading | action |
|---|---|---|
| inflation tracks `clamp_frac` / `clamp_shrink` / `lmu_range` | **H1**: `calculateMu()`'s per-gene clamp of fitted log-means at `rowMedian + 4*rowMad` biases the working weights, and is invisible to the Pearson dispersion because that statistic uses the same clamped mu | neither QL nor sandwich — it is plumbing; see the `winsor` note below |
| inflation tracks the Pearson-vs-deviance gap, or the `disp_df` over-count | **H2**: a per-gene scalar estimated the wrong way | build the QL path below |
| inflation tracks the per-patient residual-variance spread | **H3**: within-gene heteroscedasticity | build the sandwich instead; QL will disappoint |
| inflation tracks `psi` at fixed abundance | **H4**: dispersion moderation is wrong at the top of the expression range | fix the dispersion, not the test |

**A `winsor` finding that reshaped this gate.** `winsor` does **not** cap counts
(CLAUDE.md says it does; that is wrong and is corrected separately). It clamps
(a) each coefficient column across genes to `median +/- k*MAD`
(`winsoriseCols()`), and (b) each **gene's** fitted log-mean at
`rowMedian(lmu) + k*rowMad(lmu)`, one-sided, inside `calculateMu()`. None of
spiDE's nine `calculateMu()` call sites passes `winsor`, so all of them use the
default 4 regardless of what the user gave `fitSpiDE()`, and `SpiDEFit` does not
store it — so `fitSpiDE(winsor = Inf)` fits unwinsorised and is then tested
against mu clamped at 4 MADs. That is a defect in its own right, and H1 is the
hypothesis that it is also the cause.

**Nothing below is built until Gate 0 returns.** Cost of getting this wrong is
roughly a GPU-day of revalidation per configuration.

## Design

### Stage A — QL-on-top (recommended first increment)

Keep the existing tagwise NB fit. Add the QL dispersion and the F reference at
the **inference stage only**. No refit, so it is isolated, cheap, and testable
against every artefact that already exists.

1. **`R/inference.R`, inside the existing block loop** (where `mub` is already
   materialised, so no extra genes x cells intermediate):
   - NB unit deviance
     `d_i = 2 * (y log(y/mu) - (y + 1/psi) log((y + 1/psi)/(mu + 1/psi)))`,
     with the `y = 0` limit handled;
   - the per-observation effective df of edgeR v4. **This must be reimplemented**
     — edgeR's is in C (`.cxx_compute_adj_vec`) and a Bioconductor package may
     not call `edgeR:::`. Read the exact adjustment from edgeR's C source; do
     not reconstruct it from the paper;
   - `s2_g = sum(d_adj) / df_adj`, accumulated per gene.
2. **New cross-gene stage.** `limma::squeezeVar(s2, df = df_adj,
   covariate = AveLogCPM, robust = TRUE)` over the length-G vector.
3. **`R/AllClasses.R`.** New `SpiDEFit` slots: `s2`, `s2.post`, `df.prior`,
   `df.residual.adj`, `ql.trend`. All `NULL` on a non-QL fit, as `tau2`/`penalty`
   already are.
4. **Reference distribution.** `@df` gains a **third** shape — it is currently a
   scalar (`between`, or a fixed-effects `Inf`) or one value per tested column
   (`satterthwaite`); QL adds a per-gene value, and per-gene x per-column when
   combined with Satterthwaite. `.dfFor()` and `.ptByCol()` are the two places
   that must learn it; CLAUDE.md's "anything reading `@df` must handle both
   shapes" becomes three.
   Keep the **signed** statistic: the tested contrasts are 1-df, so `F = t^2`
   and a signed `t` on `df.total = df_adj + df.prior` is equivalent. This
   preserves `DirectionNiche` and the up/down split in `.hierarchicalFDR()`
   unchanged.
5. **`R/mixed.R`.** `tau2` is estimated from NB working weights, so it is
   unchanged under Stage A (the fit is untouched) — but the df composition is
   not: `.satterthwaiteDF()` and `df.prior` must be combined. edgeR/limma add
   them; adopt that and extend the lmerTest oracle in
   `tests/testthat/test-satterthwaite.R` to cover it.
6. **`R/spiGSEA.R`.** The t-to-z conversion inherits the new df; no structural
   change, but the calibration benchmark must be re-run.

### Stage B — full v4 (only if Stage A under-delivers)

Common NB dispersion + refit at `dispersion / ave.ql.disp`. Expressible through
`SpaNorm::fitNB(psi = ...)` with **no SpaNorm change**, since `psi` is already
an argument.

**Hazard, load-bearing.** `fitNB`'s `winsor` caps counts against the *current
fitted mu*, so supplying `psi` changes the IRLS weights, hence `mu`, hence which
cells are capped: the estimating and supplied-psi paths maximise **different**
objectives (`research/notes/fitnb-offset-psi-disagreement.R`). Consequences:
- Stage B cannot be validated by comparing log-likelihoods — likelihood cannot
  rank two different objectives. Score type-I and power only.
- Stage B changes `tau2` (the PQL loop re-fits under the new weights) and
  therefore the Satterthwaite df. It is not a drop-in.

### API

`fitSpiDE()`/`spiDE()` gain `dispersion = c("nb", "ql")`, defaulting to `"nb"`
until the validation below passes. `testSpiDE()` needs no new argument — the
choice is carried on the fit, as `mode` and `random` already are.

## Computational feasibility

Real dimensions from the run logs (`research/plasmode/logs/gpu_rintercept_*`):
**G = 13,348 genes, N = 77,454 cells, S = 55 patients, p ~ 360** (136
`(index, niche)` columns pre-fix, 132 post-`.isSelfNiche`-fix).

| work | scaling | flops | vs. existing |
|---|---|---|---|
| existing per-gene Gram `X'W_gX` | `O(G N p^2)` | 1.3e14 | **100%** (16 min GPU, 69–100 min CPU) |
| **QL deviance + adjusted df** | `O(G N)` | ~1e10 | **0.008%** |
| EB moderation (`squeezeVar` over 13,348 values) | `O(G)` | ~0 | milliseconds |
| cluster-robust sandwich (the alternative) | `O(G N p) + O(G S p^2)` | 4.7e11 | **0.35%** |

**Nothing in the QL switch is expensive.** It is one extra `O(G N)` reduction
over a `mub` that is already computed. The same is true of the sandwich, which
is cheaper than the bread already being formed. The cost is entirely in (a) a
refit, if Stage B, and (b) revalidation.

### The one architectural consequence

`s2` and the abundance trend are needed across **all** genes before any gene's
SE is final, which bends the "fit whole, infer blocked" invariant. Resolve it as
**two passes**, not by abandoning blocking:

- **pass 1** accumulates three numbers per gene (`s2_g`, `df_g`,
  `AveLogCPM_g`) — this is `O(G N)` and is fused into the existing `mu` loop, so
  it is nearly free;
- `squeezeVar` on the G-vector;
- **pass 2** forms the statistics — this is the expensive Gram pass, unchanged.

Alternatively cache `s2`/`df_adj` on the `SpiDEFit` during `fitSpiDE()` so
`testSpiDE()` stays single-pass. Prefer this: it keeps the published two-stage
user API intact and puts the extra pass where a full-data sweep already happens.

### Sampling policy

The project already has the correct pattern in `.fitNBmixed()`: subsample for
the **shared/nuisance** quantity, then take a **final full-cell pass** for
anything per-gene that inference consumes. Applied here:

**Legitimate to subsample**
- the abundance trend, `df.prior`, and the `chooseLowessSpan` span — smooth
  functions of G points, negligibly changed by a stratified cell subsample;
- the common NB dispersion under Stage B — edgeR itself estimates it from a
  *gene* subset (the top-abundance genes), which is the cheaper axis anyway.

**Not legitimate to subsample**
- `s2_g` itself. It multiplies the SE directly, so subsampling injects an
  `O(1/sqrt(n_sub))` error into every `t`: at 5% of 77,454 cells that is ~2% per
  gene, the same order as the effect being corrected. Compute it on all cells —
  it costs 0.008% of the run.

**Discipline**
- stratify per (cell type x sample), as `re.prop` does, or thin index types
  vanish from the subsample entirely;
- `fitNB` sets **no seed internally**, and at 600 cells two identical calls
  differ by `max|d psi| = 0.245`. Every sampled quantity must be seeded and the
  seed recorded in the output.

### GPU policy

All new work is elementwise plus row reductions, and the fp64 primitives already
exist: `SpaNorm::rowSums_gpu`, `add_vec_mat_gpu`, `mult_vec_mat_gpu`,
`dnbinom_gpu`. Fuse the deviance accumulator into the existing `full_cov` branch
so no extra genes x cells tensor is allocated.

- SpaNorm's GPU path is **fp64**, so **h100 only** — a100 fails fp64 `digamma`
  inside NB dispersion, l40s runs fp64 at ~1:64.
- The NB **deviance needs only `log`, no `digamma`**. Under Stage B (common
  dispersion) the a100 blocker would disappear for the fit as well. Worth
  measuring, not worth assuming.
- **Standing blocker:** `testSpiDE()` currently dies in a torch CUDA kernel JIT
  (nvrtc) on this cluster, so the inference stage is CPU-only in practice today.
  Build CPU-first; gate the GPU inference path behind that fix. Note
  `research/plasmode/gpu_smoke.R` passing does **not** imply the blocked
  inference path runs.

## Validation plan

Score on the existing niche-shuffle nulls, **never on likelihood** (see the
Stage B hazard). Acceptance criteria, with current values:

| # | criterion | now | target |
|---|---|---|---|
| 1 | per-gene `sd(null t)` flat across expression deciles | 1.135 (top) → 0.946 (bottom) | flat within noise |
| 2 | fraction of genes with `sd(t) > 1.3` | 0.0091 | ~0 |
| 3 | tail ratio `P(p < 1e-6)/1e-6` | 59.9x | ~1 |
| 4 | complete-null false calls, flat BH at alpha .05 | 224 | ~0 |
| 5 | power on the injected-signal grids | see `research/fdr-ordering/summary/` | not worse |

Criteria 1–3 are precise with 2–3 grids; only 4–5 need all 18. **Budget:** one
full-cohort fit is 16 min (GPU) / 69–100 min (CPU), so 18 shuffles plus the real
run is roughly a GPU-day per configuration. That, not the QL arithmetic, is what
to schedule. Pin `SPIDE_PKG` to a frozen snapshot
(`.claude/skills/run-benchmark-arm/scripts/freeze_snapshot.sh`) — a sweep whose
tasks load the live working tree is not one experiment.

Benchmark arms land as extra **rows** in the one canonical table per scenario
under `research/reports/benchmarks/tables/`, never a parallel file.

## Explicitly out of scope

- Changing the FDR cascade, its ordering, or its alpha-spending. Those are
  separate findings in `research/fdr-ordering/FINDINGS.md` and are independent
  of the variance estimator.
- The gene filter (dropping genes with shuffle `sd(t) > 1.3`). **Evaluated
  2026-09-01 and it does not make QL unnecessary**: defined honestly on
  independent grids it is a **1.5x** reduction in null false calls, not the 6x
  first reported (that figure was circular and is withdrawn). It is worth
  applying — it costs no power, worst measured TPR loss 0.011 — but
  `P(>=1 false call)` stays 1 with it on.
- `twoStageSpiDE()`, which shares none of this machinery.

## Files touched (Stage A)

```
R/inference.R      deviance + adjusted df in the block loop; squeezeVar stage
R/AllClasses.R     s2, s2.post, df.prior, df.residual.adj, ql.trend slots
R/fitSpiDE.R       dispersion = c("nb", "ql") argument, cached s2 on the fit
R/mixed.R          combine df.prior with the Satterthwaite df
R/spiGSEA.R        df propagation only
tests/testthat/    test-inference.R, test-satterthwaite.R, a new test-ql.R
vignettes/spiDE-model.Rmd
```
