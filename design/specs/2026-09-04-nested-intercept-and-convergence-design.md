# Nested (sample × cell type) intercept and per-gene convergence: design

**Status:** approved in chat 2026-09-04; implementation plan follows.
**Evidence:** `design/specs/2026-09-03-sample-celltype-intercept.md` (the defect
and its validation), `research/fdr-ordering/REPORT.md` §5c–§5d,
`research/fdr-ordering/FINDINGS.md` (2026-09-03 entry).

## Why

Two independent defects were found in the fit that every spiDE inference is
built on.

1. **The niche design has no (sample × cell type) intercept.** With one
   ridge-penalised intercept per sample shared across cell types, the tested
   `CellType:condition:niche` slopes carry a between-sample composition effect
   (samples whose type-*k* cells sit in denser type-*n* surroundings also have
   a different mean expression in type *k*). The shuffle null, which permutes
   within (sample, cell type), preserves it. It is a bias of the estimand, not
   a variance, which is why twelve variance-side candidates failed. Making the
   niche slopes within-(sample, cell type) — the research equivalent is
   centring the niche-dependent columns — calibrates the null in every
   expression band on all twelve raw grids.
2. **`SpaNorm::fitNB` does not reach the per-gene optimum for bright genes.**
   Its multi-gene IRLS shares one cell-weight vector across genes, decides
   convergence on the aggregate log-likelihood and clamps coefficients across
   genes. The top 5% of genes by expression sit 1–4 production SEs from their
   own penalised-NB optimum, and edgeR's dispersion estimated at that point is
   1.6× too large. Converging each gene (damped Newton from a sane start,
   profile-ML ψ at the converged means) is exactly calibrated on ordinary
   genes and is the honest baseline for the bright ones.

Both go into the package. They are phased so that neither blocks the other.

## Scope

In scope:

- A `SampleCellTypeInt` random-intercept block in `.buildRandomEffects()`,
  controlled by `re.celltype = TRUE` on `fitSpiDE()` / `spiDE()`.
- A per-gene convergence stage `.polishFit()` after the final `fitNB`, in both
  the fixed-effects and mixed paths, controlled by `converge = TRUE`,
  `converge.maxit`, `converge.tol`.
- A toy fixture with a planted between-sample composition effect and zero
  within-sample slope.
- Tests, vignette and CLAUDE.md updates.

Out of scope (a second spec):

- Absorbing the indicator block inside `.blockedInference()`'s batched grams,
  or reusing the polish stage's converged standard errors there. With ~660
  extra dense columns the per-gene inference gram is ~8× slower on the real
  cohort (hours rather than an hour on CPU). Accepted for now.
- The simulation re-run and the real-cohort calibration check (research
  runs, after the package change).

## 1. The nested intercept

### Design block

`.buildRandomEffects(sample_vec, slope_base, random, cell_type_vec = NULL)`
gains a fourth argument. When it is non-`NULL`:

```r
grp <- interaction(factor(sample_vec), factor(cell_type_vec), drop = TRUE)
Zct <- stats::model.matrix(~ 0 + grp)                # one 0/1 column per non-empty (sample, cell type)
colnames(Zct) <- paste0("SampleCellType", levels(grp))
re_group <- c(re_group, rep("SampleCellTypeInt", ncol(Zct)))
```

appended after the per-sample block (and after the slope block under
`random = "slope"`). The per-sample block is **kept**: it is the patient-level
counterpart of the condition contrast, and the Satterthwaite df of
`ResponseCellType` depends on its variance component. The new block is nested
inside it.

`.buildNicheDesign()` passes `cd[[cell_type]]` through when
`re.celltype && random != "none"`. Column tagging is unchanged: the new columns
are `type = "Random"` in `coefmap`, `index = niche = NA`, so `.testedCols()`,
`.nicheTestCols()` and the inference's `disp_df` (which counts non-`Random`
columns) are untouched.

### Variance components and df

`.fitNBmixed()` estimates one τ² per `re_group`, so `SampleCellTypeInt` gets
its own component with no change to the Schall loop. `.varParamCov()` and
`.satterthwaiteDF()` iterate over groups and pick it up. The `"between"` df
formulas are unchanged (they read `SampleInt` / `SampleSlope` only).

`tau2.init` stays a scalar applied to every group.

### API

`fitSpiDE(..., re.celltype = TRUE)` and `spiDE(..., re.celltype = TRUE)`.
Ignored (with no warning) when `random = "none"`. `re.celltype = FALSE`
reproduces every earlier design column for column.

`checkSample()` documentation: with `re.celltype = TRUE`, covariates constant
within (sample, cell type) are absorbed by the new block as well; the
existing rejection of sample-constant covariates stands.

### Cost

Up to S × K extra columns (55 × 12 = 660 on the cohort, fewer where a type is
absent from a sample). `fitNB`'s own IRLS forms one shared gram per iteration,
so it slows by the column ratio squared (~8×) per iteration but remains
minutes. The polish stage absorbs the block (below). Inference does not (out
of scope).

## 2. The convergence stage

### Placement

`.fitOneBandwidth()` calls `.polishFit()` after the fit returns, in both
paths:

```r
if (converge) fit <- .polishFit(Y, W, fit, penalty, re_group, BPPARAM, backend,
                                maxit = converge.maxit, tol = converge.tol, verbose)
```

where `penalty` is the fit's per-column ridge (`fit$penalty` for the mixed
path; `lambda.a` recycled for the fixed path). The stage replaces `alpha`,
`psi`; `loglik` is recomputed afterwards as today. It runs once per
bandwidth.

### Per-gene algorithm

For gene *g* with counts `y`, design `W`, penalty vector `pen`:

1. **Start** at `fitNB`'s `alpha_g`. If the fitted range is degenerate
   (`min(log mu) < -10`) or the first line search fails, restart from the
   sane start: `CellType_k` intercepts at `log(mean(y[cells of k]) + 1e-3)`,
   every other coefficient 0.
2. **Damped Newton** on `ell(a) = sum dnbinom(y, mu = exp(W a), size = 1/psi,
   log = TRUE) - 0.5 * sum(pen * a^2)`:
   score `s = W'((y - mu)/(1 + psi mu)) - pen * a`, information
   `I = W' diag(mu/(1 + psi mu)) W + diag(pen)`, step `d = I^{-1} s`,
   halve the step until `ell` does not decrease (floor 1e-6), stop when the
   gain is below `tol * |ell|` or at `maxit`. The information is **frozen for
   up to three consecutive steps** (only the score is recomputed) and rebuilt
   when a line search halves more than twice or every fourth step.
3. **ψ** by profile ML at the converged `mu`: `optimize` on `log psi` over
   `[log 1e-3, log 1e3]`.
4. **Re-polish** at the new ψ (`maxit = 20`), re-estimate ψ once more, final
   short polish.
5. Return `alpha_g`, `psi_g`, iterations, whether the sane start was used,
   and whether the cap was hit.

Counts are used **unwinsorised** (as validated); `fitNB`'s `winsor` applies
only to its own fit.

Failure handling matches the rest of the package: a singular information
matrix for one gene leaves that gene's `fitNB` values in place with a
per-gene flag (the gene is not dropped; inference already guards its own
inversions). A gene whose polish hits `maxit` keeps the polished values and
is flagged.

### The absorbed indicator block

Inside `.polishFit()` the design is split once into `X` (every column whose
`re_group` is not `SampleCellTypeInt`, dense) and a group index `gi` for the
`SampleCellTypeInt` columns (indicators by construction). With per-cell
weights `w`:

- `A = X' diag(w) X + diag(pen_x)` — the dense gram, as today.
- `B = X' diag(w) Z = t(rowsum(X * w, gi))` — p_x × G, one pass over the cells.
- `C = diag(rowsum(w, gi) + pen_z)` — diagonal.

The Newton system `[A B; B' C] [d_x; d_z] = [s_x; s_z]` is solved by the Schur
complement `S = A - B C^{-1} B'`: `d_x = S^{-1}(s_x - B C^{-1} s_z)`,
`d_z = C^{-1}(s_z - B' d_x)`. `S^{-1}` is also the X-block of the full
penalised covariance, which is what inference needs for the tested columns.
Cost per information rebuild: one p_x-column gram plus two `rowsum` passes,
independent of G.

When there is no `SampleCellTypeInt` group (`re.celltype = FALSE` or
`random = "none"`) the helper is the plain dense Newton.

### Dispatch and backend

Genes are chunked with `.chunkGenes()` and dispatched with
`BiocParallel::bplapply(BPPARAM)`, like `.blockedInference()`. The gram is
formed with a one-argument `crossprod(X * sqrt(w))` on CPU. A GPU path is not
part of this spec: the per-gene loop is sequential in its iterations and the
CPU cost is a few hours for the genome, parallel over genes.

### ψ: a documented departure

`fitNB` moderates dispersion across genes (`edgeR::estimateDisp(robust =
TRUE)`). The polish replaces it with a per-gene profile-ML ψ at the converged
means. With tens of thousands of cells per gene this is well determined, and it
is what was validated on the cohort; it is a departure from `fitNB`'s
moderation and is documented as such in the vignette and the argument help.

### API

`fitSpiDE(..., converge = TRUE, converge.maxit = 50L, converge.tol = 1e-8)`
and the same on `spiDE()`. `converge = FALSE` reproduces the previous
behaviour.

### `SpiDEFit`

One new slot `polish` (`NULL` when `converge = FALSE`, otherwise a
`data.frame` with one row per gene: `iterations`, `psi_fitnb`, `restarted`,
`capped`, `singular`). `show()` reports how many genes were polished, restarted
and capped. `validSpiDEFit()` checks the row count when non-`NULL`.

## 3. Toy fixture: a planted composition effect

`.toySPE()` gains an argument `composition = 0`. When non-zero, for index type
`A` and niche type `B`: each sample draws a mean B density offset (the
`B` cells are placed so that samples differ in how B-rich the environment of
their `A` cells is), and the `A` cells' baseline expression of a second gene
`G2` is shifted per sample in proportion to that sample's mean B density,
with **no** within-sample dependence on the local B density. Under the old
design `G2` is called in `A` × `B`; under the nested design it is not. The
existing `G1` within-sample effect is unchanged and must still be recovered.

The pre-baked `data(toySpiDE)` is not changed.

## A checked non-issue: the dispersion degrees of freedom

`.blockedInference()` computes the Pearson working dispersion with
`disp_df = ncells - #(non-Random columns)`, so the nested block -- tagged
`Random` -- does not reduce it. That is deliberate (penalised columns
contribute little effective df) but worth quantifying, since 660 extra columns
is not obviously "little".

Measured on the toy (240 cells, 43 columns, `random = "intercept"`,
`re.celltype = TRUE`): the `SampleInt` block consumes 2.26 effective df of 6
nominal (38%) and `SampleCellTypeInt` 5.74 of 18 (32%), against a `disp_df` of
221. The honest residual df is therefore ~213, so the dispersion is ~3.7% too
small and `t` ~1.8% too large -- on a fixture with 13 cells per nested group.
Scaled to the real cohort (77,454 cells, 660 nested columns at ~32% effective),
the gap is ~210 df against a `disp_df` of 77,164: **0.3%**.

Two reasons not to "fix" it. The effect is third-order at real cohort sizes,
and the validated configuration -- the centred converged null that measured
`sd(null t)` at 0.96-0.99 -- used exactly this convention, so changing it would
invalidate the number the fix is justified by. Revisit only if a future design
puts a large penalised block on few cells.

## 4. Tests (written before the code)

`tests/testthat/test-design.R`

- `.buildNicheDesign(..., random = "intercept", re.celltype = TRUE)` has one
  `SampleCellType*` column per non-empty (sample, cell type), all tagged
  `Random` with `re_group == "SampleCellTypeInt"`, after the `Sample*` block.
- `re.celltype = FALSE` gives a design identical (`identical()`) to the
  current one.
- Under `random = "slope"` the block follows the slope block.

`tests/testthat/test-mixedEffects.R`

- `fitSpiDE(random = "intercept")` returns `tau2` with names
  `c("SampleInt", "SampleCellTypeInt")`; `random = "slope"` adds
  `SampleSlope`; `df` is still a named per-tested-column vector.
- `re.celltype = FALSE` reproduces the previous `tau2` names.

`tests/testthat/test-polish.R` (new)

- For every toy gene the polished penalised log-likelihood is ≥ `fitNB`'s
  (tolerance 1e-6 relative).
- The penalised score at the polished point has `max(abs(s)) < 1e-4 * max(abs(s at fitNB))`
  or an absolute `1e-6`, whichever is larger.
- A gene fitted alone by `fitNB` to convergence (a single-gene `fitNB` on the
  toy) is changed by less than `1e-3` in every coefficient.
- `psi` is finite and positive for every gene; `polish$iterations >= 1`.
- Block dispatch (`block.size = 3`) and serial give identical `alpha`, `psi`.
- The Schur-complement solve equals the dense solve on a design with an
  indicator block (`all.equal`, tolerance 1e-8), and the X-block of the
  covariance matches `solve(I)[x, x]`.
- The sane-start restart is exercised: a start with `min(log mu) < -10`
  triggers `restarted = TRUE` and converges.
- `converge = FALSE` leaves `alpha`, `psi` identical to `fitNB`'s and
  `polish` `NULL`.

`tests/testthat/test-spiDE-e2e.R`

- The planted `G1`/`A`/`B` within-sample effect is still the strongest niche
  association for `G1` in `A` with `re.celltype = TRUE, converge = TRUE`.
- New: with `composition = 2`, `G2` in `A` × `B` has `fdr.niche < 0.05` under
  `re.celltype = FALSE` and `> 0.2` under `re.celltype = TRUE`, at the default
  bandwidth set (seeded).

`tests/testthat/test-satterthwaite.R`: the lmerTest oracle check still passes
with `re.celltype = FALSE` (its reference model has no nested term) and a
nested-term variant is added with `lme4` `(1 | sample) + (1 | sample:cell_type)`.

All existing suites pass unchanged where their arguments do not touch the new
defaults; where a snapshot depends on the default design (`_snaps`), the
snapshot is regenerated and the diff reviewed.

## 5. Documentation

- `vignettes/spiDE-model.Rmd`: a section on the nested intercept (why the
  niche slope must be within-(sample, cell type)) and one on the convergence
  stage (what `fitNB` leaves unconverged, what the polish does, the ψ
  departure).
- Argument help for `re.celltype`, `converge`, `converge.maxit`,
  `converge.tol`.
- `CLAUDE.md`: the two defaults, the ψ departure, the deferred inference cost.
- `NEWS.md` entry.

## 6. Order of work

1. Tests for the design block → `.buildRandomEffects()` / `.buildNicheDesign()`
   / `fitSpiDE()` argument → green.
2. Tests for the mixed-fit τ² names → green (should need no code).
3. Tests for `.polishFit()` (dense Newton, Schur solve, ψ, restart, dispatch)
   → `R/polish.R` → green.
4. Wire `.polishFit()` into `.fitOneBandwidth()`, slot, `show()`, validity →
   tests green.
5. Toy composition fixture → e2e tests → green.
6. Satterthwaite oracle variant.
7. Docs, NEWS, CLAUDE.md, `devtools::document()`, `devtools::test()`,
   `devtools::check()`.
8. Research: run the fixed design on the real cohort against `block`
   shuffles (`calibration-check`), and the simulation study arm.
