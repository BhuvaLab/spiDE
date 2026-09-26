# Move the polish machinery into SpaNorm, and polish SpaNorm's own model

Date: 2026-09-25. Status: proposed. Plan: `design/plans/2026-09-25-polish-to-spanorm.md`.

## 1. Why

`polishSpiDE()` converges each gene of a shared `SpaNorm::fitNB` fit to that gene's own
penalised-NB optimum. The fit is off optimum because of how `fitNB` works: it uses one
gene-averaged weight vector, an aggregate convergence criterion and a cross-gene
coefficient clamp. That is a property of SpaNorm's engine, not of spiDE's design, so
the code that corrects it belongs next to the engine. SpaNorm's own normalisation and
its SVG test run on the same engine. YTMACosMxWTAv2 (`claude/code/_polish.R`,
`31_polished_svg.R`, `claude/reports/v2_progress.md`, 2026-09-23) already runs a
temporary SpaNorm polish, built on `spiDE:::` internals.

The same pattern has been used before. The quasi-likelihood dispersion moved into
SpaNorm 1.7.10 as generic NB-GLM machinery (`qlDispersion()`), and spiDE only wires it
in.

## 2. What moves and what stays

**Moves to SpaNorm** (generic: needs only `Y, W, alpha, psi, penalty`):

- The penalised NB likelihood and mean floor: `.nbPenLoglik`, `.MU_FLOOR`.
- The per-gene Newton engine: `.polishGene` and its solver
  (`.newtonSolver`, `.absorbBlocks`, `.newtonSolverBlocked`).
- The batched engine: `R/polish-batch.R` in full, and from `R/inference-batch.R` the
  batched linear algebra that both the polish and the inference use (`.gramBatch`,
  `.absorbBatch`, `.cholBatch`, `.cholSolveBatch`, `.newtonSolverBatch`,
  `.segmentSum`) plus the backend helpers the polish uses.
- The driver without spiDE semantics: `.polishFit` minus the spiDE argument names,
  and `.reprofilePsi`.
- `.bplapplySingleBLAS` (copied, because spiDE's inference also uses it).

**Stays in spiDE** (it is spiDE's model):

- `polishSpiDE()`, `.polishSpiDEFit`, `.tau2Iterate`, `.polishPenalty`: the
  variance-component loop, the Satterthwaite and patient-bound df, and clearing the
  inference slots.
- `.absorbSpec()`: turns `re_group`/`re_sample` into an absorption spec.
- Deriving the start columns from `covtype`.
- Inference sizing (`.inferenceBlockSize`, `.covBatchSize`), `.subsetBatch`,
  `.batchQuad`, `.batchDiag`, and spiDE's own copies of the three-line backend helpers
  its inference uses (`.rowsOf`, `.colsOf`, `.asHost`).
- `longtests/testthat/test-mixed-numerics.R`, `tests/testthat/test-polish-stage.R`, and
  every test of the tau2 window. These are properties of spiDE.

## 3. New SpaNorm API

| Function | Kind | Purpose |
|---|---|---|
| `polishNB(Y, W, alpha, psi, lambda.a, offset, absorb, start.cols, psi.method, ...)` | exported | Generic per-gene polish of any `fitNB`-shaped fit |
| `nbNewtonSolver(W, pen, absorb)` | exported, `@keywords internal` | Per-gene solver: `factor()`, `solve()`, `xcov()` |
| `nbNewtonSolverBatch(W, pen, absorb)` | exported, `@keywords internal` | Batched solver on the CPU or a torch tensor |
| `nbGramBatch(W, wt_block, penalty_diag, backend, cell.tile)` | exported, `@keywords internal` | Batched penalised gram |
| `nbAbsorbGramBatch(W, pen, absorb, wt_block, cell.tile, parts)` | exported, `@keywords internal` | Batched Schur complement |
| `polishSpaNorm(spe, ...)` | exported generic, SPE + Seurat | Polishes the SpaNorm model |
| `SpaNormSVG(spe, ...)` | changed | Polishes the null whenever the full fit is polished |

`psi.method = c("profile", "fixed")` in `polishNB`. spiDE's user-facing `"moderated"`
maps to `"fixed"`: the generic function cannot assume the supplied dispersion was
moderated.

## 4. The SpaNorm model (what `polishSpaNorm` converges)

    log mu_gi = gmean_g + a1 * W_i1 + sum_{j >= 2} W_ij alpha_gj

- `W_i1` is the log library size. **`a1` is one coefficient shared by every gene**
  (`.irlsSolveAlphaFromB`: `a1 = mean_g` of the per-gene unregularised estimate). It is
  unpenalised.
- `gmean_g` is a per-gene intercept held **outside `W`** (the gmean-fold update) and is
  unpenalised.
- Columns `2..p` are penalised by `wtype`: `biology` → `lambda.a[1]`,
  `ls` → `lambda.a[2]`, `batch` → 0. The penalty is multiplied by `fit$ncells`,
  matching `fitSpaNorm` (`R/mainSpaNorm.R:278`) and, from 1.7.12, `fitSpaNormTechnical`.

The per-gene problem therefore takes `Waug = cbind(1, W[, -1])`,
`pen = c(0, lambda.vec[-1] * ncells)`, `offset = a1 * W[, 1]`, and a start of
`cbind(gmean, alpha[, -1])`. This is what the temporary `_polish.R` does. The v2
report lists four ways the first attempt got it wrong: `a1` freed per gene, no
intercept, a uniform penalty ~1,000× too weak, and an unpolished null.

**The shared coefficient: `ls = c("fixed", "joint")`.**

- `"fixed"` (the default) holds `a1` at the fit's value. This is the measured YTMA path.
- `"joint"` alternates the per-gene polish with a pooled, profiled Newton step on `a1`:
  - score `U = sum_g W1' (y_g - mu_g) / (1 + psi_g mu_g)`;
  - profiled information `I = sum_g [W1' D_g W1 - c_g' H_g^-1 c_g]`, where
    `c_g = Waug' D_g W1`;
  - line search on the total penalised log-likelihood after a warm re-polish.

The two answer different questions. SpaNorm's `a1` is an unweighted mean over genes;
the joint MLE is information-weighted, so bright genes dominate it. The joint mode
ships opt-in, and **the default is decided by a measurement** (plan, Task 9). It is
not assumed.

**Dispersion default: `"fixed"` for SpaNorm.** Profiling sent up to 46% of genes to
the `psi` floor on small cores (v2 report, 2026-09-23). The tau2 hazard that made
`"profile"` spiDE's default does not exist in SpaNorm, which has no variance
component.

**Cells: `cells = c("all", "fit")`, default `"all"`.** The penalty is `lambda * ncells`,
that is `lambda` per cell over all cells, and the normalisation and `svgTest()` both
evaluate on all cells. `"fit"` restricts to `sampling != "all"`, the cells
`fitSpaNorm` used.

## 5. SVG consistency rule

`svgTest()` compares `2(ll_full - ll_null)`. Both fits must be estimated the same way,
so polishing only the full fit inflates `F` (v2 report, 2026-09-23). The rule:

- `polishSpaNorm()` polishes `SpaNormNull` too when it is present.
- It removes stale SVG columns from `rowData`, with a warning.
- `SpaNormSVG()` fits a missing null and then polishes it with the full fit's recorded
  polish settings.
- `SpaNormSVG()` refuses a mixed pair (one fit polished, one not).

This requires the penalty fix on `fix/technical-lambda-scaling` (1.7.12). Without it
the null carries a penalty `ncells` times smaller.

## 6. Recording the polish on the fit

`SpaNormFit` gains a `polish` slot (`list`). It is empty for an unpolished fit, and
otherwise holds `settings` (psi.method, ls, cells, maxit, tol, the penalty vector and
the SpaNorm version) and `genes` (the per-gene diagnostics data.frame). Objects saved
before the slot existed are read through `methods::.hasSlot()` (spiDE's `re_sample`
pattern), and `isPolished()` is the one predicate. `polishSpaNorm()` keeps the input
fit as `SpaNormUnpolished` (and `SpaNormNullUnpolished`), the names the YTMA objects
already use.

## 7. Acceptance

1. **The move changes no number in spiDE.** Golden outputs captured before the move
   (`polishSpiDE` + `testSpiDE` on `data(toySpiDE)` and `.toyClustered()`: `alpha`,
   `psi`, `loglik`, `tau2`, `df`, `t_stat`, `se`, results) are reproduced after it to
   `tolerance = 0` on the CPU. `test-mixed-numerics.R` passes unchanged.
2. **`polishNB` is at the optimum.** The per-gene penalised score is below `1e-6`
   relative, with and without an offset, on both engines. The two engines agree to
   `1e-8`.
3. **`polishSpaNorm(ls = "fixed")` reproduces `_polish.R`** on one YTMA core to `1e-6`
   in `gmean`/`alpha`.
4. **`polishSpaNorm(ls = "joint")`** is at a stationary point of the joint objective:
   the pooled `a1` score is below `1e-6` relative, and the result matches an `optim()`
   oracle on a 5-gene problem.
5. **The SVG pair is consistent.** No mixed pair can reach `svgTest()`.
6. BiocCheck is clean on both packages. spiDE requires `SpaNorm (>= 1.7.13)`.

## 8. Out of scope

- The unpenalised `W'W` inverse in `.irlsSolveAlphaFromB`, which fails on
  rank-deficient spline designs (v2 report, df.tps ladder). This is a separate SpaNorm
  fix. `polishSpaNorm` does not invert that matrix.
- `svgTest()`'s conservative nominal-df F (92–95% of p > 0.5 in YTMA). This is a
  separate question.
- Absorbing the nested block in `.blockedInference()`. It is deferred in CLAUDE.md
  and stays in spiDE.
