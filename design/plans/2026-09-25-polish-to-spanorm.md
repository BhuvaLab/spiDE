# Polish → SpaNorm Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move the generic per-gene polish (penalised-NB Newton engine, batched linear algebra, dispersion profiling) from spiDE into SpaNorm as `polishNB()`. Add `polishSpaNorm()`, which polishes SpaNorm's own model: shared library-size coefficient, `gmean` outside `W`, `wtype` penalties. Extend the polish to the SVG null. spiDE's numbers must not change.

**Architecture:** SpaNorm owns the engine: it knows about `Y, W, alpha, psi, penalty, offset` and an absorbable indicator block. spiDE keeps the loop around it (tau2, df, `re_group` → absorb spec). `polishSpaNorm()` maps a `SpaNormFit` onto the engine: `Waug = cbind(1, W[, -1])`, `offset = a1 * W[, 1]`. It can optionally add a pooled profiled-Newton step for the shared `a1`.

**Tech Stack:** R ≥ 4.5, Bioconductor (S4, BiocParallel, SummarizedExperiment), testthat 3e, optional torch (fp64) backend, RhpcBLASctl (Suggests).

**Spec:** `design/specs/2026-09-25-polish-to-spanorm.md`. Read it first. §4 (the model) and §5 (SVG consistency) are the parts that are easy to get wrong.

## Global Constraints

- **The move is behaviour-preserving for spiDE.** The golden outputs from Task 1 must match at `tolerance = 0` on the CPU after Task 5. Do not "fix" anything in moved code during the move. File a follow-up instead.
- SpaNorm versions: 1.7.12 = the `fix/technical-lambda-scaling` merge (Task 0). 1.7.13 = the engine (Tasks 2–4). 1.7.14 = `polishSpaNorm` + SVG + joint (Tasks 6–9). spiDE `DESCRIPTION`: `SpaNorm (>= 1.7.13)`, bumped with Task 5.
- `torch` stays in `Suggests` in both packages, and nothing hard-depends on it. `RhpcBLASctl` goes to SpaNorm `Suggests`.
- Argument names are dot-separated where they mirror `fitNB` (`lambda.a`, `psi.method`, `block.size`, `batch.size`, `gpu.mem.budget`, `BPPARAM`).
- `.MU_FLOOR = exp(-30)` must remain **one** value used by both the polish and spiDE's inference. SpaNorm exports it as `nbMuFloor()` (Task 3), and spiDE calls that. spiDE never uses `SpaNorm:::` (BiocCheck flags it).
- **Never install a development SpaNorm into `~/R/x86_64-conda-linux-gnu-library/4.5`.** Running HPC jobs lazy-load from it, and overwriting a package under a running R process corrupts its lazy-load DB. Every worktree installs into its own library: `R_LIBS=$WT/.Rlib Rscript -e 'devtools::install("../SpaNorm-<wt>", upgrade = "never")'`.
- Only the **integrator** (the main session) runs `devtools::document()` on a shared branch, edits `DESCRIPTION`/`NAMESPACE`/`NEWS.md`, bumps versions and merges. Task agents run `document()` in their own worktree to test, but commit only their `R/`, `man/` and `tests/` files. The integrator regenerates `NAMESPACE` at merge.
- Every numerical task ends with a `numerical-robustness-reviewer` pass. Every spiDE `R/` change also gets `design-invariant-reviewer`.

## Review Focus

1. **A fit whose `alpha[, 1]` is not shared.** This happens with a `fitNB` result or a hand-built fit passed to `polishSpaNorm`. Expected: a clear error naming the column. Never silently polish `a1` per gene; that was the v2 mis-specification. Test in Task 7.
2. **Polishing the full fit while the saved `SpaNormNull` is unpolished, then calling `SpaNormSVG()`.** Expected: the null is polished with the same settings or the call refuses. `svgTest()` must never see a mixed pair. Test in Task 8.
3. **A non-integer counts assay** (e.g. `2^logcounts - 1`). Expected: the integer guard stops with the existing message and never returns `psi` pinned at the upper bound. Test in Task 4 (engine) and Task 7 (entry point).
4. **Genes that are zero in every cell** (49 of 305 YTMA nuclear cores had up to 12). Expected: a finite result for the other genes. The zero gene is flagged `polished = FALSE`, keeps its input coefficients, and does not abort the core. Test in Task 7.
5. **A `SpaNormFit` saved before the `polish` slot existed**, read from RDS. Expected: `isPolished()` is `FALSE` and `polishSpaNorm()` works on it. Test in Task 6.

---

## Multi-agent execution: the assessment

**The dependency graph.** "Polish" below is `polishSpaNorm`.

```
T0 (gate: land slope branch; merge SpaNorm 1.7.12) ─┐
T1 spiDE golden fixtures ───────────────────────────┤
T2 SpaNorm linear algebra ── T3 engine + offset ── T4 polishNB ──┬── T5 spiDE rewire
T6 SpaNormFit polish slot + model mapping ───────────────────────┴── T7 polishSpaNorm (fixed) ──┬── T8 SVG consistency
                                                                                                 └── T9 joint a1
                                                                        T8 + T9 ── T10 measurement (HPC) ── T11 docs + release
```

**Waves and agent count.**

| Wave | Parallel tasks | Agents | Repos / files touched | Collision risk |
|---|---|---|---|---|
| 1 | T1, T2, T6 | 3 | spiDE `tests/testthat/_golden/`; SpaNorm `R/nbSolver*.R`; SpaNorm `R/AllClasses.R`, `R/polishSpaNormModel.R` | none: disjoint files. `NAMESPACE` is resolved by the integrator |
| 2 | T3 → T4 | 1 | SpaNorm `R/polishEngine*.R`, `R/polishNB.R` | sequential by necessity: T4 wraps the T3 engine |
| 3 | T5, T7 | 2 | spiDE `R/polish*.R`, `R/inference*.R`; SpaNorm `R/polishSpaNorm.R` | none: different repos. Both install SpaNorm 1.7.13 into their own lib |
| 4 | T8, T9 | 2 | SpaNorm `R/mainSpaNormSVG.R`; SpaNorm `R/polishSpaNormJoint.R` | low. T7 defines the `ls` switch and the `.polishSharedLS()` hook, so T9 fills a stub and does not edit T7's file |
| 5 | T10 | 1 + you | HPC | a human gate: it decides a default |
| 6 | T11 | 1 + `evidence-auditor` | docs in both repos | none |

**Is it worth it?** Only moderately. The critical path T2→T3→T4→T7→T8/T9→T10 is serial and holds about 70% of the work. Parallelism saves roughly a third of the wall clock (waves 1, 3 and 4). It does not remove the real risk, which is in the **seams**: a renamed argument, a lost `.MU_FLOOR`, the per-block moments table, a dropped singular guard. That argues for:

- **At most three concurrent implementers**, each in its own worktree (`superpowers:using-git-worktrees`): `spiDE-wt-<task>` and `SpaNorm-wt-<task>`.
- **A fresh reviewer per task** (`superpowers:subagent-driven-development`), plus the project's `numerical-robustness-reviewer`, and `design-invariant-reviewer` on T5.
- **One integrator** (this session). It merges worktrees in wave order, owns `NAMESPACE`/`DESCRIPTION`/versions, installs each merged SpaNorm into a shared *dev* library (`$SCRATCH/Rlib_polishmove`) that the next wave's worktrees copy, and runs the golden test after every spiDE merge.
- **Not** a many-agent Workflow fan-out. The tasks are small in number and interface-coupled, and a Workflow needs your explicit opt-in anyway.

**What the agents need from each other.** Each task's **Interfaces** block below is the whole contract. An implementer sees only its own task, so any name used across tasks is defined once there and copied verbatim.

---

### Task 0: Gates (integrator + you)

**Files:** none new.

- [ ] **Step 1: Land `feature/slope-schur-absorption` in spiDE.** The absorption code is what moves, and its working tree has uncommitted edits to `R/polish.R`, `R/inference.R`, `R/design.R`, `R/AllClasses.R` and `R/fitSpiDE.R`. Finish, review and merge that branch to `main` first. Nothing below starts from a tree with those edits pending.
- [ ] **Step 2: Merge SpaNorm `fix/technical-lambda-scaling` (f07f723, Version 1.7.12) to `master` and push.**

```bash
cd /scratch/project_mnt/S0249/R_projects/SpaNorm
git checkout master && git merge --no-ff fix/technical-lambda-scaling
Rscript -e 'devtools::test()'
git push origin master
```
Expected: all tests pass, including the regression test that fails on 1.7.11.

- [ ] **Step 3: Create the shared dev library.**

```bash
mkdir -p $SCRATCH/Rlib_polishmove
R_LIBS=$SCRATCH/Rlib_polishmove Rscript -e 'devtools::install("/scratch/project_mnt/S0249/R_projects/SpaNorm", upgrade = "never")'
```

---

### Task 1: spiDE golden fixtures (the safety net for the move)

**Files:**
- Create: `tests/testthat/helper-golden.R` (`golden_polish_fits()`)
- Create: `tests/testthat/_golden/make_golden_polish.R`
- Create: `tests/testthat/_golden/golden_polish.rds` (generated, committed)
- Create: `tests/testthat/test-polish-golden.R`

**Interfaces:**
- Produces: `golden_polish.rds`, a named list with `toy`, `clustered` and `nichemode`. Each holds `alpha`, `psi`, `loglik`, `tau2`, `penalty`, `df`, `t_stat`, `se`, `polish` (per bandwidth) and `results` (data.frame).

- [ ] **Step 1: Write the generator.**

```r
# tests/testthat/_golden/make_golden_polish.R
# Golden outputs of polish + inference, captured BEFORE the polish machinery moved
# to SpaNorm (design/plans/2026-09-25-polish-to-spanorm.md). Regenerate ONLY on a
# deliberate numerical change, never to make the move pass.
devtools::load_all(quiet = TRUE)
snap <- function(res) {
  fits <- res@fits
  list(
    fits = lapply(fits, function(f) list(
      alpha = f@alpha, psi = f@psi, loglik = f@loglik, tau2 = f@tau2,
      penalty = f@penalty, df = f@df, t_stat = f@t_stat, se = f@se,
      polish = f@polish)),
    results = results(res))
}
set.seed(20260925)
data(toySpiDE)
toy <- spiDE(toySpiDE, condition = "condition", sigma = c(30, 50),
             BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
set.seed(20260925)
cl <- .toyClustered(sd_nested = 0.3)
clustered <- spiDE(cl, condition = "condition", sigma = 30,
                   BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
set.seed(20260925)
niche <- spiDE(toySpiDE, condition = NULL, sigma = 30,
               BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
saveRDS(list(toy = snap(toy), clustered = snap(clustered), nichemode = snap(niche)),
        "tests/testthat/_golden/golden_polish.rds")
```

Before running it, check the `spiDE()` argument names (`condition`, `sigma`, the `.toyClustered()` arguments) against `R/spiDE.R` and `R/toydata.R` on the merged `main`, and fix any mismatch in the generator. The accessor is `@fits` if `R/AllClasses.R` names it so; use the actual slot name.

- [ ] **Step 2: Write the test.**

Split the generator so the test can call it: move everything above `saveRDS(...)`
into `tests/testthat/helper-golden.R` as `golden_polish_fits()` (it returns the list).
`make_golden_polish.R` then reduces to:

```r
devtools::load_all(quiet = TRUE)
saveRDS(golden_polish_fits(), "tests/testthat/_golden/golden_polish.rds")
```

(The `devtools::load_all()` line leaves the helper; testthat loads helpers itself.)

```r
# tests/testthat/test-polish-golden.R
test_that("polish + inference reproduce the pre-move golden outputs exactly", {
  skip_on_cran()
  g <- readRDS(test_path("_golden", "golden_polish.rds"))
  now <- golden_polish_fits()
  for (nm in names(g)) {
    expect_identical(names(now[[nm]]$fits), names(g[[nm]]$fits), info = nm)
    for (bw in names(g[[nm]]$fits)) {
      for (s in names(g[[nm]]$fits[[bw]])) {
        expect_equal(now[[nm]]$fits[[bw]][[s]], g[[nm]]$fits[[bw]][[s]],
                     tolerance = 0, info = paste(nm, bw, s))
      }
    }
    expect_equal(now[[nm]]$results, g[[nm]]$results, tolerance = 0, info = nm)
  }
})
```

- [ ] **Step 3: Generate and run.**

Run: `Rscript tests/testthat/_golden/make_golden_polish.R && Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-polish-golden.R")'`
Expected: PASS. Record the wall time; the test must stay under ~2 min or it is moved to `longtests/`.

- [ ] **Step 4: Check determinism.** Run the test twice more. Any failure means an unseeded path (CLAUDE.md: `fitNB` subsamples for dispersion above a size threshold). Fix by seeding in the generator, not by loosening the tolerance.

- [ ] **Step 5: Commit.**

```bash
git add tests/testthat/_golden tests/testthat/test-polish-golden.R
git commit -m "tests: golden polish+inference outputs, the safety net for moving the polish to SpaNorm"
```

---

### Task 2: SpaNorm, the batched and per-gene linear algebra

**Files:**
- Create: `SpaNorm/R/nbSolver.R` ← spiDE `R/polish.R`: `.absorbBlocks`, `.newtonSolverBlocked`, `.newtonSolver`, and their roxygen/comment blocks verbatim.
- Create: `SpaNorm/R/nbSolverBatch.R` ← spiDE `R/inference-batch.R`: `.gramBatch`, `.segmentSum`, `.absorbBatch`, `.cholBatch`, `.cholSolveBatch`, `.newtonSolverBatch`, and the backend helpers `.rowsOf`, `.setRows`, `.mulRows`, `.scaleCols`, `.asHost`, `.matmulB`, `.rowsFinite`, `.rowMeansB`, `.colsOf`, `.maskedRowMin`, `.rowMaxB`, `.asHostMat`, `.asLike`, `.subsetState`.
- Create: `SpaNorm/tests/testthat/test-nbSolver.R` ← spiDE `test-solver-batch.R`, `test-solver-slope-absorb.R`, `test-absorb-batch.R`, keeping only the tests that build `W` directly. Tests that go through `fitSpiDE` stay in spiDE.
- Modify: `SpaNorm/R/gpuFunctions.R` only if a moved helper duplicates an existing SpaNorm helper. In that case keep SpaNorm's version and adapt the call.

**Interfaces:**
- Consumes: SpaNorm's own `is_torch_tensor`, `toRMatrix`, `toGPUMatrix`, `invert_mat_batched`.
- Produces (exported, `@keywords internal`, one man page `nbSolver.Rd`):
  - `nbNewtonSolver(W, pen, absorb = NULL)`: the old `.newtonSolver(W, pen, nested)`. It returns `list(factor = function(w), solve = function(state_or_w, s), xcov = function(state_or_w))`. `solve`/`factor` return `NULL` on a singular system.
  - `nbNewtonSolverBatch(W, pen, absorb = NULL)`: the old `.newtonSolverBatch`.
  - `nbGramBatch(W, wt_block, penalty_diag = NULL, backend = "cpu", cell.tile = NULL)`: the old `.gramBatch`.
  - `nbAbsorbGramBatch(W, pen, absorb, wt_block, cell.tile = NULL, parts = FALSE)`: the old `.absorbBatch`.
  - The internal names stay as they are (`.newtonSolver` etc.). The exported names are thin aliases (`nbNewtonSolver <- function(W, pen, absorb = NULL) .newtonSolver(W, pen, absorb)`), so the moved bodies change only `SpaNorm::` → local calls.

- [ ] **Step 1: Copy the functions verbatim** into the two new files. Then, in those files only, replace every `SpaNorm::` prefix with nothing. Rename nothing else.

- [ ] **Step 2: Add the exported aliases and roxygen.**

```r
#' Penalised NB Newton solvers and batched grams (for package developers)
#'
#' Low-level building blocks of \code{\link{polishNB}}, exported so that a
#' downstream package (e.g. spiDE) can form the penalised covariance of a
#' polished fit from the SAME factorisation the polish used. The information is
#' \code{W' diag(w) W + diag(pen)}. \code{absorb} marks an indicator block whose
#' part of the information is block-diagonal and is eliminated by a Schur
#' complement: \code{NULL} (dense), a logical over \code{W}'s columns (each
#' marked column its own 1x1 block), or an integer block id per column
#' (\code{NA} = dense).
#'
#' @param W a cells x p design (base matrix, or torch tensor for the batch forms).
#' @param pen a length-p ridge penalty.
#' @param absorb see Description.
#' @param wt_block a genes x cells weight matrix (one row per gene).
#' @param penalty_diag,backend,cell.tile,parts see \code{polishNB()}.
#' @return \code{nbNewtonSolver()}: a list of closures \code{factor(w)},
#'   \code{solve(state, s)} (the Newton step, \code{NULL} if singular) and
#'   \code{xcov(state)} (the dense-block covariance). The batch forms return
#'   the batched grams or Schur complements.
#' @examples
#' W <- cbind(1, seq(-1, 1, length.out = 20))
#' s <- nbNewtonSolver(W, pen = c(0, 1))
#' st <- s$factor(rep(1, 20))
#' s$solve(st, c(1, 0))
#' @keywords internal
#' @export
nbNewtonSolver <- function(W, pen, absorb = NULL) .newtonSolver(W, pen, absorb)
```

Add the same roxygen `@rdname nbNewtonSolver` pattern for the other three aliases.

- [ ] **Step 3: Move the tests.** Rewrite `spiDE:::` → the SpaNorm internal names. Delete any test that calls `fitSpiDE`/`toySpiDE`. List those in the commit message; they stay in spiDE and are rewired in Task 5.

- [ ] **Step 4: Run.**

Run: `cd SpaNorm && Rscript -e 'devtools::document(); devtools::test(filter = "nbSolver")'`
Expected: PASS, including the torch-backed tests, which skip without libtorch (the existing SpaNorm helper pattern `helper-gpu.R`).

- [ ] **Step 5: Review, then commit on the worktree branch.**

Dispatch `numerical-robustness-reviewer` on `R/nbSolver*.R`. The singular-guard sites must all have survived. `invert_mat_batched` is the one a grep for `invert_mat(` misses.

```bash
git add R/nbSolver.R R/nbSolverBatch.R tests/testthat/test-nbSolver.R man/nbSolver.Rd
git commit -m "Penalised NB Newton solvers and batched grams, moved from spiDE"
```

---

### Task 3: SpaNorm, the polish engine with an offset

**Files:**
- Create: `SpaNorm/R/polishEngine.R` ← spiDE `R/polish.R`: `.MU_FLOOR`, `.nbPenLoglik`, `.polishGene`.
- Create: `SpaNorm/R/polishEngineBatch.R` ← spiDE `R/polish-batch.R`: all of it, with `SPIDE_POLISH_GENE_CELL_MATS` renamed `POLISH_GENE_CELL_MATS`; plus spiDE's `.covBatchSize` and the constants it reads, copied as `.polishCovBatchSize` (spiDE keeps its own for inference).
- Create: `SpaNorm/tests/testthat/test-polishEngine.R` ← spiDE `test-polish-batch.R`, `test-psi-batch.R`, `test-nb-batch.R` (the direct-`W` tests), plus the new offset tests below.

**Interfaces:**
- Consumes: Task 2's internal `.newtonSolver`, `.newtonSolverBatch`, and the backend helpers.
- Produces:
  - `.polishGene(y, W, a0, psi0, pen, solver, maxit, tol, start.cols = NULL, psi.range, psi.method = c("profile", "fixed"), warm = FALSE, offset = NULL)`. `ct_cols` is renamed `start.cols`, and `"moderated"` is renamed `"fixed"`. `offset` is `NULL` or a length-`length(y)` numeric.
  - `.polishBatch(Yb, W, A0, psi0, pen, solver, maxit, tol, start.cols = NULL, psi.range, psi.method = c("profile", "fixed"), warm = FALSE, shared.factor = FALSE, nested = NULL, offset = NULL)`. `offset` is `NULL`, a length-`ncol(Yb)` vector (the same for every gene), or a `nrow(Yb) x ncol(Yb)` matrix.
  - `nbMuFloor()` (exported, returns `exp(-30)`) with the rationale comment moved from `R/polish.R:30-38`.

- [ ] **Step 1: Write the failing offset tests first.**

```r
# tests/testthat/test-polishEngine.R (new block)
.score <- function(y, W, a, psi, pen, off) {
  mu <- pmax(as.numeric(exp(W %*% a + off)), nbMuFloor())
  as.numeric(crossprod(W, (y - mu) / (1 + psi * mu))) - pen * a
}

test_that("polishGene reaches a zero penalised score with an offset", {
  set.seed(1)
  n <- 300
  W <- cbind(1, rnorm(n), rnorm(n))
  off <- log(runif(n, 0.5, 2))               # a known per-cell log-scale factor
  y <- rnbinom(n, mu = exp(1 + 0.4 * W[, 2] + off), size = 5)
  pen <- c(0, 1, 1)
  s <- .newtonSolver(W, pen)
  r <- .polishGene(y, W, a0 = c(0, 0, 0), psi0 = 0.2, pen = pen, solver = s,
                   psi.method = "fixed", offset = off)
  sc <- .score(y, W, r$alpha, r$psi, pen, off)
  expect_lt(max(abs(sc)) / sum(y), 1e-6)
  expect_true(r$polished)
})

test_that("the sane start subtracts the mean offset", {
  set.seed(2)
  n <- 200
  W <- cbind(1, rnorm(n))
  off <- rep(3, n)                            # a large constant offset
  y <- rnbinom(n, mu = exp(0.5 + off), size = 10)
  s <- .newtonSolver(W, c(0, 0))
  # a hopeless start forces the restart to the sane start
  r <- .polishGene(y, W, a0 = c(40, 0), psi0 = 0.1, pen = c(0, 0), solver = s,
                   psi.method = "fixed", offset = off)
  expect_equal(r$alpha[1], 0.5, tolerance = 0.1)
})

test_that("batch and per-gene engines agree with a vector and a matrix offset", {
  set.seed(3)
  G <- 6; n <- 150
  W <- cbind(1, rnorm(n))
  off <- log(runif(n, 0.5, 2))
  Y <- t(vapply(seq_len(G), function(g)
    rnbinom(n, mu = exp(0.5 + 0.2 * g * W[, 2] + off), size = 4), numeric(n)))
  pen <- c(0, 0.5)
  s <- .newtonSolver(W, pen)
  per <- t(vapply(seq_len(G), function(g)
    .polishGene(Y[g, ], W, c(0, 0), 0.25, pen, s, psi.method = "fixed",
                offset = off)$alpha, numeric(2)))
  bv <- .polishBatch(Y, W, matrix(0, G, 2), rep(0.25, G), pen, s,
                     psi.method = "fixed", offset = off)
  bm <- .polishBatch(Y, W, matrix(0, G, 2), rep(0.25, G), pen, s,
                     psi.method = "fixed",
                     offset = matrix(off, G, n, byrow = TRUE))
  expect_equal(bv$alpha, per, tolerance = 1e-8)
  expect_equal(bm$alpha, bv$alpha, tolerance = 0)
})

test_that("offset = NULL is bit-identical to the pre-offset engine", {
  set.seed(4)
  n <- 120
  W <- cbind(1, rnorm(n))
  y <- rnbinom(n, mu = exp(1 + 0.3 * W[, 2]), size = 3)
  s <- .newtonSolver(W, c(0, 0.1))
  a <- .polishGene(y, W, c(0, 0), 0.3, c(0, 0.1), s, psi.method = "profile")
  b <- .polishGene(y, W, c(0, 0), 0.3, c(0, 0.1), s, psi.method = "profile",
                   offset = NULL)
  expect_identical(a, b)
})
```

Check the return field names (`alpha`, `psi`, `polished`) of `.polishGene`/`.polishBatch` against the moved code and adjust the tests to them, not the reverse.

- [ ] **Step 2: Run them to see them fail.**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-polishEngine.R")'`
Expected: FAIL with `unused argument (offset = off)`.

- [ ] **Step 3: Move the code verbatim, then add the offset.** In `.polishGene`:
  - Normalise once at the top: `off <- if (is.null(offset)) 0 else as.numeric(offset)`, with a length check that stops with `"'offset' must have one value per cell"`.
  - Every linear-predictor site `W %*% a` becomes `W %*% a + off`. Find them with `grep -n '%\*%' R/polishEngine.R`. The `newton()` closure's initial `mu`, the line-search candidate `mu`, and the restart check are all sites.
  - The sane start becomes `a[j] <- if (any(cells)) log(mean(y[cells]) + 1e-3) - mean(off0[cells]) else 0` and `a[1] <- log(mean(y) + 1e-3) - mean(off0)`, where `off0 <- rep_len(off, length(y))`.
  - `psi.method` values: `"moderated"` → `"fixed"` everywhere in the moved files.

  In `.polishBatch`, add a helper and use it at every `Eta <- .matmulB(...)` / `.muBatch(...)` site:

```r
# The offset for the genes in `kk` (row indices into the batch), shaped like
# Eta: NULL, a per-cell vector (identical for every gene) or a genes x cells
# matrix. Every linear predictor in the batch engine goes through this one
# function so a mis-sliced offset cannot hide at one site.
.offsetRows <- function(offset, kk, like) {
  if (is.null(offset)) return(NULL)
  if (is.null(dim(offset))) {
    return(.asLike(matrix(offset, length(kk), length(offset), byrow = TRUE), like))
  }
  .asLike(.rowsOf(offset, kk), like)
}
```
  Then write `Eta <- .matmulB(.rowsOf(A, kk), tW); o <- .offsetRows(offset, rows[kk], Eta); if (!is.null(o)) Eta <- Eta + o`. `.muBatch(A, W)` gains `offset = NULL` with the same addition. The sane start in the batch engine subtracts the row means of the offset, as in the per-gene engine.

- [ ] **Step 4: Run all engine tests.**

Run: `Rscript -e 'devtools::document(); devtools::test(filter = "polishEngine|nbSolver")'`
Expected: PASS.

- [ ] **Step 5: Review, then commit.** Run `numerical-robustness-reviewer` on `R/polishEngine*.R`. Focus on the offset at every mu site, the NaN guards and the singular fallbacks.

```bash
git commit -am "Per-gene and batched polish engine with an optional offset, moved from spiDE"
```

---

### Task 4: SpaNorm, the exported `polishNB()`

**Files:**
- Create: `SpaNorm/R/polishNB.R` ← spiDE `R/polish.R`: `.polishFit` (the driver body), `.reprofilePsi`, and `.bplapplySingleBLAS`/`.singleBLAS`/`.workerBLAS` from spiDE `R/inference.R`.
- Modify: `SpaNorm/DESCRIPTION`: `RhpcBLASctl` → `Suggests`; `Version: 1.7.13`. **Integrator only.**
- Create: `SpaNorm/tests/testthat/test-polishNB.R`

**Interfaces:**
- Consumes: Task 3's `.polishGene`, `.polishBatch`, `.polishBatchSize`; Task 2's solvers.
- Produces (exported):

```r
polishNB(Y, W, alpha, psi, lambda.a = 0, offset = NULL, absorb = NULL,
         start.cols = NULL, psi.method = c("profile", "fixed"),
         psi.range = c(1e-3, 1e3), warm = FALSE, maxit = 50L, tol = 1e-8,
         engine = c("batch", "gene"), batch.size = NULL, block.size = NULL,
         backend = c("cpu", "auto", "gpu"), gpu.mem.budget = NULL,
         BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
# returns list(alpha = <genes x p>, psi = <genes>, loglik = <genes>,
#              polish = data.frame(iterations, restarted, capped, singular,
#                                  psi_bound, polished, row.names = rownames(Y)))
```
  `lambda.a` is a scalar or a length-`ncol(W)` vector, **applied as given with no scaling**. The caller owns the scaling; say so in the docs. `start.cols` is a logical over the columns of `W`, or `NULL`.

- [ ] **Step 1: Write the failing tests.**

```r
# tests/testthat/test-polishNB.R
test_that("polishNB converges a fitNB fit to each gene's own optimum", {
  set.seed(10)
  G <- 20; n <- 400
  W <- cbind(1, rnorm(n), rnorm(n))
  Y <- t(vapply(seq_len(G), function(g)
    rnbinom(n, mu = exp(0.5 + g / 10 + 0.3 * W[, 2]), size = 3), numeric(n)))
  fit <- fitNB(Y, W, lambda.a = 0.1, verbose = FALSE, backend = "cpu")
  pol <- polishNB(Y, W, fit$alpha, fit$psi, lambda.a = 0.1, psi.method = "fixed")
  for (g in seq_len(G)) {
    mu <- exp(as.numeric(W %*% pol$alpha[g, ]))
    sc <- crossprod(W, (Y[g, ] - mu) / (1 + pol$psi[g] * mu)) - 0.1 * pol$alpha[g, ]
    expect_lt(max(abs(sc)) / sum(Y[g, ]), 1e-6)
  }
  expect_true(all(pol$polish$polished))
  expect_equal(pol$psi, fit$psi)                  # "fixed" keeps the dispersion
})

test_that("gene and batch engines agree", {
  set.seed(11)
  G <- 8; n <- 250
  W <- cbind(1, rnorm(n))
  Y <- t(vapply(seq_len(G), function(g) rnbinom(n, mu = exp(1 + 0.2 * W[, 2]), size = 2), numeric(n)))
  a0 <- matrix(0, G, 2); p0 <- rep(0.5, G)
  a <- polishNB(Y, W, a0, p0, lambda.a = c(0, 1), engine = "gene")
  b <- polishNB(Y, W, a0, p0, lambda.a = c(0, 1), engine = "batch")
  expect_equal(a$alpha, b$alpha, tolerance = 1e-8)
  expect_equal(a$psi, b$psi, tolerance = 1e-6)
})

test_that("polishNB refuses non-integer counts", {
  Y <- matrix(c(1.5, 2, 3, 4), 1)
  expect_error(polishNB(Y, cbind(rep(1, 4)), matrix(0, 1, 1), 0.1),
               "needs integer counts")
})

test_that("a gene with no counts is flagged, not fatal", {
  set.seed(12)
  n <- 100
  W <- cbind(1, rnorm(n))
  Y <- rbind(rnbinom(n, mu = 5, size = 3), rep(0, n))
  r <- polishNB(Y, W, matrix(0, 2, 2), c(0.3, 0.3), lambda.a = c(0, 1))
  expect_true(r$polish$polished[1])
  expect_true(all(is.finite(r$alpha[1, ])))
  expect_false(isTRUE(r$polish$polished[2]) && any(!is.finite(r$alpha[2, ])))
})
```

- [ ] **Step 2: Run them to see them fail.** Expected: `could not find function "polishNB"`.

- [ ] **Step 3: Implement.** Move `.polishFit` into `polishNB()`:
  - Rename `pen` → `lambda.a` in the signature and keep `pen` internally.
  - `absorb` replaces `re_group`/`nested`. The body uses only the absorb spec.
  - `start.cols` replaces `covtype`. The spiDE driver derived `ct_cols` from `covtype == "CellType"`; that derivation moves to spiDE in Task 5.
  - `offset` is threaded to both engines. The gene-block slicing uses the same `.rowsOf`, so a matrix offset is sliced with `Y`.
  - The integer guard and the GPU/`engine` checks move verbatim, message text included.
  - Keep `.reprofilePsi` internal.

  If the zero-gene test exposes a crash in the moved code, **do not fix it in this task**. Record it and fix it in a separate commit after the golden test has passed in Task 5. The move must stay pure.

- [ ] **Step 4: Run.**

Run: `Rscript -e 'devtools::document(); devtools::test()'`
Expected: the whole SpaNorm suite PASSES.

- [ ] **Step 5: Review, then commit.** `numerical-robustness-reviewer` on `R/polishNB.R`. The integrator bumps the version to 1.7.13, merges, pushes and installs into `$SCRATCH/Rlib_polishmove`.

---

### Task 5: spiDE, rewire onto SpaNorm and delete the moved code

**Files:**
- Modify: `R/polish.R`. Delete the moved functions. `.polishSpiDEFit`'s `run_polish()` calls `SpaNorm::polishNB(...)`.
- Modify: `R/polish-batch.R`. Delete it; everything in it moved.
- Modify: `R/inference-batch.R`. Delete the moved functions, keeping `.inferenceBlockSize`, `.covBatchSize`, `.subsetBatch`, `.batchQuad`, `.batchDiag`, the `SPIDE_*` constants, and local copies of `.rowsOf`, `.colsOf`, `.asHost`.
- Modify: `R/inference.R:410-420, 640-680`: `.absorbBatch` → `SpaNorm::nbAbsorbGramBatch`, `.gramBatch` → `SpaNorm::nbGramBatch`, `.newtonSolver` → `SpaNorm::nbNewtonSolver`, `.MU_FLOOR` → `SpaNorm::nbMuFloor()`.
- Modify: `R/design.R:88` (a comment only) and any other `.newtonSolver` reference.
- Modify: `DESCRIPTION`: `SpaNorm (>= 1.7.13)`; version bump. **Integrator only.**
- Modify/Delete tests: the moved tests are deleted. The spiDE-coupled tests in `test-polish.R`, `test-polish-batch.R`, `test-inference.R`, `test-gpuPolish.R`, `test-slope-absorb-wiring.R` and `longtests/testthat/test-nested-intercept.R` are rewired to call `SpaNorm::` names.

**Interfaces:**
- Consumes: `SpaNorm::polishNB`, `SpaNorm::nbNewtonSolver`, `SpaNorm::nbNewtonSolverBatch`, `SpaNorm::nbGramBatch`, `SpaNorm::nbAbsorbGramBatch`, `SpaNorm::nbMuFloor` (Tasks 2–4).
- Produces: no new API. `polishSpiDE()` keeps its signature, and `psi = "moderated"` maps to `psi.method = "fixed"`.

- [ ] **Step 1: Confirm the golden test passes before editing** (installed SpaNorm 1.7.13 in the task lib):
Run: `R_LIBS=$WT/.Rlib Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-polish-golden.R")'`
Expected: PASS.

- [ ] **Step 2: Rewire `run_polish()`.**

```r
  run_polish <- function(alpha0, psi0, pen_now, warm = FALSE) {
    SpaNorm::polishNB(
      Yf, f@W, alpha0, psi0, lambda.a = pen_now,
      absorb = .absorbSpec(f),
      start.cols = .testedStartCols(f),
      psi.method = if (psi == "moderated") "fixed" else "profile",
      warm = warm, maxit = maxit, tol = tol, engine = engine,
      batch.size = batch.size, block.size = block.size, backend = backend,
      gpu.mem.budget = gpu.mem.budget, BPPARAM = BPPARAM, verbose = verbose)
  }
```

Add `.testedStartCols()` next to `.absorbSpec()` in `R/polish.R`. It is the single place that derives the start columns from `covtype`, replacing the `covtype = as.character(f@covtype)` → `ct_cols` derivation that lived inside `.polishFit`. Copy that derivation's exact expression from the pre-move `.polishFit`. Any difference breaks the golden test. The `psi` argument name above is `polishSpiDE`'s own, forwarded through `.polishSpiDEFit`; follow the actual variable name there.

- [ ] **Step 3: Rewire inference and delete the moved code.** Delete in one commit and rewire in the next, so a reviewer can read the rewire diff on its own.

- [ ] **Step 4: Run the golden test, the full suite and the numerics longtest.**

Run:
```bash
R_LIBS=$WT/.Rlib Rscript -e 'devtools::document(); devtools::test()'
R_LIBS=$WT/.Rlib Rscript -e 'devtools::load_all(); testthat::test_file("longtests/testthat/test-mixed-numerics.R"); testthat::test_file("longtests/testthat/test-nested-intercept.R")'
```
Expected: all PASS, and the golden test at `tolerance = 0`. A failure here means a seam; bisect by reverting the rewire commit.

- [ ] **Step 5: Review, then commit.** Dispatch `design-invariant-reviewer` (fit whole / infer blocked, the batched inversion site, the polish clearing inference) and `numerical-robustness-reviewer`. Then run `R CMD check` + `BiocCheck::BiocCheck()`. `:::` into SpaNorm must not appear anywhere: `grep -rn "SpaNorm:::" R/` returns nothing.

---

### Task 6: SpaNorm, the `polish` slot and the SpaNorm-model mapping

Runs in wave 1, parallel to T2.

**Files:**
- Modify: `SpaNorm/R/AllClasses.R`: add the slot `polish = "list"` (prototype `list()`), a `show` line, and `isPolished()`.
- Create: `SpaNorm/R/polishSpaNormModel.R`: `.spaNormPolishProblem()`.
- Create: `SpaNorm/tests/testthat/test-polishSpaNormModel.R`

**Interfaces:**
- Produces:
  - `isPolished(fit)` (exported): `TRUE` iff `methods::.hasSlot(fit, "polish") && length(fit@polish) > 0`.
  - `.polishSlot(fit)` (internal): `fit@polish` if the slot exists, else `list()`.
  - `.spaNormPolishProblem(fit, cells = c("all", "fit"))` (internal) returns
    `list(X, pen, a1, offset, A0, psi, cells_idx, w1)`:
    - `X = cbind("(gmean)" = 1, W[cells_idx, -1])`
    - `pen = c(0, lambda.vec[-1] * fit$ncells)` (see below)
    - `a1` = the shared LS coefficient (scalar)
    - `offset = a1 * W[cells_idx, 1]`
    - `A0 = cbind(fit$gmean, fit$alpha[, -1])`
    - `w1 = W[cells_idx, 1]`
    - `cells_idx` is logical over all cells

```r
# The penalty SpaNorm's fitters use (fitSpaNorm, R/mainSpaNorm.R:258-278, and
# from 1.7.12 fitSpaNormTechnical): biology -> lambda.a[1], ls -> lambda.a[2],
# batch -> 0, times ncells; column 1 (the shared logLS) is not penalised and is
# not in the per-gene problem at all.
.spaNormPenalty <- function(fit) {
  lam <- rep_len(fit$lambda.a, 2L)
  v <- numeric(length(fit$wtype))
  v[fit$wtype == "biology"] <- lam[1]
  v[fit$wtype == "ls"] <- lam[2]
  c(0, v[-1] * fit$ncells)          # leading 0: the unpenalised gmean column
}

.spaNormPolishProblem <- function(fit, cells = c("all", "fit")) {
  cells <- match.arg(cells)
  a1s <- fit$alpha[, 1]
  if (diff(range(a1s)) > 1e-8 * max(1, abs(a1s[1]))) {
    stop("column 1 of the fit's alpha (the library-size coefficient) is not ",
         "shared across genes, so this is not a SpaNorm-model fit; use ",
         "polishNB() for a generic fitNB() fit", call. = FALSE)
  }
  idx <- if (cells == "all") rep(TRUE, nrow(fit$W)) else fit$sampling != "all"
  W <- fit$W[idx, , drop = FALSE]
  a1 <- a1s[1]
  list(X = cbind(`(gmean)` = 1, W[, -1, drop = FALSE]),
       pen = .spaNormPenalty(fit), a1 = a1, offset = a1 * W[, 1], w1 = W[, 1],
       A0 = cbind(fit$gmean, fit$alpha[, -1, drop = FALSE]),
       psi = fit$psi, cells_idx = idx)
}
```

- [ ] **Step 1: Write the failing tests.**

```r
# tests/testthat/test-polishSpaNormModel.R
mk_fit <- function(shared = TRUE) {
  set.seed(5)
  n <- 50; G <- 4
  W <- cbind(logLS = rnorm(n), b1 = rnorm(n), b2 = rnorm(n), l2 = rnorm(n))
  al <- matrix(rnorm(G * 4), G, 4)
  al[, 1] <- if (shared) 1.02 else seq(0.9, 1.1, length.out = G)
  SpaNormFit(ngenes = G, ncells = n, gene.model = "nb", df.tps = 2L,
             sample.p = 0.5, lambda.a = c(1e-4, 2e-4), batch = NULL, W = W,
             alpha = al, gmean = rnorm(G), psi = rep(0.2, G),
             wtype = factor(c("ls", "biology", "biology", "ls")),
             loglik = 0, sampling = factor(rep(c("glm", "all"), each = n / 2)))
}

test_that("the problem maps SpaNorm's model: intercept, offset, wtype penalty", {
  f <- mk_fit()
  p <- .spaNormPolishProblem(f)
  expect_equal(unname(p$X[, 1]), rep(1, 50))
  expect_equal(p$pen, c(0, 1e-4 * 50, 1e-4 * 50, 2e-4 * 50))
  expect_equal(p$offset, 1.02 * f$W[, 1])
  expect_equal(p$A0, cbind(f$gmean, f$alpha[, -1]))
  expect_equal(sum(.spaNormPolishProblem(f, "fit")$cells_idx), 25)
})

test_that("a per-gene library-size coefficient is refused", {
  expect_error(.spaNormPolishProblem(mk_fit(shared = FALSE)), "not .*shared")
})

test_that("a fit saved before the polish slot existed reads as unpolished", {
  f <- mk_fit()
  expect_false(isPolished(f))
  # simulate an old object: drop the slot from the attributes
  old <- f; attr(old, "polish") <- NULL
  expect_false(isPolished(old))
  expect_identical(.polishSlot(old), list())
})
```

If R refuses to build a `SpaNormFit` without the `polish` slot, construct the "old" object with `unclass`-free attribute removal as shown. If `attr<-` restores the prototype, use a real RDS written by SpaNorm 1.7.12, committed as `tests/testthat/fixtures/spanormfit_1712.rds` and generated from `mk_fit()` under the installed 1.7.12.

- [ ] **Step 2: Run to fail**, **Step 3: implement**, **Step 4: run to pass.**

Run: `Rscript -e 'devtools::document(); devtools::test(filter = "polishSpaNormModel|SpaNormFit")'`
Expected: PASS. The existing `test-SpaNormFit.R` must still pass: the slot has a prototype, so `SpaNormFit(...)` calls without `polish` keep working.

- [ ] **Step 5: Commit** (`R/AllClasses.R`, `R/polishSpaNormModel.R`, the tests, `man/`).

---

### Task 7: SpaNorm, `polishSpaNorm()` with `ls = "fixed"`

**Files:**
- Create: `SpaNorm/R/polishSpaNorm.R`
- Modify: `SpaNorm/R/AllGenerics.R` if it exists; otherwise put the `setGeneric` at the top of `polishSpaNorm.R`, as `mainSpaNorm.R` does.
- Create: `SpaNorm/tests/testthat/helper-polish.R` (a cached small fit) and `test-polishSpaNorm.R`

**Interfaces:**
- Consumes: `polishNB` (T4), `.spaNormPolishProblem`, `isPolished`, `.polishSlot` (T6); SpaNorm's `getSpaNormFit`, `getAdjustmentFun`, `normaliseBlocked`.
- Produces:

```r
setGeneric("polishSpaNorm", function(spe, adj.method = c("auto", "logpac", "pearson", "medbio", "meanbio"),
  scale.factor = 1, psi.method = c("fixed", "profile"), ls = c("fixed", "joint"),
  cells = c("all", "fit"), null = TRUE, maxit = 50L, tol = 1e-8,
  engine = c("batch", "gene"), batch.size = NULL, backend = c("cpu", "auto", "gpu"),
  BPPARAM = BiocParallel::SerialParam(), verbose = TRUE, assay = NULL, ...)
  standardGeneric("polishSpaNorm"))
# methods: SpatialExperiment, Seurat.
.polishSpaNormFit(fit, Y, psi.method, ls, cells, ...)   # SpaNormFit -> polished SpaNormFit
.polishSharedLS(Y, prob, pol, psi.method, ...)          # T9's hook; errors in T7
```
  After the call: `metadata(spe)$SpaNorm` is polished, `metadata(spe)$SpaNormUnpolished` is the input, and `logcounts` is rewritten. If `null = TRUE` and `SpaNormNull` exists, the null is polished too, the same way (`SpaNormNullUnpolished` is kept). The SVG columns in `rowData` are dropped with a warning.

- [ ] **Step 1: Write the fixture helper and the failing tests.**

```r
# tests/testthat/helper-polish.R
# One small real SpaNorm fit, built once per test run (~15 s).
.polish_spe <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      data(HumanDLPFC, package = "SpaNorm", envir = environment())
      set.seed(20260925)
      top <- order(-Matrix::rowSums(SummarizedExperiment::assay(HumanDLPFC, "counts")))[1:40]
      spe <- HumanDLPFC[top, 1:600]
      cache <<- SpaNorm(spe, df.tps = 2L, sample.p = 0.25, verbose = FALSE,
                        backend = "cpu")
    }
    cache
  }
})
```

```r
# tests/testthat/test-polishSpaNorm.R
test_that("the polished fit is at each gene's optimum with a1 held", {
  spe <- polishSpaNorm(.polish_spe(), verbose = FALSE)
  f <- S4Vectors::metadata(spe)$SpaNorm
  f0 <- S4Vectors::metadata(spe)$SpaNormUnpolished
  expect_true(isPolished(f)); expect_false(isPolished(f0))
  expect_equal(f$alpha[, 1], f0$alpha[, 1])                  # a1 untouched
  p <- .spaNormPolishProblem(f)
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  A <- cbind(f$gmean, f$alpha[, -1])
  for (g in seq_len(nrow(Y))) {
    mu <- pmax(exp(as.numeric(p$X %*% A[g, ] + p$offset)), nbMuFloor())
    sc <- crossprod(p$X, (Y[g, ] - mu) / (1 + f$psi[g] * mu)) - p$pen * A[g, ]
    expect_lt(max(abs(sc)) / sum(Y[g, ]), 1e-6, label = rownames(Y)[g])
  }
  expect_equal(f$psi, f0$psi)                                # psi.method = "fixed"
})

test_that("polishing raises every gene's penalised likelihood", {
  spe <- .polish_spe()
  f0 <- S4Vectors::metadata(spe)$SpaNorm
  f <- S4Vectors::metadata(polishSpaNorm(spe, verbose = FALSE))$SpaNorm
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  pl <- function(fit) {
    p <- .spaNormPolishProblem(fit)
    A <- cbind(fit$gmean, fit$alpha[, -1])
    vapply(seq_len(nrow(Y)), function(g) {
      mu <- pmax(exp(as.numeric(p$X %*% A[g, ] + p$offset)), nbMuFloor())
      sum(dnbinom(Y[g, ], size = 1 / fit$psi[g], mu = mu, log = TRUE)) -
        0.5 * sum(p$pen * A[g, ]^2)
    }, numeric(1))
  }
  expect_true(all(pl(f) >= pl(f0) - 1e-8 * abs(pl(f0))))
})

test_that("logcounts are rewritten from the polished fit", {
  spe <- .polish_spe()
  out <- polishSpaNorm(spe, verbose = FALSE)
  expect_false(identical(SummarizedExperiment::assay(out, "logcounts"),
                         SummarizedExperiment::assay(spe, "logcounts")))
  # the same result as SpaNorm() re-normalising with the polished fit cached
  again <- suppressWarnings(SpaNorm(out, df.tps = 2L, sample.p = 0.25, verbose = FALSE))
  expect_equal(SummarizedExperiment::assay(again, "logcounts"),
               SummarizedExperiment::assay(out, "logcounts"))
})

test_that("a non-integer assay is refused at the entry point", {
  spe <- .polish_spe()
  SummarizedExperiment::assay(spe, "counts") <- SummarizedExperiment::assay(spe, "counts") + 0.5
  expect_error(polishSpaNorm(spe, verbose = FALSE), "integer counts")
})

test_that("an all-zero gene does not abort the polish", {
  spe <- .polish_spe()
  cnt <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  cnt[1, ] <- 0
  SummarizedExperiment::assay(spe, "counts") <- cnt
  out <- polishSpaNorm(spe, verbose = FALSE)
  pol <- S4Vectors::metadata(out)$SpaNorm@polish$genes
  expect_true(all(pol$polished[-1]))
  expect_true(all(is.finite(S4Vectors::metadata(out)$SpaNorm$alpha[-1, ])))
})

test_that("ls = 'joint' is not available until Task 9", {
  expect_error(polishSpaNorm(.polish_spe(), ls = "joint", verbose = FALSE),
               "not implemented")
})
```

The all-zero test edits the counts after the fit, which is valid: the polish reads the counts, and the fit's `gmean` for that gene is just a bad start. Delete the last test in Task 9.

- [ ] **Step 2: Run to fail.** Expected: `could not find function "polishSpaNorm"`.

- [ ] **Step 3: Implement.**

```r
.polishSpaNormFit <- function(fit, Y, psi.method = "fixed", ls = "fixed",
                              cells = "all", ...) {
  prob <- .spaNormPolishProblem(fit, cells)
  Yc <- Y[, prob$cells_idx, drop = FALSE]
  pol <- polishNB(Yc, prob$X, prob$A0, prob$psi, lambda.a = prob$pen,
                  offset = prob$offset, start.cols = c(TRUE, rep(FALSE, ncol(prob$X) - 1)),
                  psi.method = psi.method, ...)
  a1 <- prob$a1
  if (ls == "joint") {
    j <- .polishSharedLS(Yc, prob, pol, psi.method = psi.method, ...)
    pol <- j$pol; a1 <- j$a1
  }
  fit@gmean <- as.numeric(pol$alpha[, 1])
  fit@alpha <- cbind(a1, pol$alpha[, -1, drop = FALSE])
  colnames(fit@alpha) <- colnames(fit$W)
  fit@psi <- as.numeric(pol$psi)
  fit@polish <- list(
    settings = list(psi.method = psi.method, ls = ls, cells = cells,
                    pen = prob$pen, a1.input = prob$a1, a1 = a1,
                    SpaNorm = as.character(utils::packageVersion("SpaNorm"))),
    genes = pol$polish)
  fit
}

.polishSharedLS <- function(Y, prob, pol, ...) {
  stop("ls = \"joint\" is not implemented yet", call. = FALSE)
}
```

  - `start.cols = c(TRUE, FALSE, ...)` makes the sane start put the (offset-corrected) overall log mean on the `gmean` column.
  - `fit@loglik` is left as the fit's iteration trace. The per-gene log-likelihood lives in `@polish$genes`.
  - The SPE method: `checkSPE`, fetch the fit with `getSpaNormFit(spe)` (it errors if absent), polish it, polish `getSpaNormFit(spe, null = TRUE, validate = FALSE)` if `null` and present, then normalise with `normaliseBlocked(getAdjustmentFun(fit$gene.model, adj.method), emat, scale.factor, fit, BPPARAM, block.size)`. Use the same `block.size` rule as `.spaNormCore` and put it in a shared helper instead of copying it.
  - Write the metadata (`SpaNorm`, `SpaNormUnpolished`, `SpaNormNull`, `SpaNormNullUnpolished`) and drop the SVG columns: `getSVGResults(spe, stop_if_missing = FALSE)` names them, as in `SpaNormSVG`.
  - The Seurat method mirrors it through `spe@misc` and `SetAssayData`, like the `SpaNorm` Seurat method.
  - Refuse to re-polish an already polished fit unless `overwrite = TRUE` is passed via `...`. Otherwise it would polish from the polished fit, and `SpaNormUnpolished` would be overwritten by a polished fit.

- [ ] **Step 4: Run.**

Run: `Rscript -e 'devtools::document(); devtools::test()'`
Expected: PASS.

- [ ] **Step 5: Review, then commit.** `numerical-robustness-reviewer` on `R/polishSpaNorm.R`.

---

### Task 8: SpaNorm, polished SVG consistency

**Files:**
- Modify: `SpaNorm/R/mainSpaNormSVG.R` (`SpaNormSVG` method, `svgTest`)
- Modify: `SpaNorm/tests/testthat/test-mainSpaNormSVG.R`

**Interfaces:**
- Consumes: `isPolished`, `.polishSlot`, `.polishSpaNormFit` (T6/T7); `fitSpaNormTechnical` (1.7.12, scaled penalty).
- Produces: `svgTest(Y, fit.spanorm, fit.technical)` stops with `"the full and null SpaNorm fits must be polished alike"` when `isPolished()` differs. `SpaNormSVG()` fits a missing null and polishes it with `fit.spanorm@polish$settings` (psi.method, ls, cells).

- [ ] **Step 1: Write the failing tests.**

```r
test_that("SpaNormSVG polishes the null when the full fit is polished", {
  spe <- polishSpaNorm(.polish_spe(), null = FALSE, verbose = FALSE)
  out <- SpaNormSVG(spe, verbose = FALSE)
  nul <- S4Vectors::metadata(out)$SpaNormNull
  expect_true(isPolished(nul))
  expect_identical(nul@polish$settings$psi.method,
                   S4Vectors::metadata(out)$SpaNorm@polish$settings$psi.method)
})

test_that("svgTest refuses a mixed pair", {
  spe <- SpaNormSVG(.polish_spe(), verbose = FALSE)          # unpolished pair
  full <- S4Vectors::metadata(spe)$SpaNorm
  nul <- S4Vectors::metadata(spe)$SpaNormNull
  Y <- SummarizedExperiment::assay(spe, "counts")
  fullP <- .polishSpaNormFit(full, as.matrix(Y))
  expect_error(svgTest(Y, fullP, nul), "polished alike")
})

test_that("polishSpaNorm polishes an existing null and clears stale SVG results", {
  spe <- SpaNormSVG(.polish_spe(), verbose = FALSE)
  expect_warning(out <- polishSpaNorm(spe, verbose = FALSE), "SVG")
  expect_true(isPolished(S4Vectors::metadata(out)$SpaNormNull))
  expect_false("svg.F" %in% colnames(SummarizedExperiment::rowData(out)))
})

test_that("the polished LRT statistic is non-negative for every gene", {
  out <- SpaNormSVG(polishSpaNorm(.polish_spe(), verbose = FALSE), verbose = FALSE)
  rd <- SummarizedExperiment::rowData(out)
  expect_true(all(rd$svg.F >= 0))
})
```

  The last test is the nesting check. The polished full fit maximises a strictly larger model than the polished null, so the unclamped `2(ll_full - ll_null)` should be ≥ 0 up to tolerance. To check it before the clamp, have the implementation record `attr(df.svg, "F.raw")`, and assert `min(F.raw) > -1e-6 * max(abs(F.raw))`. A clearly negative raw F means the pair was polished inconsistently.

  **Caveat, record it and do not paper over it:** `svgTest` scores with `winsorisePsi(psi)`, not the polish's `psi`. Under `psi.method = "fixed"` the two share the same fit `psi` per model, but they differ between the full and null models, so small negative raw F values are possible and are not a polish bug. If the check fails, measure how many genes and by how much before deciding.

- [ ] **Step 2: Run to fail.** **Step 3: Implement** the three branches in `SpaNormSVG`:
  - null absent: fit it, and if the full fit is polished, polish the null with the full fit's settings;
  - null present and pairing consistent: use it;
  - null present and pairing inconsistent: if the full fit is polished and the null is not, polish the null with a message; otherwise error.

  Add the guard to `svgTest`. **Step 4: Run** `devtools::test()`. **Step 5: Commit.**

---

### Task 9: SpaNorm, `ls = "joint"`: the pooled profiled-Newton step on `a1`

**Files:**
- Create: `SpaNorm/R/polishSpaNormJoint.R` (it replaces the `.polishSharedLS` stub; delete the stub from `polishSpaNorm.R`)
- Create: `SpaNorm/tests/testthat/test-polishSpaNormJoint.R`

**Interfaces:**
- Consumes: `polishNB`, `nbNewtonSolver`, `nbMuFloor`, `.spaNormPolishProblem`.
- Produces: `.polishSharedLS(Y, prob, pol, psi.method, maxit.ls = 10L, tol.ls = 1e-6, block = 500L, ...)` → `list(pol, a1, iterations, score)`.

**The step** (spec §4). With `D_g = mu/(1 + psi mu)` and `r_g = (y - mu)/(1 + psi mu)` at the current polished per-gene fit:
`U = sum_g w1' r_g`,
`I = sum_g [ (w1^2)' D_g - c_g' H_g^{-1} c_g ]` with `c_g = X' (D_g * w1)` and `H_g = X' D_g X + diag(pen)`.
`I` is the information for `a1` after profiling out each gene's own coefficients, and `solver$solve(solver$factor(D_g), c_g)` gives `H_g^{-1} c_g` from the same factorisation the polish uses.

- [ ] **Step 1: Write the failing tests, including an `optim()` oracle.**

```r
# tests/testthat/test-polishSpaNormJoint.R
test_that("joint a1 matches a full optim() of the joint objective", {
  set.seed(30)
  n <- 120; G <- 5
  w1 <- log(runif(n, 0.5, 2)); x2 <- rnorm(n)
  X <- cbind(1, x2); pen <- c(0, 2)
  Y <- t(vapply(seq_len(G), function(g)
    rnbinom(n, mu = exp(0.3 + g / 5 + 0.8 * w1 + 0.3 * x2), size = 5), numeric(n)))
  psi <- rep(0.2, G)
  negobj <- function(par) {
    a1 <- par[1]; A <- matrix(par[-1], G, 2)
    -sum(vapply(seq_len(G), function(g) {
      mu <- exp(as.numeric(X %*% A[g, ] + a1 * w1))
      sum(dnbinom(Y[g, ], size = 1 / psi[g], mu = mu, log = TRUE)) - 0.5 * sum(pen * A[g, ]^2)
    }, 0))
  }
  opt <- optim(c(1, rep(0, 2 * G)), negobj, method = "BFGS",
               control = list(maxit = 5000, reltol = 1e-14))
  prob <- list(X = X, pen = pen, a1 = 1, offset = 1 * w1, w1 = w1,
               A0 = matrix(0, G, 2), psi = psi, cells_idx = rep(TRUE, n))
  pol <- polishNB(Y, X, prob$A0, psi, lambda.a = pen, offset = prob$offset,
                  psi.method = "fixed")
  j <- .polishSharedLS(Y, prob, pol, psi.method = "fixed")
  expect_equal(j$a1, opt$par[1], tolerance = 1e-4)
  expect_equal(j$pol$alpha, matrix(opt$par[-1], G, 2), tolerance = 1e-3)
  expect_lt(abs(j$score), 1e-6 * sum(Y))
})

test_that("joint polish never lowers the total penalised likelihood", {
  spe <- .polish_spe()
  fx <- S4Vectors::metadata(polishSpaNorm(spe, ls = "fixed", verbose = FALSE))$SpaNorm
  jt <- S4Vectors::metadata(polishSpaNorm(spe, ls = "joint", verbose = FALSE))$SpaNorm
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  tot <- function(fit) {
    p <- .spaNormPolishProblem(fit)
    A <- cbind(fit$gmean, fit$alpha[, -1])
    sum(vapply(seq_len(nrow(Y)), function(g) {
      mu <- pmax(exp(as.numeric(p$X %*% A[g, ] + p$offset)), nbMuFloor())
      sum(dnbinom(Y[g, ], size = 1 / fit$psi[g], mu = mu, log = TRUE)) -
        0.5 * sum(p$pen * A[g, ]^2)
    }, 0))
  }
  expect_gte(tot(jt), tot(fx) - 1e-8 * abs(tot(fx)))
  expect_true(length(unique(jt$alpha[, 1])) == 1)        # still shared
})
```

- [ ] **Step 2: Run to fail.** Expected: the stub's `"not implemented"` error.

- [ ] **Step 3: Implement.**

```r
.polishSharedLS <- function(Y, prob, pol, psi.method = "fixed", maxit.ls = 10L,
                            tol.ls = 1e-6, block = 500L, ...) {
  X <- prob$X; w1 <- prob$w1; pen <- prob$pen
  solver <- nbNewtonSolver(X, pen)
  G <- nrow(Y)
  blocks <- split(seq_len(G), ceiling(seq_len(G) / block))
  total <- function(A, psi, a1) {
    sum(vapply(blocks, function(b) {
      Mu <- pmax(exp(A[b, , drop = FALSE] %*% t(X) +
                     matrix(a1 * w1, length(b), length(w1), byrow = TRUE)), nbMuFloor())
      sum(stats::dnbinom(Y[b, , drop = FALSE], size = 1 / psi[b], mu = Mu, log = TRUE)) -
        0.5 * sum(sweep(A[b, , drop = FALSE]^2, 2, pen, `*`))
    }, 0))
  }
  score_info <- function(A, psi, a1) {
    U <- 0; I <- 0
    for (b in blocks) {
      Mu <- pmax(exp(A[b, , drop = FALSE] %*% t(X) +
                     matrix(a1 * w1, length(b), length(w1), byrow = TRUE)), nbMuFloor())
      den <- 1 + psi[b] * Mu                      # psi[b] recycles down columns
      R <- (Y[b, , drop = FALSE] - Mu) / den
      D <- Mu / den
      U <- U + sum(R %*% w1)
      C <- (D * matrix(w1, nrow(D), length(w1), byrow = TRUE)) %*% X
      for (k in seq_along(b)) {
        h <- solver$solve(solver$factor(D[k, ]), C[k, ])
        I <- I + sum(D[k, ] * w1^2) - (if (is.null(h)) 0 else sum(C[k, ] * h))
      }
    }
    list(U = U, I = I)
  }
  a1 <- prob$a1
  A <- pol$alpha; psi <- pol$psi
  ll <- total(A, psi, a1)
  it <- 0L
  si <- score_info(A, psi, a1)
  while (it < maxit.ls && abs(si$U) / sqrt(max(si$I, .Machine$double.eps)) > tol.ls) {
    it <- it + 1L
    if (!is.finite(si$I) || si$I <= 0) {
      stop("the profiled information for the shared library-size coefficient ",
           "is not positive (", format(si$I), "); the joint polish cannot step",
           call. = FALSE)
    }
    step <- si$U / si$I
    accepted <- FALSE
    for (h in 0:10) {
      a1n <- a1 + step / 2^h
      pn <- polishNB(Y, X, A, psi, lambda.a = pen, offset = a1n * w1,
                     psi.method = psi.method, warm = TRUE, ...)
      lln <- total(pn$alpha, pn$psi, a1n)
      if (is.finite(lln) && lln >= ll - 1e-9 * abs(ll)) { accepted <- TRUE; break }
    }
    if (!accepted) break
    a1 <- a1n; pol <- pn; A <- pn$alpha; psi <- pn$psi; ll <- lln
    si <- score_info(A, psi, a1)
  }
  list(pol = pol, a1 = a1, iterations = it, score = si$U)
}
```

  - Under `psi.method = "profile"`, `psi` moves with each re-polish. The objective is then the profile likelihood, and the line search still guarantees ascent in it.
  - Record `iterations` and `score` in `@polish$settings` (`ls.iterations`, `ls.score`).
  - The per-gene `solver$factor` inside `score_info` is the dense path: SpaNorm has no absorbed block. It costs one `p x p` factorisation per gene per outer iteration, which is small next to the re-polish.

- [ ] **Step 4: Run.**

Run: `Rscript -e 'devtools::test(filter = "polishSpaNorm")'`
Expected: PASS. Delete Task 7's "not implemented" test.

- [ ] **Step 5: Review, then commit.** Run `numerical-robustness-reviewer`. The one non-guarded inverse here goes through `solver$solve`, which returns `NULL` when singular; the code treats that as a zero correction. The reviewer must judge whether that is conservative enough or whether the gene should be excluded from `U` as well.

---

### Task 10: Measurement on YTMA cores (HPC; you decide the default)

**Files:**
- Create: `YTMACosMxWTAv2/claude/code/80_polishSpaNorm_validate.R` and its sbatch (in the YTMA repo, following `31_polished_svg.R`).

- [ ] **Step 1: Reproduce `_polish.R`.** On cores 709_51 (395 cells) and 30 (1,161 cells), run `polish_fit(Y, full, nc, scale_pen = TRUE, psi_method = "moderated")` against `.polishSpaNormFit(full, Y, psi.method = "fixed", ls = "fixed", cells = "all")`. Expect `gmean`/`alpha` equal to `1e-6` (spec acceptance 3). A difference means `_polish.R`'s convergence rule (`gain < tol * |ll|`, no restart) stopped earlier. Report the per-gene score norm of both before choosing which is right.
- [ ] **Step 2: Compare `ls = "fixed"` and `"joint"` on four cores** (two small, two large). Measure:
  - `a1` input vs joint;
  - `cor(logcounts)`;
  - residual Moran's I and `dev_sd`, the v2 report's metrics;
  - polished SVG counts at FDR 0.05;
  - wall time.

  Seed and record: `set.seed(20260925)`, the SpaNorm commit and the job ids.
- [ ] **Step 3: Decide the default and record it.** Use the `record-finding` skill to update the v2 report, SpaNorm `NEWS.md` and the memory file. If joint moves `a1` by less than its own profiled SE, `sqrt(1/I)`, on every core, keep `"fixed"`: joint then only adds cost. If it moves more, show the v2 metrics before choosing.

The submission follows the `hpc-job-submission`/`hpc-job-sizing` skills. Give an explicit `--export=ALL,R_LIBS=$SCRATCH/Rlib_polishmove` (the memory note: sbatch drops the environment inside an allocation), and copy the driver to a task-private file at task start (CLAUDE.md: the incremental-parse hazard).

---

### Task 11: Docs, release, downstream

- [ ] **Step 1: SpaNorm.**
  - Add a vignette section "Polishing the fit" covering when to polish, the SVG pairing rule and the `ls` choice, with numbers from Task 10.
  - Add `NEWS.md` entries for 1.7.12–1.7.14.
  - Add `polishSpaNorm` and `polishNB` to `_pkgdown.yml`.
  - Run `BiocCheck`.
- [ ] **Step 2: spiDE.**
  - Update CLAUDE.md: the pipeline stage 2b paragraph ("the per-gene engine lives in SpaNorm (`polishNB`, 1.7.13); spiDE keeps the tau2 loop"), the "Key invariant" paragraph's list of batching helpers, and the `R/inference-batch.R` references.
  - Add a `NEWS.md` entry.
  - Run `evidence-auditor` on the changed CLAUDE.md paragraphs and the vignette.
- [ ] **Step 3: YTMA.** Replace `claude/code/_polish.R` with calls to `SpaNorm::polishSpaNorm`, and keep the old file under `_superseded/`. Do not rerun completed arms; the measurement in Task 10 covers the equivalence.
- [ ] **Step 4: Release.** The integrator bumps SpaNorm to 1.7.14 and spiDE's minimum, rebuilds both pkgdown sites locally before pushing (memory: pkgdown failure modes), and pushes.

---

## Self-review notes

- **Spec coverage:**
  - §2 (move/stay): T2–T5.
  - §3 (API): T2, T4, T7, T8.
  - §4 (model, `ls`, psi and cells defaults): T6, T7, T9, T10.
  - §5 (SVG): T8.
  - §6 (slot): T6.
  - §7 acceptance 1: T1 + T5.
  - Acceptance 2: T3/T4.
  - Acceptance 3: T10 step 1.
  - Acceptance 4: T9.
  - Acceptance 5: T8.
  - Acceptance 6: T5 step 5 and T11.
- **Names used across tasks:** `polishNB`, `nbNewtonSolver`, `nbNewtonSolverBatch`, `nbGramBatch`, `nbAbsorbGramBatch`, `nbMuFloor`, `isPolished`, `.polishSlot`, `.spaNormPolishProblem` (fields `X, pen, a1, offset, w1, A0, psi, cells_idx`), `.polishSpaNormFit`, `.polishSharedLS`, `.testedStartCols`. Each is defined in exactly one task's Interfaces block.
- **Known soft spots, called out in the tasks rather than hidden:**
  - the golden test's generator arguments must be checked against the merged `main` (T1);
  - the old-object fixture may need a real 1.7.12 RDS (T6);
  - the raw-F nesting check can be disturbed by `winsorisePsi` (T8);
  - the singular-gene handling in the joint score (T9).
