# Fitting stage v2: a batched converged fitter: design

**Status:** approved in chat 2026-09-15 (full replacement, fresh recalibration,
both statistical defects in scope). Phases 0 and 0b landed on
`feature/batched-fitter`; Phase 1 follows.

## Why

The fit is two stages with a bad cost split. `fitSpiDE()` calls
`SpaNorm::fitNB`, which shares one gene-averaged cell weight vector across
genes, decides convergence on the *aggregate* log-likelihood and clamps
coefficient columns across genes — a different estimator, not an
under-converged per-gene one. `polishSpiDE()` then converges each gene on its
own penalised likelihood, which is what every downstream statistic is built on.

On the v11 ICI arm (13,348 genes, 77,454 cells, 1,107 columns) the polish is
430–830 min on 64 cores against 20 min for the fit and 37 for the inference:
**>98% of the pipeline**, 2,800–5,400 core-hours per six-grid arm. That cost is
why calibration stops at one bandwidth and five seeds; the campaign the cohort
reports ask for (20 seeds × 4 bandwidths) is 37,000–71,000 core-hours.

## What the measurement says (2026-09-15, job 28494601)

`notes/polish_cost/microbench.R` at the true shape (n = 77,454, px = 398 dense
after the SampleCellTypeInt block G = 709 is absorbed):

| | 1 thread | 8 threads |
|---|---|---|
| cached factorisation vs rebuild | **3.0×** | 2.9× |
| `dsyrk` vs `dgemm` for the gram | 1.7× | **16.4×** |
| batched matvec, B = 256 | 36.8× | **135×** |

Three consequences, and they set the order of work:

1. **The information matrix was rebuilt every step**, though `R/polish.R:141-143`
   documents reuse for up to three. A reuse costs 1% of a rebuild. Landed in
   `3d78d1c`, bit-identical.
2. **The gram form was the wrong one.** `crossprod(X, X*w)` *degrades* with
   threads (0.43 → 2.44 s) while `crossprod(X*sqrt(w))` scales (0.25 → 0.15 s).
   This is the mechanism behind the recorded "one worker gains only 1.5× from
   four BLAS threads", so the 64 × 1 worker/thread policy is a workaround for a
   fixable problem. Landed in `e10a7cc`.
3. **Batching amortises the design read**, 90 → 0.8 ms/gene. This is the
   premise of Phase 1 and it holds.

It also **cancelled** a planned task: one `dnbinom` sum is 0.007 s against
0.46 s of gram, so lifting the ψ-only `lgamma` constants out of the line search
is worth ~1% of an iteration and is not worth the numerical risk.

## Scope

In: a block-batched per-gene Newton as the primary fitter; a device path; one
variance-component loop; `fitNB` demoted to starting values and the cross-gene
dispersion moderation; the Satterthwaite between-sample reference df.

Out: the dispersion *rule* (`psi = "moderated"` keeping `fitNB`'s cross-gene
value). Phase 4 makes that value degrade, so the rule errors there rather than
silently changing meaning; replacing it by squeezing the converged per-gene
profile dispersions is a separate change with its own gate.

## The estimator does not change

`.polishGene()` already computes the per-gene penalised MLE, and
`.newtonSolver()` already absorbs the nested block analytically (Schur
complement, diagonal `C`). Phase 1 is a restructuring, not a new estimator: the
same objective, the same convergence rule, the same restart and fallback
semantics. `.polishGene()` stays in the tree as the reference implementation
and the test oracle, reachable through `engine = "gene"`.

## Order of work

| phase | content | gate |
|---|---|---|
| 0 | measure at the true shape | recorded in FINDINGS ✅ |
| 0b | memoise the factorisation; symmetric gram | bit-identical / objective not worse ✅ |
| 1 | `.polishBatch()`, batched on CPU | batched == unbatched to 1e-10; blocking invariance |
| 2 | device path via `.gramBatch()`'s torch branch | CPU == GPU at `gpu_tol()` |
| 2a-2d | the inference side of the device path | landed, see below |
| 3+4 | one τ² loop; `fitNB` to starting values | the `longtests` τ² window and df anchors |
| 5 | the between-sample reference df | the `lmerTest` and `S - 2` anchors ✅ |
| 6 | revalidation and recalibration | benchmark + null grids ✅ |

Phase 3 is a verification rather than a phase: `.tau2Iterate()` already
iterates to a tolerance with Steffensen acceleration (`polish.R:534-611`). What
remains capped is the *fit's* loop (`re.maxit = 2L`), which dissolves when
`fitNB` is demoted.

## Phase 2, as built (2026-09-16)

The inference half landed first because it needed no new numerics: the budget
fix (`9a7e93d`), `.absorbBatch()` (`14652b5`), the guard flip in
`.blockedInference()` (`6c5bdfc`) and the cell tile on `.gramBatch()`
(`f9d6630`). The device path now absorbs the nested block instead of inverting
a dense `p x p` gram, and both batched grams are bounded by a cell tile rather
than by `ncells`.

**The hard part of the remaining half is not the gram.** `.gramBatch()` and
`.absorbBatch()` are already tensor-capable, so `.polishBatch()`'s information
matrix is solved. What is not is its **factorisation state**: `newton()` keeps
`fac`, a *list of per-gene factorisations*, refreshed per gene under a per-gene
staleness counter, and calls `solver$solve(fac[[k]], S[k, ])` one gene at a
time. That list is the thing that cannot go to a device -- it is R objects, one
per gene, and the per-gene `solve()` is a kernel launch per gene per iteration.

So Phase 2e is not "add a backend argument to `.polishBatch()`". It is:

1. Replace `.newtonSolver()`'s per-gene closures, *inside the batched engine
   only*, with a batched factor/solve: `.absorbBatch()` for the stack, one
   batched Cholesky for the factorisation, one batched triangular solve for the
   step. `.polishGene()` keeps the closures; it is the reference implementation.
2. ~~Decide what staleness means when the factorisation is one tensor.~~
   **Measured, and the prediction here was wrong** (FINDINGS, 2026-09-16, job
   28511603). This paragraph said a shared tensor factorisation "throws away
   most of the memoisation the 3x in Phase 0b came from" and that it had to be
   settled before any code. It costs **11% more factorisations at a 128-gene
   batch and 6% at 64**, on 400 genes of the top-5 bandwidth-30 checkpoint.
   Refreshes cluster: genes in a block share the design and converge in similar
   numbers of iterations, so the iterations where any gene is stale are largely
   the iterations where most are. **Carry one shared factorisation for the
   active stack** -- it is the simpler implementation and the cheaper one.
3. The line search and the profile-psi bisection are already batch-shaped and
   are elementwise over genes x cells, which is what the device wants.

**Both halves are now measured** (FINDINGS, 2026-09-16, jobs 28511603 and
28527447), and they split the design cleanly:

- *Numerically* the shared factorisation is free. On 300 real genes at batch 64
  no gene's objective is worse (worst relative -6.5e-13), the convergence flags
  are identical gene for gene, alpha agrees to 8e-6 and the iteration count
  falls slightly, because a fresher information matrix is a marginally better
  Newton direction.
- *On the CPU it is 1.7x slower* (3.3 to 5.5 min). The shared state buys kernel
  launches, and a CPU has none to buy: what is left is a per-slice Cholesky and
  a per-gene right-hand side against a cached LU under Phase 0b's memoisation.

So **`shared.factor` is device-only**. `engine = "batch"` on the CPU keeps the
per-gene list; the device path takes the shared stack; and the two are the same
estimator reaching the same optimum by different paths -- measured, not
asserted. The remaining Phase 2e work is the rest of the loop (mu, the score,
the line search's `dnbinom` column-sums, the profile-psi bisection) on tensors,
with `shared.factor` switched on by the backend rather than by the caller.

What point 2 does NOT settle is the numerical path. A synchronous refresh gives
genes that were not stale a fresher information matrix than they would have
had. That cannot move the fixed point and a fresher matrix is not a worse
Newton direction, but the batch/per-gene comparison stops being an equality
test -- the same trade already recorded for the line search, and it is measured
on the objective, not assumed. The ratio also has not been checked on the
unrestricted arm (1,107 columns, 764 nested against 662/492 here).

## Phases 5 and 6, as built (2026-09-18/21)

**Phase 5 found an inversion, not a scale error.** The Satterthwaite df for a
between-patient contrast was not merely too large, it was anticorrelated with
the evidence: Spearman cor(df, cells) = **-0.978** over the 13 compartments,
Tumor (27,764 cells) getting df 306 where Smooth muscle (1,073) got 34,544. A
rare compartment's (patient x cell type) groups hold few cells, their random
effects shrink toward the pooled ones, the variance-component gradients
collapse, `varse2` shrinks, and `df = 2 vjj^2 / varse2` runs away toward the
cell scale -- anti-conservative exactly where the data is thinnest.

The fix is `.boundPatientDF()` with `.patientsPerTested()`: no between-patient
contrast may exceed that compartment's OWN patient count minus two. It is
applied at the fit site AND the polish site, because the polish recomputes
every quantity the fit computed and a correction applied at one site only is
inert where it matters (1d046e9).

**Which columns it binds was wrong at first, and the grids caught it.** The
first version asked whether a column is CONSTANT inside every (patient x cell
type) group -- true of `CellType_k:Responder`, false of
`CellType_k:Responder:Niche_n` because niche density varies cell to cell inside
a group. That question is about how the REGRESSOR varies, not about what
replicates the CONTRAST. Responder status is a property of the patient either
way, so the three-way term compares groups of patients exactly as the two-way
one does. The `_dfb` grids showed the consequence directly: cell-type layer
calibrated (FDP 0.288 -> 0.280), niche layer still anti-conservative. The bound
now covers every tested column carrying the condition factor (8d2f4c9).

**That made `satterthwaite` the wrong default.** Once every tested column is
condition-bearing and every one of them is capped, the per-column vector IS the
bound -- it binds on 168 of 168 columns on the cohort design -- so it carries
nothing the scalar does not, while costing a variance-parameter covariance per
gene and inviting callers to read structure into a constant. `df.method`
defaults to `"between"`; `"satterthwaite"` remains available and is bounded the
same way, so that arm cannot mislead either.

Three tests were re-anchored rather than relaxed, and the re-anchors assert
mechanism so they cannot pass vacuously: `test-polish-stage.R` used to require
that the polish CHANGES `@df`, which it no longer does because both sites land
on the same cap, so it now poisons the input `@df` and requires the polish to
reproduce the bound -- which a copy could not; `test-mixedEffects.R` moved its
per-column claim onto an explicit `satterthwaite` fit; the longtest df anchors
moved to the bound.

**Phase 6 recalibrated on two null axes, not one.** The harness driver gained a
`SPIDE_SHUF=perm` mode (a patient-level permutation of the condition), so the
block shuffle (niche association destroyed, condition kept) and the response
permutation (condition destroyed, niche kept) come out of one chain. Reported
side by side, because a layer whose FDP differs sharply between them is riding
on the axis the other null preserved. On the 22 x 8 fine arm: block FDP
0.300 / 0.254 / 0.273 by layer against permutation 0.672 / 0.702 / 0.673. The
gap is a property of a 55-patient 25/30 cohort -- v11 already measured 0.624
and 0.835 -- not of the estimator.

## NEWS owed at merge

The branch does not touch `NEWS.md` or `DESCRIPTION` -- Phase 1 shipped a new
default engine without them, deliberately, since the version this lands under
is a release decision and every edit there is merge-conflict surface. So the
entries are recorded here instead, to be written once:

- **`polishSpiDE(engine = "batch")` is the default** and converges a block of
  genes together; `"gene"` is the per-gene reference implementation. 4.2x on
  the production arm, 6.1x on the top-5 arm, alpha to 1.6e-13 and psi bit-
  identical. `batch.size` and `options(spiDE.polish.mem.budget)` bound it.
- **`spiDE()` forwards `engine` and `batch.size`**, so the reference engine is
  reachable from the top-level entry point.
- **The nested block is absorbed on the GPU too**, not only the CPU: the
  covariance sub-batch inverts px x px rather than p x p (398 against 1,107 on
  the cohort design).
- **The covariance memory budget is divided among workers.** Behaviour change
  for anyone running a wide design under many workers: sub-batches get smaller
  and peak memory stops scaling with worker count.
- Phase 0b, if it is mentioned at all, is 4% end to end and NOT 3x -- the
  microbenchmark's 3x does not survive contention (FINDINGS, 2026-09-15).
- **`df.method` now defaults to `"between"`.** USER-VISIBLE BEHAVIOUR CHANGE:
  `@df` comes back a scalar where it used to be a named per-tested-column
  vector, so any caller indexing it by name must be checked. The reason is in
  the bound below -- once every tested column is capped, the vector is the
  bound everywhere and carries nothing the scalar does not.
- **A between-patient contrast is capped at its own compartment's patient
  count minus two**, at the fit site and again at the polish site, under BOTH
  `df.method` values. This covers every tested column carrying the condition
  factor, the three-way `CellType:Response:Niche` terms included. Without it
  the Satterthwaite df was anticorrelated with the evidence (Spearman -0.978
  against cell count) and anti-conservative where the data is thinnest.
- The `.satterthwaiteDF()` arm is kept and documented, bounded the same way,
  so an opt-in user is not misled either.

## The hard part of Phase 1

Every per-gene branch is a partition of the gene index set; batching turns
control flow into set operations over index vectors.

- **Active set**: an integer index, compacted when the active fraction falls
  below ~0.75, never per iteration. Compaction is performance, not semantics —
  tested by running with the threshold at 0 and at 1 and asserting identity.
- **Line search**: a per-gene `step` vector and a trial-round loop. One batched
  `dnbinom` column-sum per round evaluates every pending gene, so a block costs
  `max(halvings)` rounds rather than `sum(halvings)`. The accept test, the 1e-9
  slack, the 1e-6 floor and the 20-halving cap are copied verbatim.
- **Staleness stays per gene.** A batch-synchronous policy has the same fixed
  point but a different path, and at a 1e-8 relative stopping rule that turns
  the comparison against `engine = "gene"` from an equality test into a
  "both converged somewhere near" test.
- **A singular gene must not fail its batch.** Per-slice factorisation with a
  per-slice `ok`, so it drops to `fallback()` as today. This is the failure mode
  `.waldCauchyBlock()` has (`inference.R:405-414`), and it must not be repeated.

## Tests (written before the code)

1. `B = 1` reproduces `.polishGene()` to 1e-12 on every field and flag.
2. Batch-size invariance: 1, 2, 5, 13 on a 13-gene fixture, identical.
3. Compaction-threshold invariance: 0 and 1, identical.
4. **Mixed control flow in one batch** — a converged gene, a degenerate start,
   an all-zero gene, a gene needing four halvings, a gene whose ψ sits on the
   bound — reproduces the five run singly, flag for flag. Write this first.
5. A singular gene does not poison its batch.
6. Absorbed batched gram == dense batched gram on the tested columns.
7. GPU set (Phase 2) in `test-gpuPolish.R`, `skip_if_no_gpu()` plus a skip on
   MPS: float32 is refused, not warned about.

## Acceptance

Fresh recalibration was chosen, so parity with the current arm is not the gate.
The ladder is: unit suite with no tolerance loosened → `longtests` (the τ²
window 0.3–0.75 against a planted 0.49, the df anchors, the deflation test) →
the simulation benchmark as a *paired* arm through `polish_variants`, so the
engine is the only thing that differs → the cohort shuffle nulls scored against
the pass criteria in `package_fixed_design.R`'s header → then, and only then,
re-derive the |z| thresholds and re-render the cohort reports.

## Provenance

All work on `feature/batched-fitter`, in a worktree, off a clean `main`. Every
benchmark or cohort arm runs from a `freeze_snapshot.sh` snapshot, never the
live tree.
