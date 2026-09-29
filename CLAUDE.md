# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this package does

spiDE finds **neighbourhood-dependent differential expression** in spatial transcriptomics. Within an
*index* cell type it asks, for every (gene, index type, niche type) *triplet*, whether the gene's
expression changes with the local density (the *niche*) of the niche type **within patients**
(the **pooled** test), and whether that dependence differs between two patient-level conditions
(the **condition-specific** test). Patients are the unit of replication: a triplet's evidence is
how its niche slope varies between patients, never the number of cells. Built on the author's
**SpaNorm** package, whose `polishNB()` is the per-gene negative binomial solver.

From 0.99.30 the package is a rebuild around two per-patient engines. The joint mixed-effects
model of spiDE <= 0.99.22 is archived (see "The archive"); it failed calibration on every cohort it
was measured on.

## Names (user rule, 2026-09-29)

Every engine, arm, scenario and gate carries a human-readable name in code, docs, reports and
commit messages, never a code. The simplification study's scripts keep their historical file names
(`m1.R`, `h1.R`, ...); this is the mapping:

| study label | name |
|---|---|
| M0 | the **mixed model**: spiDE <= 0.99.22, archived as `spiDEmixed` |
| M1 | the **sandwich engine** (`engine = "sandwich"`) |
| M3 (candidate A) | the **slopes engine** (`engine = "slopes"`, default): effective-patient df + trended prior |
| candidate B | trended prior with the effective-patient df on the pooled test only (rejected) |
| M2 | **two-stage OLS** (historical; simplification study only) |
| N1 / N3 / N3c / N4 | depth `"loglib"` / `"spanorm_offset"` / `"spanorm_covariate"` / SpaNorm PAC (reference) |
| LS spline / nonlinear LS | depth `"spatial_spline"` / `"nonlinear"` (being benchmarked; not yet in the package) |
| D0-D5 | depth scenarios **iid**, **smooth technical**, **niche-aligned technical**, **niche-driven content**, **aligned + driven**, **strong gene response** |
| G0-G3 | the **mechanism**, **power-ceiling**, **calibration** and **power** gates |

Function names must be meaningful in the spiDE context; clashes with the archived API are fine.

## The archive (read before touching old code or old objects)

- **`spiDEmixed`** (`research/mixed/`, research submodule) is spiDE 0.99.22 verbatim, renamed:
  every function, test, long test, fixture, vignette and design document, plus the old CLAUDE.md
  (`CLAUDE-spiDE-0.99.22.md`), NEWS and `.claude/` (`claude/`). Installed with spiDE 0.99.22 and
  SpaNorm 1.7.14 in the pinned library `research/libs/mixed-final` (provenance in
  `mixed-final.SNAPSHOT`); the frozen source tree is `spiDE_snapshots/legacy-mixed-0.99.22`
  (tag `v0.99.22-mixed-final` in this repo, `pre-spiDEmixed` in research).
- **Any code that calls the old functions** (`fitSpiDE(random = ...)`, `polishSpiDE()`,
  `nicheDesign()`, `computeSizeFactors()`, `compositionTest()`, `fits()`, `bandwidths()`, the old
  `testSpiDE()`/`spiDE()`) calls them as `spiDEmixed::` or loads the pinned library first:
  `.libPaths(c("/scratch/project_mnt/S0249/R_projects/spiDE/research/libs/mixed-final", .libPaths()))`.
  Every research script was pinned this way on 2026-09-29 (`research/mixed/move/research_repoint.tsv`).
- **Saved objects** from <= 0.99.22 are refused by this version (`.isLegacy()` / `.assertCurrent()`
  in `R/AllClasses.R`) with a pointer to `spiDEmixed::readSpiDE(path)`, which re-homes them.
- **`research/mixed/move/check_move.R`** (run from this repo's root) proves nothing was left behind:
  every file at the tag, every symbol body, every `test_that()` and Rd name, a denylist grep of this
  tree, and the research pins. Re-run it after any change that touches either side.
- `spiDEtwostage` (`research/twostage/`) is the two-stage estimator archived at 0.99.19; its test
  calls `spiDEmixed::fitSpiDE()`.

## Commands

Standard R/Bioconductor package (roxygen2, testthat edition 3). SpaNorm >= 1.7.14 is required; on
Bunya it is in `research/libs/simplify` (put it first on `.libPaths()`).

```r
devtools::load_all(); devtools::document(); devtools::test()
testthat::test_file("tests/testthat/test-engine.R")
devtools::check(); BiocCheck::BiocCheck()
```

Build the site with `Rscript -e 'pkgdown::build_site()'`. It carries only the three vignettes; the
validation reports are published by the research repo (`research/reports/benchmarks/build_site.R`
into `research/docs/`), which also holds the archive of the mixed model's reports. **The package
site never links an archived report** (user rule); a current report goes into the Validation menu
of `_pkgdown.yml` as an absolute research-site URL.

Regenerate `data/toySpiDE.rda` with `source("data-raw/make_toySpiDE.R")` (16 patients x 150 cells,
seed 7). `longtests/` holds slow checks, not run by `devtools::test()` or CI.

Project automations live in `.claude/` (allowlisted in `.gitignore`): hooks `r-parse-check.sh`
(parses every edited `.R` file), `protect-running-drivers.sh` (blocks in-place writes to scripts
while SLURM jobs run; write `<file>.new` and `mv` it), `protect-canonical-tables.sh` (only
installers write canonical benchmark tables), `roxygen-drift.sh`; agents
`design-invariant-reviewer`, `numerical-robustness-reviewer`, `evidence-auditor`,
`slurm-run-triage`; skills `calibration-check`, `run-benchmark-arm`, `hpc-job-sizing`,
`record-finding`, `build-site`, `deploy-bioc`. CI: `.github/workflows/check-bioc.yml` (R CMD check +
BiocCheck, four configurations) and `pkgdown.yaml`. Dot-separated argument names (`min.cells`,
`min.detect`) follow SpaNorm's.

## Architecture

**Pipeline:** `buildNiches()` -> `fitSpiDE()` -> `testSpiDE()`; `spiDE()` chains them. One
bandwidth per fit. `testNicheAbundance()` and `spiGSEA()` sit beside it.

1. **Niches** (`R/buildNiches.R`): per sample and bandwidth, a Gaussian KDE of every cell type at
   every cell (`spatstat.explore::densityfun`), stored as `reducedDim(spe, "Niche<sigma>")`;
   `mergeNiches()` pools columns and records the groups in `metadata(spe)$spiDE_niche_groups`
   keyed by reducedDim name, which `.isSelfNiche()` honours (an index type is never tested against
   its own niche, merged or not).
2. **Fit** (`R/fitSpiDE.R`, `R/fit-design.R`, `R/fit-glm.R`): per index type, the type's cells
   in patients with >= `min.cells`, genes detected in >= `min.detect` of them. Design
   `[patient intercepts | centred covariates (loglib) | log1p niche densities | (condition x niche)]`
   (`.indexDesign()`); one NB GLM per gene with `SpaNorm::polishNB()`, the patient block absorbed
   by Schur complement, profile-ML dispersion. By Frisch-Waugh-Lovell every niche coefficient is a
   within-patient slope; the patient intercepts absorb the condition main effect and every
   between-patient composition effect.
   - **slopes engine** (default): one condition-free fit, then `.patientSlopes()`: each patient's
     one-step NB slope from the shared fit, with `v_model` and a within-patient spatial tile
     sandwich `v_tile` (tiles of `3 * sigma`); `.patientFactor()` = per-(patient, niche)
     `max(1, median_g v_tile / v_model)` for the pooled test.
   - **sandwich engine**: the pooled model and, with a condition, the condition x niche model;
     `.sandwichCR2()` = low-rank CR2 + Bell-McCaffrey df (`strata` enters as strata x niche
     nuisance columns).
3. **Test** (`R/testSpiDE.R`, `R/test-slopes.R`): slopes engine: per (index, niche) column, limma
   across patients weighted by `1 / (v + tau2_DL)`, `eBayes(trend = TRUE, robust = TRUE)` on log
   mean expression, df = `min(df.total, Kish n_eff - p)` (`.slopeColumnTest()`). Pooled test uses
   `v_model * factor`; condition test uses `v_tile`. `.bhFamilies()`: pooled BH over every
   triplet; condition BH over the triplets whose pooled q < `fdr` (`procedure = "filtered"`,
   default) or all. `results(test = )` reads one table; `patientSlopes()` the per-patient slopes.

**Classes** (`R/AllClasses.R`): `SpiDEFit` (engine, sigma, condition, patients, index, params) and
`SpiDEResults` (table, condition, contrast, procedure, fdr, fit). Generics in `R/AllGenerics.R`,
input checks in `R/checkers.R` (extend them; do not duplicate validation).

### Invariants (the `design-invariant-reviewer` agent enforces them)

1. Per-gene fit work is exact and may be blocked over genes (`.bpGenes()`, `polishNB(BPPARAM)`).
2. Anything shared across genes (the per-patient factor, the eBayes prior and trend, the BH
   families, the filter set) is computed on the whole family, never inside a gene block.
3. The slopes engine's fit never sees the condition, so a permutation reruns only `testSpiDE()`;
   the pooled test and the filter are invariant to relabelling (tested).
4. Absorbed blocks nest within patients.
5. Counts are never densified whole (per-gene rows, or genes x patients).
6. One BLAS thread per forked worker (`.bpGenes()` via RhpcBLASctl).
7. A bad gene or patient drops out (`NA`); it never aborts a run. Every per-gene inversion is
   guarded.

## Evidence and decisions

Most defaults were chosen on measurement. **The calibration vignette
(`vignettes/spiDE-calibration.Rmd`) is the only place in the package that quotes benchmark
numbers.** Before changing a default, read:

- `research/simplify/` — the simplification study (2026-09-28/29) that chose the engines:
  `README.md` (pre-registered gates and rules), `FINDINGS.md`, the report
  `report/simplify_top5.qmd` (data `out/report_data.rds`, built by `R/24_report_data.R`),
  published on the research site. `tests/package_vs_prototype.R` shows the packaged slopes engine
  reproduces the prototype (GSE282639: identical calls, max |dt| 6e-4).
- **Calibration** is judged on two real-data nulls per cohort: the **block null** (each section's
  niche field shifted toroidally, `.blockShuffle()`) and the **permutation null** (condition
  permuted across patients within slide, `.permLabels()`), against the calibration gate: per
  expression band median null RMS of z in [0.90, 1.10], block tail <= 0.10, permutation tail
  <= 0.07. The slopes engine passes in all five comparisons; the sandwich engine all but one
  borderline filtered-procedure tail (GSE289194).
- **YTMA is down-weighted** (user directive 2026-09-29: suspected data issues). The three public
  cohorts (GSE250346, GSE282639, GSE289194) carry decisions; YTMA enters sensitivity analyses at
  weight 0.25. Under the pre-registered equal-weight rule the sandwich engine would have stayed the
  recommendation; the slopes engine is the default under the directive. Decisions are scored on
  several cohorts, never YTMA alone.
- **Depth.** `depth = "loglib"` (one slope per gene on centred log library size) is the only
  option shipped. Normalisation is judged **against the estimand, never by residual library-size
  association**: library size confounds biology in spatial data, which is why SpaNorm balances LS
  removal against biology with competing splines. The simplification study's YTMA-only verdict
  (SpaNorm offset equivalent to raw + loglib) is withdrawn. `research/bench2/` holds the new
  benchmark (spline-based, gene-specific LS simulator fitted to real templates; depth options
  `spatial_spline`, `nonlinear`, `spanorm_offset`, `spanorm_covariate`); its README is
  pre-registered before any run. SpaNorm's block design for the spatial LS spline is on SpaNorm
  branch `feature/block-design`.
- **Composition.** Between-patient association of niche abundance with expression is real signal
  that the patient intercepts absorb by design (it inflated the mixed model's slopes before 0.99.17;
  archived CLAUDE.md, "The cause"). `testNicheAbundance()` tests it on pseudobulk; never report it
  as niche-dependent DE, or the reverse. `.toySPE(composition = k)` needs `n_per >= 200` and
  bandwidth 50 for the covariate to carry it.

### Open items

- `spiGSEA()` is experimental: no gene-set null has been run.
- The niche covariate's own spatial autocorrelation is unmodelled; the block null keeps it.

## Data hazards

- **Counts must be raw integers** (`checkCounts(integer.only = TRUE)`). The YTMA v10 objects set
  `counts(spe) <- 2^logcounts(spe) - 1` (`batch_nichede_v10.R:224`), infeasible values that also
  destroy sparsity; `research/fdr-ordering/R/make_raw_spe.R` rebuilds a raw twin.
- One sample per patient in `sample_id`; several sections per patient go in `section` (tiles are
  laid per section).
- The condition must be constant within a patient (`.conditionCoding()`); where it is confounded
  with slide, pass `strata`.

## HPC practice

**A warning about the benchmark harness.** A sweep whose tasks load `SPIDE_HOME` from the **live
working tree** is not one experiment: tasks start at different times, so edits mid-run give different
tasks different code. A previously reported result ("two-stage null inflation grows with S, 0.078 →
0.127") came from such a run and **did not survive re-measurement on a frozen snapshot** — the paired
ablation shows no trend with S in either configuration (p = 0.82 / 0.62). Always pin `SPIDE_PKG` to a
snapshot (`.claude/skills/run-benchmark-arm/scripts/freeze_snapshot.sh`). The same applies to **every
`Rscript` driver**, not only the package: R parses a script file incrementally through an 8 KB stdio
buffer, so a task that is inside its multi-hour fit reads the *next* expression from whatever the
live file holds by then, at the old byte offset. Two losses on 2026-09-07: 255 benchmark tasks wrote
their result and then died on parse garbage (the outputs were intact; the array was marked FAILED and
its `afterok` aggregation left at `DependencyNeverSatisfied`), and at 22:26 an in-place edit of
`research/fdr-ordering/R/package_fixed_design.R` killed the eight cohort tasks whose fit ended after
it, six to seven hours each. Every sbatch therefore now copies its driver to a task-private temp
file at task start and runs the copy, and `freeze_snapshot.sh` puts frozen copies of both drivers
under `$SPIDE_PKG/drivers/`, which the sbatch scripts prefer. Two facts that made the rescue
possible: a `git checkout` or `sed -i` writes a **new inode** and cannot touch a running reader,
while the Edit tool and a shell `>` rewrite **in place** and can; and a reader that has not yet
refilled its buffer (position 8192 in `/proc/<pid>/fdinfo`, visible through `srun --overlap
--cpu-bind=none` on the task's raw job id) is saved by truncating and rewriting the file in place
with the bytes it expects. Edit the harness between arrays, never during one.

Also: set `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` in every sbatch that forks workers (a forked
worker inherits the parent's BLAS threads, and oversubscribed workers ran the old polish 9x
slow); inside a Claude session pass `sbatch --export=ALL,...` (the session sets
`SBATCH_EXPORT=NONE`) and read the first task log before calling a submission done; pilot one task
before an array and size from `sstat`/`seff` (`hpc-job-sizing`). Commit, never push, unless asked.

## Where the evidence lives

- `vignettes/spiDE-model.Rmd` (the model), `vignettes/spiDE-calibration.Rmd` (the numbers).
- `research/simplify/` (the engines), `research/bench2/` (the depth benchmark, in progress),
  `research/public/` (the three public cohorts: builders, nulls, `FINDINGS.md`).
- `research/mixed/` (the archived mixed model and its whole evidence trail, including
  `CLAUDE-spiDE-0.99.22.md`), `research/reports/benchmarks/` and `research/docs/` (its six reports,
  now the research site's archive), `research/fdr-ordering/`, `research/fdr-triplet/`,
  `research/plasmode/` (the mixed-model era investigations; pinned to the legacy library).
- `.claude/skills/calibration-check/` — run it before quoting any discovery count.
