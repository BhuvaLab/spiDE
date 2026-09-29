---
name: calibration-check
description: Score whether spiDE (>= 0.99.30) test results are calibrated enough for their FDR to mean anything - the spread of null z per test, per expression band and per index type, the block and permutation null tails against the pre-registered limits, and condition dropout (informative missingness). Use after any real-cohort run or null grid set, before quoting any discovery count, or when deciding which index types to report.
---

# Calibration check

Run this **before quoting any discovery count**. `testSpiDE()` always
produces q-values; whether they mean anything depends on calibration that the
results table cannot show about itself. Calibration is judged on **null
grids** of the same configuration (cohort, bandwidth, engine, index types,
covariates); a real run on its own can only be checked for dropout and
described.

Decisions rest on the three public cohorts (GSE250346, GSE282639, GSE289194).
YTMA results (the ICI and the LUAD stage cohorts) enter at weight 0.25 (user
directive 2026-09-29: suspected data issues); a verdict that holds on YTMA
alone is not a verdict.

## Usage

```bash
Rscript .claude/skills/calibration-check/scripts/calibration_check.R \
  [--expr fit.rds] real.rds [...] [--block grid01.rds ...] [--perm perm001.rds ...]
```

- Each `.rds` holds a `SpiDEResults` (`testSpiDE()`/`spiDE()`), or the
  data.frame of `results(res, test = "both")`, or a list with either under
  `$results`/`$res`.
- Files before any flag are **real** runs, each scored on its own.
- Files after `--block` or `--perm` form one **null set**, scored together.
  `--null` is a generic null set scored without pass/fail limits.
- `--expr` names a `SpiDEFit` or `SpiDEResults` whose per-index mean
  expression defines the bands when the inputs are data.frames. Every block
  and perm grid of one fit has the same mean expression, so one object serves
  the whole set. Without it the first real `SpiDEResults` is used.
- The script reads slots with `attr()` and never loads spiDE. A stray
  `is()`/`inherits()`/`is.data.frame()` on the S4 object makes R load
  whichever spiDE the default library holds, which on this cluster is a
  0.99.19 mixed-model build.
- A spiDE <= 0.99.22 mixed-model object or table is refused with a pointer to
  `spiDEmixed::readSpiDE()` and the archived script at
  `research/mixed/claude/skills/calibration-check/`.

## Making the null grids

The two nulls of the simplification study (`research/simplify/R/common.R`):

- **block** (`.blockShuffle()`): each section's niche field is shifted
  toroidally by a random 20-80% of its extent (`set.seed(7000 + grid)`), each
  cell taking the niche row nearest its shifted position. Spatial smoothness
  survives, the niche-expression link does not. Every grid is a **refit**
  (`fitSpiDE()` on the shifted `reducedDim`) and scores both tests.
- **perm** (`.permLabels()`): the condition label is permuted across patients
  **within slide**, so a slide effect cannot pass for condition. For the
  slopes engine the fit never sees the condition, so a permutation reruns only
  `testSpiDE()`: add the permuted label as a column of `fit@patients` and test
  it (the package's own relabelling test in `tests/testthat/test-engine.R` does
  exactly this). The sandwich engine uses the condition at fit time, so each
  permutation needs a refit.

Save a perm grid as `results(res, test = "both")` (small), and pass one real
`SpiDEResults` via `--expr` for the bands. The study's validation sets were
ten block grids and 1,000 permutations per cohort. With fewer than ten grids
the script says so, because a tail fraction that coarse cannot resolve a
0.07 limit.

## What it reports

`z = sign(t) * qnorm(p / 2, lower.tail = FALSE)`, a standard normal under the
null whatever the test's df (the slopes engine's df differ per gene, capped
by the Kish effective patient count; the sandwich engine's are
Bell-McCaffrey).

**(a) Global spread per test**: `n`, `sd(z)`, `RMS(z)`, fraction `|z| > 1.96`.
On a null these should be about 1, 1 and 0.05. On real data signal inflates
them, so the script labels them descriptive.

**(b) Per-band RMS of z.** RMS per (gene, index) over niches (and over
grids), then the median within each quintile of mean expression within the
index type. This is the study's band criterion. The overall median hides
the failure it exists to catch: the mixed model's null inflation was
monotone in expression in every cohort.

**(c) Per index type**, the spread of (a) by index and test.

**(d) Null tails** (null sets only): the fraction of grids with at least one
BH call at 0.05 for the pooled test, the condition test over every triplet
and the condition test over the filtered family. The families are
recomputed from `p` with `.bhFamilies()`'s rule (pooled q < the fit's `fdr`
admits a triplet), so a set run with either `procedure` scores alike. With a
real run alongside, the script also prints real against mean null calls
(empirical FDP). That comparison is meaningful only when the grids are that
run's own nulls.

**(e) Condition dropout** (real runs only). For each index type, the number
of patients per condition carrying a usable slope, and a Fisher test of
inclusion against condition.
- **Slopes engine:** usable means the patient passed `min.cells` for the
  index type and has a finite slope for at least half of its (gene, niche)
  pairs. A patient's slope is dropped for a gene when its niche columns are
  collinear within its cells.
- **Sandwich engine:** it stores no per-patient slopes, so only the index-level
  inclusion (`min.cells`) is checked.

If inclusion depends on condition, the contrast is confounded at its root,
and no threshold or variance correction fixes informative missingness.

## The limits

These are pre-registered in `research/simplify/README.md` and quoted in
`vignettes/spiDE-calibration.Rmd`:

| criterion | limit | scored on |
|---|---|---|
| per-band median RMS of null z | in [0.90, 1.10] in every band | condition test on both nulls; pooled test on the block null |
| block tail | <= 0.10 | pooled and condition tests (and the filtered family) |
| perm tail | <= 0.07 | condition test over every triplet and over the filtered family |

The pooled test never sees the labels, so it is not scored on the
permutation null. The script leaves the pooled rows out of a perm set's
spread and flags them as a label leak if their calls differ across
permutations.

For reference, the shipped slopes engine on the validation nulls, in every
cohort (`research/simplify/FINDINGS.md`, 2026-09-29): bands 0.909-0.984,
worst perm tail 0.034, and worst block tail 0.10 (GSE250346, pooled; at the
limit).

## Interpreting the result

- **Quote a discovery count only from a configuration whose nulls pass.**
  Quote it with the null it passed and its grid count. Otherwise, quote real
  calls against the null's mean calls (empirical FDP), never the raw BH count.
- **The filtered family is the primary condition procedure.** Its filter is
  label-free, so it is exact under relabelling, but it is still scored on the
  perm null. In the study the sandwich engine's filtered procedure failed
  there once: 9 of 100 permutations on GSE289194, a 95% interval of roughly
  0.04-0.16.
- **A slide-confounded condition leaks through the perm null.** On the YTMA
  stage cohort, stage is confounded with TMA. The two-stage OLS arm's filtered
  procedure called on 0.31 of null grids until TMA entered the design, which
  brought it to 0.065. `strata` is that lever in the package: `fitSpiDE()`
  uses it for the sandwich engine and `testSpiDE()` for the slopes engine.
- **Dropout confounding outranks everything above.** Report an index type
  with condition-dependent inclusion as confounded, whatever its bands.

## If nothing passes

That is a real answer, not a failure of the script. It means this cohort and
configuration cannot support the claimed FDR. The levers, in order:
1. `strata` for a slide- or batch-confounded condition.
2. The pre-registered, label-free restrictions of the study: index types
   covering at least 90% of patients with at least 10 cells, and gene
   detection of at least 0.10 within the type (`min.detect`).
3. The other engine. The slopes engine is stronger on effects that vary
   between patients; the sandwich engine is stronger on shared effects in
   sparse data.

Any restriction chosen after seeing the nulls must be stated as post-hoc.
