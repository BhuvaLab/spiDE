---
name: calibration-check
description: Score whether a spiDE results table is calibrated enough for its FDR to mean anything, and identify which index cell types are usable. Reports sd(t), per-index breakdown, dropout-vs-condition confounding, and re-FDRs the valid subset. Use after any real-cohort or permutation run, before quoting discovery counts, or when deciding which index types to report.
---

# Calibration check

Run this **before quoting any discovery count**. A spiDE results table always
produces q-values; whether they mean anything depends on calibration that the
table itself does not report.

## Usage

```bash
Rscript .claude/skills/calibration-check/scripts/calibration_check.R \
  <results.rds> [fit.rds]
```

The fit object is optional but strongly preferred — the dropout-vs-condition
test and the valid-subset selection both need `@diagnostics$inclusion`.

## What it reports, and how to read it

**1. Global calibration.** `sd(t)` should be ~1.0 and `frac|t|>1.96` ~0.05
under the null. Above ~1.2, BH is not valid on that table and the discovery
count is an artefact. Measured on YTMA at sigma=30: sd(t) = 1.363,
frac = 0.145, which corresponded to ~5,400 "significant" genes of a 13,348 gene
panel — roughly 40% of the panel, which is the tell.

**2. Per-index calibration.** Inflation is driven by **cells per (patient,
index) subset** — measured correlation −0.838, versus −0.659 against patient
count. Ruled out as drivers: bandwidth (sd 1.363 at sigma=30 vs 1.364 at
sigma=70) and patient count (restricting to >=30 patients moved 1.363 only to
1.277). Mechanism: ~12 niche columns fit to as few as 30-40 cells gives a
near-singular design whose Fisher variances understate uncertainty.

**3. Dropout vs condition.** A separate and more serious problem. If which
patients carry enough cells of an index type depends on their condition, the
contrast is confounded at its root and **no threshold or variance correction
fixes it**. Measured on YTMA: B cell p=0.027, DC p=0.006, Monocyte p=0.016.

**4. The valid subset**, re-FDR'd on its own. On YTMA this was Tumor +
Fibroblast (55/55 patients, sd(t) 1.05 and 1.12), giving 4-6 genes at q<=0.05
against ~5,400 in the unrestricted table.

## Interpreting the result

- **Do not** report the unrestricted discovery count when global sd(t) > 1.2.
- **Do** state the restriction as pre-specifiable criteria (complete patient
  inclusion; calibrated variance) — both are independent of the outcome, which
  is what makes a post-hoc filter defensible.
- **Prefer hits replicating across bandwidths.** Bandwidths are not independent
  tests, but a hit surviving several is less likely to be a thin-subset
  artefact. On YTMA, 5 of 30 genes at q<=0.10 replicated in >=2 bandwidths.
- A raw `p<0.05` fraction near 0.05 in the valid subset is the sanity check
  that the restriction worked.

## If nothing passes

That is a real answer, not a failure of the script. It means the cohort cannot
support this estimator at the chosen `min.cells` and niche resolution. The
levers are: raise `min.cells` (costs index types — on YTMA, reaching 100 drops
the entire T cell / DC / Monocyte / Mast compartment), or reduce niche columns
for thin index types so the p/n ratio recovers without deleting them.
