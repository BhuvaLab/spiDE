# Research-site reports reorganisation — plan for review

Date: 2026-09-07. Status: **awaiting approval**; nothing has been changed.
Scope: the reports rendered to https://bhuvalab.github.io/spiDE-research/
(`research/reports/benchmarks/*.Rmd`, `build_site.R`, `docs/`). The
repository's directory layout is not in scope.

## 1. What is on the site now, and what is wrong with it

| report | lines | state |
|---|---|---|
| Simulation study | 1,265 | coherent in its first half (the simulator) and second half (calibration, power, FDR, robustness, scalability), but its design comparison is a trailing section that reads as the history of a promotion ("which model design does the package ship", "what the arm actually varies", the corrected reading), and the 0.99.17 switch arms were bolted on after it. Calibration evidence is split between a null-calibration section, that trailing section, and the switch-arm subsection. |
| The niche-shuffle null | 199 | **superseded.** Written when the verdict was "no detectable niche-dependent DE for either estimator"; the cause of that verdict has since been found and fixed, and the verdict reversed. Nothing on the site says so. The reader-facing account of the whole investigation exists only as `fdr-ordering/REPORT.md` (discovery-ordered) and an artifact page that is in `docs/` but not in the index. |
| Two-stage benchmark | 640 | sound but carries two generations of measurement side by side ("the `nb` reference arm on the published simulator", "the structured library-size re-run", withdrawn claims); the estimator's real-cohort calibration (cells per subset, inclusion confounding) sits in the shuffle-null report instead of here. |
| Mixed-fit speedups (`re.prop`) | 200 | documents a **retired** trade-off: subsampling biases `tau2` downward and the default is off. Its final section duplicates the simulation report's df calibration. |
| Combining niche p-values | 361 | fine; light trim. |
| Gene-set inference | 368 | fine; light trim. |
| Between-sample stratum | 321 | a negative result that is not on the index at all. |

The same faults as the vignettes had: discovery order, redundancy, and a
missing synthesis. And one fault the vignettes did not have: a report on the
site states a conclusion that is now known to be wrong.

## 2. Proposed set: six reports, in reading order

### 1. Simulation study — calibration, power, FDR (restructured, ~1,000 lines)

Keep Part 1 (the simulator) as is. Restructure Part 2 so that **calibration is
one section with one subsection per default**, in the order the model vignette
introduces them: the random modes (`none` / `intercept` / `slope`); the
reference df (`between` / `satterthwaite`); the design term
(`CellType:condition`, the two-design comparison and the `ctresp` family,
with the corrected reading and *no* promotion history); the 0.99.17 switches
(nested intercept, convergence; the caveat that the simulator plants no
composition effect). Then power and FDR, per-category recall, spatial
robustness, scalability, summary. Remove: "which model design does the
package ship", "what the arm actually varies", and every sentence dated by a
promotion or a re-run.

### 2. The real cohort — the null tail, its cause, and the fix (NEW, ~450 lines; replaces *The niche-shuffle null*)

Written for a reader, by argument. The cohort and the shuffle null (`free`
and `block`, and why the patient-label permutation cannot detect a spurious
slope); the symptom — the real data indistinguishable from their null, the
heavy tail per-gene and confined to bright genes in populous types; what it
was not, as one table (the refutation ledger, one line per candidate); the
cause — the composition association the shuffle preserves, demonstrated with
the fit held fixed, and the Frisch–Waugh statement; the fix and its
validation on all grids; the reversed verdict and its stress tests (three
artefact arms, the count is not a set of findings, the one robust triplet);
the composition test on the same cohort (the immunoglobulin-in-tumour signal
is composition); the substrate (raw counts and the library-size slope); what
remains open. Reads its figures from `fdr-ordering/figures/` (15 tracked
PNGs) and its numbers from `fdr-ordering/summary/` (11 CSVs) plus the
canonical tables; every quoted number lands as an exactly-one-row lookup, the
site's existing convention. `fdr-ordering/REPORT.md` becomes a pointer;
`FINDINGS.md` stays the notebook. The artifact page comes out of `docs/`.

### 3. The two-stage estimator (condensed, ~450 lines)

One grid: the structured-library-size run with all three stage-1 arms, and
the `ols` recommendation stated once. Keep methods, null, raw power,
delivered discoveries, cells per sample, layouts and the negative control,
the niche-independent response, the adversarial scenarios. **Add** the
real-cohort section moved from the shuffle-null report: calibration governed
by cells per subset (the $r = -0.84$ relation, Tumor vs Mast), the cost of a
flat `min.cells`, dropout confounded with condition. Remove the
published-simulator `nb` reference section and the withdrawn claims; keep a
two-line methods note that the earlier grid existed and why it was replaced.

### 4. Combining niche p-values — Brown vs Cauchy (trimmed, ~300 lines)

As is, minus the "in the full pipeline at realistic scale" section's
repetition of the simulation report's numbers (one paragraph + link).

### 5. Gene-set inference — calibrating spiGSEA (trimmed, ~320 lines)

As is; the cost section becomes a paragraph.

### 6. What was tried and rejected (NEW, ~350 lines)

One section per refuted alternative, each with the question, the measurement
that decided it, and the number — the page every "do not rebuild" note can
point to. Absorbs two existing reports: *Between-sample stratum* (its
factorial and decomposition figures, from `reports/data/`) and *Mixed-fit
speedups* (the `re.prop` trade-off and the downward `tau2` bias). Adds, from
the notebooks: the nested fixed design; the QL dispersion; the gene filter's
circular 6×; cell-type-specific size factors; the random slope as a fix for
the composition confound; the niche-only design's recall (cross-linking
report 1); `winsor` and the log1p niche transform (from `niche-transform/`).

## 3. Cross-cutting rules

- Reading order on the index is 1–6, each with a one-line description that
  says what question it answers, and the two package vignettes linked above
  them.
- A conclusion appears in the report that owns the question and is linked
  from elsewhere; the switch-arm numbers live in report 1, the cohort numbers
  in report 2, the two-stage numbers in report 3.
- No promotion or re-run history in reports 1–5; anything superseded that is
  worth remembering goes to report 6.
- Numbers are read from tables (exactly-one-row lookups), never typed, as
  now.

## 4. Execution and verification

1. Report 2 first (the largest gap), from `REPORT.md` and the artifact page;
   render.
2. Report 1 restructure; report 3 condense with the moved cohort section;
   render both.
3. Report 6 from the stratum and mixed-benchmark reports plus the notebooks;
   then retire those two Rmd files (git history keeps them); render.
4. Trim reports 4 and 5; update `build_site.R` (order, titles, descriptions),
   rebuild `docs/`, remove the artifact page from `docs/`; update
   `_pkgdown.yml`'s validation menu and CLAUDE.md's "Where the evidence
   lives"; commit the submodule and the package pointer; push.

Estimated effort: report 2 is a day of writing; the rest is a second day.
Rendering the site takes minutes on the shipped tables.

## 5. Decisions for the author

1. Six reports as above (recommended) versus keeping the shuffle-null report
   alongside the new cohort report with a "superseded" banner.
2. Absorb the `re.prop` and stratum reports into *What was tried and
   rejected* (recommended) versus leaving them as standalone appendices.
3. Whether report 2 should also carry the plan of the remaining open work
   (the two candidate cures, the transcriptome run) or leave that to the
   package NEWS. Recommendation: a short closing section, since it is where
   the evidence for those items lives.
