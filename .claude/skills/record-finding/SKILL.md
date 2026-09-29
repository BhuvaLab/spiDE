---
name: record-finding
description: Record a measured result in every place this project keeps evidence - the study's findings log, NEWS.md, the calibration vignette, the relevant CLAUDE.md section, and a memory file - with consistent wording, readable names, the date, the job ids, the code snapshot or pinned library, and the seed. Use after a benchmark arm, a null-grid set, a cohort run or an ablation has produced a number that will be quoted or that changes a default, and before that number is quoted anywhere.
---

# Record a finding

The package's discipline is that every non-obvious default was chosen on
measurement and the measurement is written down. That decays in two ways:
a number lands in one place and not the others, or it is written without
the run that produced it. The `evidence-auditor` agent then finds the drift
months later. This skill writes the same finding to each location once, in
that location's voice.

## 1. Collect the provenance first

Refuse to write a number without all of these; ask for the missing one.

| item | where it comes from |
|---|---|
| date | today, absolute (YYYY-MM-DD), never "yesterday" |
| job ids | `sacct -u $USER -S <date> -X --format=JobID,JobName,State` |
| code | the pinned library's `research/libs/<name>.SNAPSHOT`, or `SNAPSHOT_PROVENANCE.txt` in `$SPIDE_PKG` (branch, HEAD, uncommitted diff); "live tree" is a defect to record, not a provenance; spiDE and SpaNorm versions (`fit@params$version` holds spiDE's) |
| seed | what the driver wrote alongside the output (block grid numbers, permutation seed, plasmode or simulator seed) |
| the arm's settings | engine (`slopes`/`sandwich`), `sigma`, `depth`, `covariates`, `strata`, `procedure`, `fdr`, `min.cells`/`min.patients`/`min.detect`, the index types and gene filter, number of patients, grids / permutations / replicates |
| cohorts | every cohort it was scored on; the three public cohorts carry decisions, YTMA enters at weight 0.25 |
| the canonical table it lives in | `research/bench2/tables/<scenario>.rds` with the method/engine/depth row it added, or the study's `out/*.rds` it was read from (e.g. `research/simplify/out/report_data.rds`) |
| what it is compared with | the arm, the prior measurement or the pre-registered criterion it is being read against |

A count of discoveries additionally needs the `calibration-check` skill's
output on that configuration's own nulls. Do not record a count that has not
been through it.

## 2. Decide the claim

Write one sentence: what was measured, on what, and what it decides. Be exact
about what it does **not** show. The record is full of numbers that failed in
one of these ways:
- **Circular:** the mixed model's 6x gene filter, scored on its own grid.
- **Measured on a predecessor estimator:** every mixed-model figure, now that
  spiDE 0.99.30 has replaced the model.
- **Read from a summary that could not tell two hypotheses apart:** the
  overall median null RMS passing while an expression band failed.
- **Chosen on YTMA alone:** the normalisation comparison.

If the finding has a caveat of that kind, the caveat goes in every location,
not only the long one.

Where the finding decides a pre-registered gate or rule (the study's
`README.md`), quote **that rule's** verdict unchanged. Put any re-weighted
reading (e.g. with YTMA down-weighted) next to it, never in its place.

Name engines, arms, scenarios and gates in words: the mixed model, the
sandwich engine, the slopes engine, two-stage OLS; the mechanism,
power-ceiling, calibration and power gates. No letter-and-number codes.

## 3. Write to each location

Read the existing text around the insertion point in each file first; match
its voice and do not repeat what it already says.

**The study's findings log** (`research/simplify/FINDINGS.md`,
`research/bench2/README.md`'s dated log, or the study's own file under
`research/<study>/`): the long form. It has a dated `##`/`###` heading, the
question, the setup with the full provenance table, the numbers as a table,
the interpretation, and a "what this does not show" line. This is the only
place the raw numbers need to appear in full. It is committed in the research
submodule, so the package's submodule pointer must move afterwards. A finding
about the archived mixed model goes to its own records (`research/mixed`,
`research/fdr-ordering/FINDINGS.md`), dated and scoped to spiDE <= 0.99.22.

**`vignettes/spiDE-calibration.Rmd`**: only if the finding changes a number
it quotes or adds one the package should quote. It is the **only** place in
the package that quotes benchmark numbers. Never paste the number into
another vignette, a roxygen block or `README.md`.

**`NEWS.md`**: only if the finding changes a user-visible default, argument
or behaviour. Add it under the current version's heading in the existing
style (bold lead, then the change), with a pointer to the calibration
vignette or the study rather than the raw numbers.

**`CLAUDE.md`**: only if the finding changes how future work should be done
(a default, a warning, a measured cost, a refuted hypothesis). Put it in the
section it belongs to, dated in parentheses as the others are. If it
supersedes an earlier statement, edit that statement rather than adding a
contradiction below it. Refuted hypotheses go into the list of things not to
re-run. CLAUDE.md may be open in another session: re-read the section
immediately before editing it.

**Memory** (`/home/uqdbhuva/.claude/projects/-scratch-project-mnt-S0249-R-projects-spiDE/memory/`):
- one file with the frontmatter the memory system expects (`name`,
  `description`, `metadata.type: project`);
- the fact, then **Why** and **How to apply** lines;
- links to related memories;
- a one-line pointer added to `MEMORY.md`.

Update an existing memory that covers the same question instead of adding a
second.

**The report**, if one reads the data:
- `research/simplify/report/simplify_top5.qmd` reads `out/report_data.rds`
  (assembled by `R/24_report_data.R`);
- bench2's reports read `research/bench2/tables/`.

Check that every one-row lookup still returns one row with the new arm
present (a second arm broke a report's lookup on 2026-09-08). Rebuild the
sites with the `build-site` skill.

## 4. Check before you stop

- The same number appears in every location it was written to, with the same
  units, the same cohorts and the same caveat.
- A canonical table was regenerated by its installer, never written directly;
  the `protect-canonical-tables` hook blocks the direct write for
  `research/bench2/tables/` and the archived
  `research/reports/benchmarks/tables/`.
- `vignettes/spiDE-calibration.Rmd` is still the only vignette quoting
  benchmark numbers.
- No study code appears in the new text.
- Run the `evidence-auditor` agent on the edited files when the finding
  changes or supersedes a quoted number.
