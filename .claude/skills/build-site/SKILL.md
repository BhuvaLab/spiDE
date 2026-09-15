---
name: build-site
description: Build the spiDE package site (pkgdown) and the research site (research/reports/benchmarks/build_site.R) locally in the required order - sync the articles from the research submodule first, then build - and verify every report rendered. Use before pushing main, after editing a benchmark report or a canonical table, after `_pkgdown.yml` changes, or when either site build fails with "vignette missing from index" or "expected 1 row".
disable-model-invocation: true
---

# Build both sites

Two sites render the same benchmark reports and have each failed at least once
after CI was already green. Run this before pushing `main`; both builds are
local and touch nothing tracked.

## Run it

```bash
bash .claude/skills/build-site/scripts/build_site.sh            # both sites
bash .claude/skills/build-site/scripts/build_site.sh package    # pkgdown only
bash .claude/skills/build-site/scripts/build_site.sh research   # research only
```

The script (1) runs `vignettes/articles/pkgdown-sync-articles.R`, which reads
the article list from `_pkgdown.yml`, copies exactly those reports and deletes
any stale copy; (2) builds pkgdown into a scratch destination (`SITE_DEST`,
default under `$TMPDIR`), so `docs/` and the tracked tree are untouched;
(3) runs `build_site.R` from the research repo root into `research/docs/`;
(4) checks that every report's `.html` mtime is newer than the build start,
because a render that aborts part-way leaves the earlier pages fresh and the
later ones stale with no error in the summary line.

Never run `pkgdown::build_site()` without the sync: that is the first
failure below.

## The known failures and their remedies

| Error text | Cause | Remedy |
|---|---|---|
| `In _pkgdown.yml, 1 vignette missing from index: "articles/spiDE-…"` | A stale gitignored copy under `vignettes/articles/` from a sync made before that report was retired (the two-stage benchmark, 0.99.19). `git pull` cannot remove it and pkgdown cannot ignore it. Per clone, not repo state. | Run the sync (the script does); it deletes copies the index does not list. |
| `indexed articles without a report source: …` (from the sync) | `_pkgdown.yml` lists a report the submodule no longer has. | Fix the yaml, or `git submodule update` if the submodule is behind. |
| `research submodule not checked out` | `research/` is empty on this clone. | `git submodule update --init`. |
| `gv(): expected 1 row, got 2` or any one-row lookup failing during a knit | A canonical table gained a second `design` arm (every `gsea_*.rds` since 2026-09-08 carries `celltype-response` and `nested-polished-ql`) and a report's lookup does not filter on it. | Add a design filter in that report's `bench()` (`.design_current`, falling back to the first design present), as the GSEA and simulation reports do. Any new arm must be checked against every report reading that table. |
| Report `.html` older than the build start with no error | The render stopped part-way; the summary line only reports the last report. | Read `build_site.R`'s output for the report after the last fresh one. |
| `pkgdown::build_site()` fails on `install = TRUE` | The package does not install (usually stale `NAMESPACE`/`man/`). | `Rscript -e 'devtools::document()'`, then rerun. |

## After it passes

Push. Poll the GitHub Actions API no faster than every three minutes: the
cluster shell has no `GITHUB_PAT`, the unauthenticated limit is 60 requests
per hour, and a 30 s poll locks you out for an hour.

Do not commit `research/docs/` from the package repo; it belongs to the
research submodule and is committed there, with the submodule pointer
updated in the package afterwards.
