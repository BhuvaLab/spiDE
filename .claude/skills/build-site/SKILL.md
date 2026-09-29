---
name: build-site
description: Build the spiDE package site (pkgdown, the three vignettes) and the research site (research/reports/benchmarks/build_site.R, current reports plus the mixed-model archive) locally, and verify every page was written this run. Use before pushing main, after editing a vignette, `_pkgdown.yml` or a research report, after a study's report data change, or when either site build fails.
disable-model-invocation: true
---

# Build both sites

The package site (`bhuvalab.github.io/spiDE`) carries only the package's
three vignettes: `spiDE`, `spiDE-model` and `spiDE-calibration`. The
validation reports live on the research site
(`bhuvalab.github.io/spiDE-research`). That site is built from the `research/`
submodule by `research/reports/benchmarks/build_site.R` into `research/docs/`,
in two sections:
- **current** reports for spiDE >= 0.99.30, e.g. `simplify-top5.html`;
- the **"Archive: spiDE <= 0.99.22 (mixed model)"**.

The package site links only to current reports, from its `_pkgdown.yml`
Validation menu, never to the archive. The article sync of spiDE <= 0.99.22
(`vignettes/articles/`, `pkgdown-sync-articles.R`) is gone. Both sites have
failed at least once after CI was already green, so build them locally
before pushing `main`.

## Run it

```bash
bash .claude/skills/build-site/scripts/build_site.sh            # both sites
bash .claude/skills/build-site/scripts/build_site.sh package    # pkgdown only
bash .claude/skills/build-site/scripts/build_site.sh research   # research only
bash .claude/skills/build-site/scripts/build_site.sh research --render-archive
```

The script works in four steps:
1. **pkgdown.** It runs `pkgdown::build_site()` (the equivalent of
   `Rscript -e 'pkgdown::build_site()'`, with `install = TRUE`) into a scratch
   destination (`SITE_DEST`, default under `$TMPDIR`), so `docs/` and the
   tracked tree are untouched. It prepends `research/libs/simplify` to
   `R_LIBS` (override with `SPIDE_LIB`), because the default user library
   holds an older SpaNorm than the package's `SpaNorm (>= 1.7.14)`.
2. **Check the current reports' renders.** `build_site.R` does not render
   them: it copies each from where it is rendered (e.g.
   `research/simplify/report/simplify_top5.html`). The script fails if a
   render is missing, and warns if it is older than its `.qmd`.
3. **The research site.** It runs `build_site.R` from the research repo root.
   Any further arguments pass through: `--render-archive` re-renders the six
   archived reports with the pinned `research/libs/mixed-final` library;
   otherwise they stay frozen as rendered and only get the new navbar and the
   archive banner.
4. **Check freshness.** Every vignette page and every `research/docs/*.html`
   must have an mtime newer than the build start. A render that aborts
   part-way leaves the earlier pages fresh and the later ones stale, with no
   error in the summary line.

## The known failures and their remedies

| Error text | Cause | Remedy |
|---|---|---|
| `MISSING render: research/<study>/report/<name>.html` | A current report listed in `build_site.R` was never rendered on this clone. | `cd research && quarto render <study>/report/<name>.qmd` (it reads the study's `out/` data; see the study README), then rerun. |
| `WARNING: ... is older than ....qmd` | The report source changed after its last render. | Re-render before publishing, or the site ships the old text. |
| `archived report missing: docs/<name>.html (run with --render-archive)` | `research/docs/` lacks an archived page on this clone. | `build_site.sh research --render-archive`. It loads `research/libs/mixed-final`, and the archive reports need their old canonical tables `research/reports/benchmarks/tables/*.rds`. |
| `STALE <page>.html` under `research/docs/` with no build error | The page is not in `build_site.R`'s lists (an orphan), or the build stopped part-way. | An orphan: remove it in the research repo. Otherwise read the build output after the last fresh page. |
| `In _pkgdown.yml, N vignette(s) missing from index`, or an error naming an indexed article that does not exist | `_pkgdown.yml`'s `articles:` list and `vignettes/*.Rmd` disagree (the index has no ignore mechanism). | List exactly the vignettes that exist. A validation report is a Validation-menu link to the research site, not an article. |
| `pkgdown::build_site()` fails on `install = TRUE`, or `SpaNorm (>= 1.7.14)` not available | The package does not install: usually a stale `NAMESPACE`/`man/`, or the library path lacks SpaNorm 1.7.14. | `Rscript -e 'devtools::document()'`; check `SPIDE_LIB`. |
| A vignette chunk errors | The vignette calls a removed function or argument (the mixed-model API: `polishSpiDE`, `random =`, `fits()`, ...) or reads a legacy object. | Update the vignette to the 0.99.30 API; an old object is read only with `spiDEmixed::readSpiDE()`. |
| `research submodule not checked out` | `research/` is empty on this clone. | `git submodule update --init`. |

## After it passes

Push. Poll the GitHub Actions API no faster than every three minutes. The
cluster shell has no `GITHUB_PAT`, the unauthenticated limit is 60 requests
per hour, and a 30 s poll locks you out for an hour.

Do not commit `research/docs/` from the package repo. It belongs to the
research submodule and is committed there, with the submodule pointer
updated in the package afterwards.
