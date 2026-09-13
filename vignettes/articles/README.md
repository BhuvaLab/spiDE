# Validation-report articles (generated — do not edit here)

The `.Rmd` files and `tables/` in this directory are BUILD ARTEFACTS, synced
from the research submodule by `pkgdown-sync-articles.R` (run automatically in
the pkgdown CI workflow; run it manually before a local `pkgdown::build_site()`).
The canonical sources live in `research/reports/benchmarks/` — edit them there.

Why copies and not symlinks: pkgdown's article discovery lists files by type
and does not see symlinked `.Rmd`s (measured: a symlinked article is absent
from `as_pkgdown()$vignettes`). CI re-syncs on every build, so the published
site cannot drift from the submodule; only a stale local preview can, which is
why this README and the sync script exist.

Everything here except this README and the sync script is gitignored, and the
whole directory is `.Rbuildignore`d so the package tarball is untouched.

## If `pkgdown::build_site()` stops with "vignette missing from index"

```
Error in `build_articles_index()`:
! In _pkgdown.yml, 1 vignette missing from index: "articles/spiDE-twostage-benchmark".
```

pkgdown insists that every `.Rmd` it finds under `vignettes/articles/` is
listed in `_pkgdown.yml`, and it has no way to ignore one. The file it names is
a STALE COPY from a sync made before that report was retired from the package
site (the two-stage benchmark left it in 0.99.19). Because the copies are
gitignored, a `git pull` can neither remove nor replace it, so it survives on
every clone that synced earlier. Run the sync once on that clone and build
again:

```
Rscript vignettes/articles/pkgdown-sync-articles.R   # removes copies not in the index
Rscript -e 'pkgdown::build_site()'
```

The sync reads the article list from `_pkgdown.yml`, copies exactly those
reports and deletes any other `spiDE-*.Rmd` here, so the index and the
directory agree by construction. The reverse error ("in index but not on
disk") means the yaml lists a report the submodule no longer has: fix the
yaml. CI runs the sync before every build, which is why the published site
never shows either error.
