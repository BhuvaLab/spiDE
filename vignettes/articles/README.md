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
