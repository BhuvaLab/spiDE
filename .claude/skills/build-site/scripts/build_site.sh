#!/usr/bin/env bash
# Build the spiDE package site and/or the research site, in the required
# order, into scratch/expected destinations, and verify every report rendered.
#
#   build_site.sh [package|research|both]     (default: both)
#
# Env: SITE_DEST   pkgdown destination (default: $TMPDIR/spiDE-site or /tmp)
#      SPIDE_ROOT  package root (default: git toplevel of the cwd)
set -euo pipefail
what="${1:-both}"
ROOT="${SPIDE_ROOT:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"
cd "$ROOT"
[ -f DESCRIPTION ] && grep -q '^Package: spiDE' DESCRIPTION || { echo "not the spiDE package root: $ROOT" >&2; exit 1; }
start=$(date +%s)
fail=0

if [ "$what" = package ] || [ "$what" = both ]; then
  echo "== [1/2] sync articles from the research submodule (index-driven)"
  Rscript vignettes/articles/pkgdown-sync-articles.R
  dest="${SITE_DEST:-${TMPDIR:-/tmp}/spiDE-site}"
  echo "== [1/2] pkgdown::build_site -> $dest"
  Rscript -e "pkgdown::build_site(override = list(destination = '$dest'), install = TRUE, lazy = FALSE, preview = FALSE)"
  echo "-- articles rendered:"
  for rmd in vignettes/articles/spiDE-*.Rmd; do
    html="$dest/articles/$(basename "${rmd%.Rmd}").html"
    if [ -f "$html" ] && [ "$(stat -c %Y "$html")" -ge "$start" ]; then
      echo "   ok    $(basename "$html")"
    else
      echo "   STALE $(basename "$html") (missing or older than the build start)"; fail=1
    fi
  done
fi

if [ "$what" = research ] || [ "$what" = both ]; then
  echo "== [2/2] research site: build_site.R -> research/docs"
  [ -f research/reports/benchmarks/build_site.R ] || { echo "research submodule not checked out" >&2; exit 1; }
  ( cd research && Rscript reports/benchmarks/build_site.R )
  echo "-- reports rendered:"
  for rmd in research/reports/benchmarks/spiDE-*.Rmd; do
    html="research/docs/$(basename "${rmd%.Rmd}").html"
    if [ -f "$html" ] && [ "$(stat -c %Y "$html")" -ge "$start" ]; then
      echo "   ok    $(basename "$html")"
    else
      echo "   STALE $(basename "$html") (missing or older than the build start)"; fail=1
    fi
  done
fi

if [ "$fail" -ne 0 ]; then
  echo "== FAILED: at least one report did not render this run; read the build output above the last fresh page" >&2
  exit 1
fi
echo "== done in $(( $(date +%s) - start )) s"
