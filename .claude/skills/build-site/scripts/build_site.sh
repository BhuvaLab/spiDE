#!/usr/bin/env bash
# Build the spiDE package site (pkgdown) and/or the research site, into
# scratch/expected destinations, and verify every page was written this run.
#
#   build_site.sh [package|research|both] [build_site.R args, e.g. --render-archive]
#
# Env: SITE_DEST   pkgdown destination (default: $TMPDIR/spiDE-site or /tmp)
#      SPIDE_ROOT  package root (default: git toplevel of the cwd)
#      SPIDE_LIB   library prepended for the package build (default:
#                  research/libs/simplify, which holds SpaNorm >= 1.7.14; the
#                  default user library holds an older SpaNorm)
set -euo pipefail
what="${1:-both}"; shift || true
ROOT="${SPIDE_ROOT:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"
cd "$ROOT"
[ -f DESCRIPTION ] && grep -q '^Package: spiDE' DESCRIPTION || { echo "not the spiDE package root: $ROOT" >&2; exit 1; }
start=$(date +%s)
fail=0
fresh() {  # fresh <file> <label>: ok if it exists and was written this run
  if [ -f "$1" ] && [ "$(stat -c %Y "$1")" -ge "$start" ]; then
    echo "   ok    $2"
  else
    echo "   STALE $2 (missing or older than the build start)"; fail=1
  fi
}

if [ "$what" = package ] || [ "$what" = both ]; then
  dest="${SITE_DEST:-${TMPDIR:-/tmp}/spiDE-site}"
  lib="${SPIDE_LIB:-$ROOT/research/libs/simplify}"
  [ -d "$lib" ] && export R_LIBS="$lib${R_LIBS:+:$R_LIBS}"
  echo "== [package] pkgdown::build_site -> $dest (R_LIBS=${R_LIBS:-unset})"
  Rscript -e "pkgdown::build_site(override = list(destination = '$dest'), install = TRUE, lazy = FALSE, preview = FALSE)"
  echo "-- vignettes rendered:"
  for rmd in vignettes/*.Rmd; do
    fresh "$dest/articles/$(basename "${rmd%.Rmd}").html" "$(basename "${rmd%.Rmd}").html"
  done
  fresh "$dest/index.html" "index.html"
fi

if [ "$what" = research ] || [ "$what" = both ]; then
  [ -f research/reports/benchmarks/build_site.R ] || { echo "research submodule not checked out (git submodule update --init)" >&2; exit 1; }
  # current reports are rendered where they live and only copied by build_site.R:
  # warn when a rendered page is older than its source (it would publish a stale render)
  for src in $(grep -o 'src = "[^"]*"' research/reports/benchmarks/build_site.R | sed 's/src = "//; s/"$//'); do
    html="research/$src"; qmd="${html%.html}.qmd"
    if [ ! -f "$html" ]; then
      echo "   MISSING render: $html -- render it first (quarto render ${qmd#research/}, from research/)"; fail=1
    elif [ -f "$qmd" ] && [ "$qmd" -nt "$html" ]; then
      echo "   WARNING: $html is older than $qmd -- re-render before publishing"
    fi
  done
  echo "== [research] build_site.R $* -> research/docs"
  ( cd research && Rscript reports/benchmarks/build_site.R "$@" )
  echo "-- pages written:"
  for html in research/docs/*.html; do fresh "$html" "$(basename "$html")"; done
fi

if [ "$fail" -ne 0 ]; then
  echo "== FAILED: at least one page was not written this run; read the build output above" >&2
  exit 1
fi
echo "== done in $(( $(date +%s) - start )) s"
