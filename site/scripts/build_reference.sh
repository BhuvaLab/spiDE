#!/bin/sh
# Build the generated reference pages (git-ignored): R pages from man/*.Rd via altdoc into
# reference/r/, Python pages via quartodoc into reference/python/. Examples of both are executed
# at render time. Run from anywhere:
#   bash site/scripts/build_reference.sh && quarto render site
# Local settings (QUARTO_PYTHON, SPIDEPY_DIR, R_LIBS, ...) come from _environment.local;
# PYTHON overrides the interpreter.
set -e
cd "$(dirname "$0")/.."

if [ -f _environment.local ]; then
  set -a
  # shellcheck disable=SC1091
  . ./_environment.local 2>/dev/null || true
  set +a
fi
PYTHON="${PYTHON:-${QUARTO_PYTHON:-python}}"

# Python tools: pip install -r requirements.txt (quartodoc pinned; needs griffe < 2)
# altdoc goes into the project-local library, never the user's
mkdir -p _rlib
# A FIXED altdoc version: the Rd conversion uses its private altdoc:::.rd2qmd.
ALTDOC_VERSION=0.7.3
export ALTDOC_VERSION
Rscript -e 'v <- Sys.getenv("ALTDOC_VERSION"); ok <- requireNamespace("altdoc", lib.loc = "_rlib", quietly = TRUE) && as.character(packageVersion("altdoc", lib.loc = "_rlib")) == v; if (!ok) { if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org"); remotes::install_version("altdoc", version = v, lib = "_rlib", repos = "https://cloud.r-project.org", upgrade = "never") }'
Rscript scripts/build-reference-r.R

# quartodoc needs griffe < 2 (quartodoc 0.11 does not support griffe 2): see requirements.txt.
# The script builds the pages and turns each docstring's Examples (spidepy's own, on toy_spide())
# into executed chunks.
"$PYTHON" scripts/build_python_reference.py

# Quarto resolves {{< include >}} before its pre-render step runs, so copy the included
# files (spidepy docs, NEWS.md) here too. SPIDEPY_DIR: see copy-assets.sh.
sh scripts/copy-assets.sh

# every R export and every spidepy.__all__ name must have a generated page
"$PYTHON" scripts/check_reference_coverage.py
