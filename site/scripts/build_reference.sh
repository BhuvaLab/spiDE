#!/bin/sh
# Build the generated R reference pages (git-ignored) from man/*.Rd via altdoc into
# reference/r/. Run from anywhere:
#   bash site/scripts/build_reference.sh && quarto render site
set -e
cd "$(dirname "$0")/.."

# altdoc goes into the project-local library, never the user's
mkdir -p _rlib
# A FIXED altdoc version: the Rd conversion uses its private altdoc:::.rd2qmd.
ALTDOC_VERSION=0.7.3
export ALTDOC_VERSION
Rscript -e 'v <- Sys.getenv("ALTDOC_VERSION"); ok <- requireNamespace("altdoc", lib.loc = "_rlib", quietly = TRUE) && as.character(packageVersion("altdoc", lib.loc = "_rlib")) == v; if (!ok) { if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org"); remotes::install_version("altdoc", version = v, lib = "_rlib", repos = "https://cloud.r-project.org", upgrade = "never") }'
Rscript scripts/build-reference-r.R

# Quarto resolves {{< include >}} before its pre-render step runs, so copy the included
# files (NEWS.md) here too.
sh scripts/copy-assets.sh
