#!/bin/sh
# Pre-render step (run by quarto from the site/ directory). Copies files that live outside the
# project into it, since Quarto only publishes files inside the project. All copies are
# git-ignored.
#
# SPIDEPY_DIR: checkout of https://github.com/BhuvaLab/spidepy holding docs/r-to-python.md,
# and docs/r-discrepancies.md. Default: a sibling of this repository
# (../spidepy relative to the repository root), which is where it lives locally. CI checks it
# out elsewhere and sets the variable.
set -e
SPIDEPY_DIR="${SPIDEPY_DIR:-../../spidepy}"   # relative to site/

cp ../man/figures/logo.png logo.png

mkdir -p _includes articles
cp ../NEWS.md _includes/NEWS.md

# Pages rendered from the package vignettes. Each entry is  <vignette>:<site page>.
# (Get started is not among them: site/get-started.qmd is written by hand, R and Python side by
# side, and follows vignettes/spiDE.Rmd and spidepy's docs/tutorial.ipynb.)
# The vignettes stay as they are; only the `output:` and `vignette:` front matter (BiocStyle
# and vignette-engine keys, meaningless to Quarto) is dropped from the copy.
# To change the article list edit ARTICLES here, the Articles menu in _quarto.yml and
# articles/index.qmd.
ARTICLES="spiDE-model:articles/spiDE-model spiDE-design:articles/spiDE-design spiDE-calibration:articles/spiDE-calibration"
for a in $ARTICLES; do
  src="../vignettes/${a%%:*}.Rmd"
  dst="${a##*:}.Rmd"
  awk '
    /^---$/ { fm++ ; print; next }
    fm == 1 && /^(output|vignette):/ { skip = 1; next }
    fm == 1 && skip && /^[ \t]/ { next }
    { skip = 0; print }
  ' "$src" > "$dst"
done

for f in docs/r-to-python.md docs/r-discrepancies.md; do
  if [ ! -f "$SPIDEPY_DIR/$f" ]; then
    echo "copy-assets.sh: $SPIDEPY_DIR/$f not found; set SPIDEPY_DIR to a spidepy checkout" >&2
    exit 1
  fi
done
cp "$SPIDEPY_DIR/docs/r-to-python.md" _includes/r-to-python.md
cp "$SPIDEPY_DIR/docs/r-discrepancies.md" _includes/r-discrepancies.md
# the two pages are shown on one site page, so the cross-link becomes an anchor
sed -i 's|(r-discrepancies.md)|(#differences-from-r-spide)|' _includes/r-to-python.md
