#!/bin/sh
# Pre-render step (run by quarto from the site/ directory). Copies files that live outside the
# project into it, since Quarto only publishes files inside the project. All copies are
# git-ignored.
set -e

cp ../man/figures/logo.png logo.png

mkdir -p _includes articles
cp ../NEWS.md _includes/NEWS.md

# Pages rendered from the package vignettes. Each entry is  <vignette>:<site page>.
# The vignettes stay as they are; only the `output:` and `vignette:` front matter (BiocStyle
# and vignette-engine keys, meaningless to Quarto) is dropped from the copy.
# To change the article list edit ARTICLES here, the Articles menu in _quarto.yml and
# articles/index.qmd.
ARTICLES="spiDE:get-started spiDE-model:articles/spiDE-model spiDE-design:articles/spiDE-design spiDE-calibration:articles/spiDE-calibration"
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
