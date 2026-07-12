# Precompile the one heavy spiDE vignette.
#
# Only `spiDE-simulation` is precomputed: it draws a large (12,000-gene) dataset
# live for its "what a realistic dataset looks like" figures, which is too slow
# for the Bioconductor build farm. We use the standard knitr `.Rmd.orig`
# precompile pattern: the live source lives in `spiDE-simulation.Rmd.orig`, which
# is knit here, once, into a static `spiDE-simulation.Rmd` with all outputs and
# figures baked in. The shipped `.Rmd` has no evaluatable chunks, so
# `R CMD build` renders it in seconds.
#
# The other vignettes (`spiDE`, `spiDE-model`, `spiDE-cauchy-vs-brown`,
# `spiDE-mixed-benchmark`) are ordinary live vignettes built at install time --
# they only load small summary tables or run a quick Monte-Carlo, so they do not
# need precompiling.
#
# Run manually from the package ROOT after editing spiDE-simulation.Rmd.orig,
# with the dev version of the package loaded so the baked results reflect it:
#
#   Rscript -e 'devtools::load_all(); source("vignettes/precompile.R")'
#
# Commit the regenerated `spiDE-simulation.Rmd` and its `-fig-*.png` images.

vignettes <- c(
  "spiDE-simulation"
)

local({
  old <- setwd("vignettes")
  on.exit(setwd(old), add = TRUE)
  for (v in vignettes) {
    message("Precompiling ", v, " ...")
    # per-vignette figure prefix so the committed PNGs don't collide
    knitr::opts_chunk$set(fig.path = paste0(v, "-fig-"))
    knitr::knit(
      input = paste0(v, ".Rmd.orig"),
      output = paste0(v, ".Rmd")
    )
  }
})
