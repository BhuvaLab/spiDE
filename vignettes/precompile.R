# Precompile the expensive spiDE vignettes.
#
# These three vignettes run genuinely expensive computation (Monte-Carlo
# calibration, timed PQL re-fits, bplapply parameter sweeps) that must NOT run on
# the Bioconductor build farm. We use the standard knitr `.Rmd.orig` precompile
# pattern: the live source lives in `<name>.Rmd.orig`, which is knit here, once,
# into a static `<name>.Rmd` with all outputs and figures baked in. The shipped
# `.Rmd` has no evaluatable chunks, so `R CMD build` renders it in seconds.
#
# Run manually from the package ROOT after editing any *.Rmd.orig, with the dev
# version of the package loaded so the baked results reflect it:
#
#   Rscript -e 'devtools::load_all(); source("vignettes/precompile.R")'
#
# Commit the regenerated `<name>.Rmd` files and their `<name>-fig-*.png` images.

vignettes <- c(
  "spiDE-cauchy-vs-brown",
  "spiDE-mixed-benchmark",
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
