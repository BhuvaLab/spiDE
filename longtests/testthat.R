# Bioconductor Long Tests builder entry point. Runs weekly with a 6-hour
# budget, separate from the nightly check, so the numerically demanding
# mixed-effects assertions live here rather than in tests/.
library(testthat)
library(spiDE)
test_check("spiDE")
