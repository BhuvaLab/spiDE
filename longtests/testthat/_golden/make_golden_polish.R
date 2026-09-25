# Golden outputs of polish + inference, captured BEFORE the polish machinery moved
# to SpaNorm (design/plans/2026-09-25-polish-to-spanorm.md). Regenerate ONLY on a
# deliberate numerical change, never to make the move pass.
#
# This lives under longtests/ (see longtests/testthat/test-polish-golden.R and
# longtests/testthat/helper-golden.R): the four production-pipeline fits take
# ~4-5 min and tolerance 0 cannot be expected to hold across CI's other
# platforms/BLAS. devtools::load_all()'s automatic helper-sourcing only covers
# tests/testthat/helper*.R, so the longtests helper is sourced explicitly here.
#
# Run from the package root, with the development SpaNorm library the move
# targets, e.g.:
#   R_LIBS=/scratch/user/uqdbhuva/Rlib_polishmove Rscript longtests/testthat/_golden/make_golden_polish.R
devtools::load_all(quiet = TRUE)
source("longtests/testthat/helper-golden.R")
saveRDS(golden_polish_fits(), "longtests/testthat/_golden/golden_polish.rds")
