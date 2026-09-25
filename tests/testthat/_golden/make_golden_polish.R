# Golden outputs of polish + inference, captured BEFORE the polish machinery moved
# to SpaNorm (design/plans/2026-09-25-polish-to-spanorm.md). Regenerate ONLY on a
# deliberate numerical change, never to make the move pass.
#
# Run with the development SpaNorm library the move targets, e.g.:
#   R_LIBS=/scratch/user/uqdbhuva/Rlib_polishmove Rscript tests/testthat/_golden/make_golden_polish.R
devtools::load_all(quiet = TRUE)
saveRDS(golden_polish_fits(), "tests/testthat/_golden/golden_polish.rds")
