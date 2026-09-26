# Safety net for moving the polish machinery to SpaNorm
# (design/plans/2026-09-25-polish-to-spanorm.md): pins the production
# fit -> polish -> test outputs at tolerance 0 against a golden snapshot
# captured with the pre-move code. Lives in longtests/ (see
# longtests/testthat/helper-golden.R's header) because the four fits take
# ~4-5 min.
#
# OPT-IN, not run by the Bioconductor long-test builder: `.BBSoptions` sets
# `RunLongTests: TRUE`, so everything under longtests/ otherwise runs weekly
# on whatever BLAS/platform the builder has, and a tolerance-0 comparison
# (fitNB's IRLS, the Newton polish, the Cholesky-based covariance -- all
# BLAS-heavy) cannot be expected to hold across BLAS kernels, let alone
# platforms: the same node and OpenBLAS version, only the CPU kernel dial
# changed (`OPENBLAS_CORETYPE=Prescott`), already moves `nbGramBatch` by up
# to 2.3e-12 and `invert_mat` by up to 2.2e-19. This is therefore a LOCAL
# move gate at tolerance 0 -- reproduce the pre-move numbers exactly on this
# machine's BLAS -- not a cross-platform correctness check, and it is gated
# on an explicit opt-in so it never runs unattended on the builder or in CI.
# Run it with:
#   SPIDE_RUN_GOLDEN=true Rscript -e 'devtools::load_all(); testthat::test_file("longtests/testthat/test-polish-golden.R")'
#
# No skip on a SpaNorm/spiDE version mismatch: Task 5 (the SpaNorm-side move)
# must be able to compare the golden across SpaNorm versions, so a mismatch
# has to fail loudly, not be swallowed by a skip. Every expectation instead
# carries an `info` string naming the golden's and the running session's
# SpaNorm/spiDE versions, so a failure reads as a diagnosis (stale golden vs.
# real regression) rather than a bare "not identical".
#
# The per-bandwidth list (`res@fits`) is iterated by POSITION, not by name:
# polishSpiDE() re-derives it with `lapply(seq_along(object@fits), ...)`
# (R/polish.R), which drops the `Niche<sigma>` names fitSpiDE() set, so
# `names(g[[nm]]$fits)` is NULL for every config here -- `for (bw in
# names(...))` would silently iterate zero times and this test would check
# nothing but the (always-empty on this fixture) `results` table. Matching by
# position is exact here because both sides build the list in the same
# `sigma` order from the same call. (The name-dropping itself looks like an
# unrelated production bug in polishSpiDE(), not fixed here -- out of this
# task's scope, flagged in the report.)
skip_if_not(identical(Sys.getenv("SPIDE_RUN_GOLDEN"), "true"),
            "golden gate is opt-in: set SPIDE_RUN_GOLDEN=true to run it")

test_that("polish + inference reproduce the pre-move golden outputs exactly", {
  g <- readRDS(test_path("_golden", "golden_polish.rds"))
  now <- golden_polish_fits()
  ver <- sprintf(
    "golden built with SpaNorm %s / spiDE %s; running SpaNorm %s / spiDE %s",
    g$meta$SpaNorm, g$meta$spiDE,
    as.character(utils::packageVersion("SpaNorm")),
    as.character(utils::packageVersion("spiDE"))
  )
  for (nm in setdiff(names(g), "meta")) {
    expect_equal(length(now[[nm]]$fits), length(g[[nm]]$fits),
                 info = paste(nm, ver))
    for (bw in seq_along(g[[nm]]$fits)) {
      for (s in names(g[[nm]]$fits[[bw]])) {
        expect_equal(now[[nm]]$fits[[bw]][[s]], g[[nm]]$fits[[bw]][[s]],
                     tolerance = 0, info = paste(nm, bw, s, ver))
      }
    }
    expect_equal(now[[nm]]$results, g[[nm]]$results, tolerance = 0,
                 info = paste(nm, ver))
  }
})
