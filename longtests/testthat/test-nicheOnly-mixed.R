# Condition-free (niche-only) mixed-effects fits, with LIVE fits.
#
# These sit here rather than in tests/ for the same reason as
# test-mixed-numerics.R: each is a real PQL fit, and between them they cost far
# more than Bioconductor's 10-minute budget for the nightly check. The Long
# Tests builder runs weekly with a 6-hour budget.
#
# .toyNiche() and .toyClustered() are internal, so they are reached via
# spiDE:::.

# --- mixed effects in niche mode -------------------------------------------
# .toyNiche()'s Age is constant within sample, so it must NOT be passed as a
# covariate here: the per-sample random intercept absorbs it and checkSample()
# rejects it by design.
spe_mix <- buildNiches(spiDE:::.toyNiche(n_samples = 6), sigma = 20)

test_that("condition-free mixed fits run with a random intercept", {
  set.seed(1)
  fi <- fitSpiDE(spe_mix, condition = NULL, sigma = 20, random = "intercept",
                 verbose = FALSE)
  expect_identical(fi@mode, "niche")
  f <- fits(fi)[[1]]
  expect_false(is.null(f@tau2))
  expect_false(is.null(f@penalty))
  expect_true(all(f@df > 0))
  # satterthwaite (the default) gives one df per tested column
  expect_equal(length(f@df), sum(spiDE:::.testedCols(f@covtype, "niche")))
})

test_that("condition-free mixed fits run with random slopes", {
  set.seed(1)
  fs <- fitSpiDE(spe_mix, condition = NULL, sigma = 20, random = "slope",
                 verbose = FALSE)
  f <- fits(fs)[[1]]
  # the slopes sit on the tested CellType:niche columns
  expect_true(any(f@re_group == "SampleSlope", na.rm = TRUE))
  expect_true(all(f@df > 0))
})

test_that("df.method='between' is cell-level under random='intercept'", {
  # No condition contrast exists, and the CellType:niche slope is identified by
  # niche density varying cell-to-cell INSIDE each sample -- so the replication
  # unit is the cell, not the patient.
  set.seed(1)
  fi <- fitSpiDE(spe_mix, condition = NULL, sigma = 20, random = "intercept",
                 df.method = "between", verbose = FALSE)
  f <- fits(fi)[[1]]
  p_fixed <- sum(is.na(f@re_group))
  expect_equal(length(f@df), 1L)
  expect_equal(f@df, max(f@ncells - p_fixed, 1))
})

test_that("df.method='between' is S - 1 under random='slope'", {
  # Per-sample random slopes on the tested columns move the contrast back into
  # the between-sample stratum.
  set.seed(1)
  fs <- fitSpiDE(spe_mix, condition = NULL, sigma = 20, random = "slope",
                 df.method = "between", verbose = FALSE)
  expect_equal(fits(fs)[[1]]@df, 5)   # 6 samples - 1
})

test_that("condition mode still uses S - 2 under df.method='between'", {
  spe_c <- buildNiches(spiDE:::.toySPE(n_samples = 6), sigma = 20)
  set.seed(1)
  fc <- fitSpiDE(spe_c, condition = "condition", sigma = 20,
                 random = "intercept", df.method = "between", verbose = FALSE)
  expect_equal(fits(fc)[[1]]@df, 4)   # 6 samples - 2
})

test_that("random effects sharply reduce false calls on null clustered data", {
  # .toyClustered() plants a per-sample random INTERCEPT and no niche effect.
  # A sample-constant shift cannot create a within-sample niche slope, so every
  # CellType:niche coefficient is null by construction.
  #
  # What this pins down is that the random-intercept correction does real work
  # in niche mode: on this fixture the fixed-effects fit makes 37 calls and the
  # random-intercept fit makes 5. It deliberately does NOT assert nominal
  # type-I control, because niche mode does not achieve it here and asserting
  # otherwise would be a false comfort:
  #
  #   random = "none"      37 calls
  #   random = "intercept"  5 calls
  #   random = "slope"      5 calls  (tau2_SampleSlope collapses to its floor,
  #                                   correctly -- this fixture has no
  #                                   between-sample slope heterogeneity)
  #
  # The residual is unmodelled SPATIAL autocorrelation, not between-sample
  # pseudo-replication: the niche covariate varies smoothly in space, so
  # neighbouring cells are not independent replicates of a within-sample slope,
  # and Satterthwaite duly hands these columns a near-cell-level df (~624 of
  # 640 cells). spiDE does not model spatial autocorrelation. See the
  # limitation noted in ?fitSpiDE.
  spe_null <- buildNiches(spiDE:::.toyClustered(n_samples = 8, n_genes = 30),
                          sigma = 20)
  set.seed(1)
  fixed <- spiDE(spe_null, condition = NULL, sigma = 20, random = "none",
                 fdr = 0.05, verbose = FALSE)
  set.seed(1)
  mixed <- spiDE(spe_null, condition = NULL, sigma = 20, random = "intercept",
                 fdr = 0.05, verbose = FALSE)

  n_fixed <- nrow(results(fixed))
  n_mixed <- nrow(results(mixed))
  expect_gt(n_fixed, 20L)          # the uncorrected fit is badly inflated
  expect_lt(n_mixed, n_fixed / 3)  # the correction removes most of it
  expect_lt(n_mixed, 10L)          # and what remains is a handful, not a flood
})
