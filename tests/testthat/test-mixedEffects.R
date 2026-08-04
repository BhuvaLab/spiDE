# Tests for the mixed-effects (random-effects via ridge) correction for
# cell-level pseudo-replication. See R/fitSpiDE.R (.fitNBmixed) and R/inference.R.

fixture <- function(nm) {
  p <- system.file("extdata", "testfits", nm, package = "spiDE")
  if (!nzchar(p)) skip(paste("fixture not installed:", nm))
  readRDS(p)
}

test_that("random='none' has no RE slots; mixed fits populate them", {
  # contract only -- slot presence and nullity -- so it runs against fixtures.
  # It used to do two live fits, one of them random='slope', which is the most
  # expensive thing the short suite did.
  fa <- fixture("fit_toyspe_none.rds")
  expect_null(fa@re_group)
  expect_null(fa@tau2)
  expect_null(fa@penalty)
  expect_null(fa@df)

  fb <- fixture("fit_toyspe_slope.rds")
  expect_true("Random" %in% levels(fb@covtype))
  expect_true(all(c("SampleInt", "SampleSlope") %in% fb@re_group))
  expect_equal(length(fb@penalty), ncol(fb@W))
  expect_length(fb@df, 1L)
  expect_true(is.finite(fb@df))
})

test_that("random-effect design targets the response-related effects", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  des <- spiDE:::.buildNicheDesign(spe, "condition", 20, random = "slope")
  reg <- des$re_group
  n_samples <- length(unique(spe$sample_id))
  # one random intercept per sample (counterpart of the response main effect)
  expect_equal(sum(reg == "SampleInt", na.rm = TRUE), n_samples)
  # random slopes = samples x (CellType:niche base terms) i.e. the non-response
  # bases of the ResponseNiche terms
  n_niche_base <- sum(des$covtype == "Niche")
  expect_equal(sum(reg == "SampleSlope", na.rm = TRUE), n_samples * n_niche_base)
})

test_that("checkSample rejects within-sample condition and sample-level covariates", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  # a covariate constant within sample is confounded with the random intercept
  spe$Age <- ave(rnorm(ncol(spe)), spe$sample_id)
  expect_error(
    fitSpiDE(spe, "condition", sigma = 20, covariates = "Age",
             random = "intercept", verbose = FALSE),
    "constant within sample"
  )
  # a condition that varies within a sample is not a patient-level label
  spe$bad <- rep(c("x", "y"), length.out = ncol(spe))
  expect_error(
    fitSpiDE(spe, "bad", sigma = 20, random = "intercept", verbose = FALSE),
    "varies within sample"
  )
})

test_that(".stratifiedCellIdx honours the proportion, floor and cap per stratum", {
  set.seed(1)
  # three samples x two cell types; stratum sizes span the floor and the prop
  cell_type <- rep(c("A", "B"), times = c(1500, 40))
  sample <- rep(c("S1", "S2", "S3"), length.out = length(cell_type))
  cell_type <- sample(cell_type) # shuffle so strata interleave
  idx <- spiDE:::.stratifiedCellIdx(cell_type, sample, prop = 0.1,
                                    min.cells = 100L)
  expect_type(idx, "logical")
  expect_length(idx, length(cell_type))

  strata <- interaction(sample, cell_type, drop = TRUE)
  for (s in levels(strata)) {
    n <- sum(strata == s)
    kept <- sum(idx[strata == s])
    expect_equal(kept, min(n, max(ceiling(0.1 * n), 100L)))
  }
  # every stratum keeps at least one cell (samples stay represented)
  expect_true(all(tapply(idx, strata, any)))

  # prop >= 1 keeps everything
  expect_true(all(spiDE:::.stratifiedCellIdx(cell_type, sample, prop = 1)))
})

test_that("re.prop is validated", {
  # the cheap half: argument validation fails fast, no fit is performed. The
  # expensive half -- that a re.prop below the per-stratum floor leaves tau2 and
  # alpha unchanged -- asserts what fitSpiDE COMPUTES, so it moved to
  # longtests/test-mixed-numerics.R rather than being frozen into a fixture.
  spe <- buildNiches(.toySPE(), sigma = 20)
  expect_error(
    fitSpiDE(spe, "condition", sigma = 20, random = "intercept",
             re.prop = 1.5, verbose = FALSE),
    "re.prop"
  )
})

# The two statistical tests that were here -- tau2 recovery and the
# anti-conservative deflation -- moved to longtests/test-mixed-numerics.R.
# Between them they cost ~1,200 s, against Bioconductor's 10-minute budget for
# the whole check. They could not be precomputed: both assert what fitSpiDE()
# COMPUTES, so freezing them would let a regression in the fitting code pass.
# The Long Tests builder runs them weekly with a 6-hour budget.
