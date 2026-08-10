# Two-stage estimator: the properties that justify its existence.
test_that("twoStageSpiDE returns a populated SpiDEResults", {
  spe <- buildNiches(spiDE:::.toySPE(n_genes = 30), sigma = 30, verbose = FALSE)
  r <- twoStageSpiDE(spe, condition = "condition", sigma = 30, min.cells = 10,
                     verbose = FALSE)
  expect_s4_class(r, "SpiDEResults")
  expect_true(is.data.frame(results(r)))
  expect_identical(r@sigma, 30)
  expect_length(r@fits, 0L)          # no per-bandwidth GLM fit exists
})

test_that("it accepts patient-level covariates that fitSpiDE rejects", {
  # This is the capability that motivates a separate estimator: checkSample()
  # refuses sample-constant covariates under random != "none" because they are
  # collinear with the per-sample random intercept. Two-stage analyses AT the
  # patient level, so they are ordinary covariates.
  spe <- buildNiches(spiDE:::.toySPE(n_genes = 30), sigma = 30, verbose = FALSE)
  expect_error(
    fitSpiDE(spe, condition = "condition", sigma = 30, covariates = "Age",
             verbose = FALSE),
    "constant within sample")
  expect_s4_class(
    twoStageSpiDE(spe, condition = "condition", sigma = 30, min.cells = 10,
                  patient.covariates = "Age", verbose = FALSE),
    "SpiDEResults")
})

test_that("a within-patient condition is refused", {
  # The contrast is patient-level by construction; silently averaging a
  # cell-varying condition would produce a meaningless slope contrast.
  spe <- buildNiches(spiDE:::.toySPE(n_genes = 20), sigma = 30, verbose = FALSE)
  set.seed(1)
  spe$condition <- sample(c("A", "B"), ncol(spe), replace = TRUE)
  expect_error(
    twoStageSpiDE(spe, condition = "condition", sigma = 30, verbose = FALSE),
    "varies within patient")
})

test_that("index and niche restriction shrinks the hypothesis space", {
  # Restriction is not cosmetic here: the full space buries the signal under
  # multiplicity, so the arguments must actually take effect.
  spe <- buildNiches(spiDE:::.toySPE(n_genes = 20), sigma = 30, verbose = FALSE)
  full <- twoStageSpiDE(spe, condition = "condition", sigma = 30, min.cells = 10,
                        fdr = 1, verbose = FALSE)
  restr <- twoStageSpiDE(spe, condition = "condition", sigma = 30, min.cells = 10,
                         index = "A", niche = "B", fdr = 1, verbose = FALSE)
  expect_lt(nrow(results(restr)), nrow(results(full)))
  expect_identical(unique(results(restr)$ct_index), "A")
  expect_identical(unique(results(restr)$ct_niche), "B")
})

test_that("patients below min.cells are dropped, not silently included", {
  spe <- buildNiches(spiDE:::.toySPE(n_genes = 20), sigma = 30, verbose = FALSE)
  strict <- twoStageSpiDE(spe, condition = "condition", sigma = 30,
                          min.cells = 10000L, fdr = 1, verbose = FALSE)
  expect_equal(nrow(results(strict)), 0L)   # nothing estimable, nothing invented
})
