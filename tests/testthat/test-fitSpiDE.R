# fitSpiDE(): input checks and the shape of the fit. The engines' numerics are
# in test-engine.R and test-sandwich.R, recovery in test-spiDE-e2e.R.

test_that("fitSpiDE errors when niches are missing", {
  spe <- .toySPE()
  expect_error(fitSpiDE(spe, condition = "condition"), "buildNiches|niche")
})

test_that("several bandwidths need an explicit sigma", {
  spe <- buildNiches(.toySPE(), sigma = c(10, 20), verbose = FALSE)
  expect_error(fitSpiDE(spe, index = "A", min.patients = 6, verbose = FALSE), "sigma")
})

test_that("non-integer counts are refused before any fitting", {
  spe <- buildNiches(.toySPE(), sigma = 20, verbose = FALSE)
  SummarizedExperiment::assay(spe, "counts") <- SummarizedExperiment::assay(spe, "counts") + 0.5
  expect_error(fitSpiDE(spe, index = "A", sigma = 20, min.patients = 6, verbose = FALSE),
               "integer counts")
})

test_that("a covariate with non-finite values is refused by name", {
  spe <- buildNiches(.toySPE(), sigma = 20, verbose = FALSE)
  SummarizedExperiment::colData(spe)$bad <- log(c(0, runif(ncol(spe) - 1)))
  expect_error(fitSpiDE(spe, sigma = 20, covariates = "bad", index = "A", min.patients = 6,
                        verbose = FALSE), "bad")
})

test_that("a patient-level covariate is refused: the patient intercepts absorb it", {
  spe <- buildNiches(.toySPE(), sigma = 20, verbose = FALSE)
  expect_error(fitSpiDE(spe, sigma = 20, covariates = "Age", index = "A", min.patients = 6,
                        verbose = FALSE), "constant within sample")
})

test_that("the fit records patient-level columns, so a condition can be named later", {
  spe <- buildNiches(.toySPE(), sigma = 20, verbose = FALSE)
  fit <- fitSpiDE(spe, index = "A", sigma = 20, min.patients = 6, verbose = FALSE)
  expect_true(all(c("patient", "ncells", "condition") %in% colnames(fit@patients)))
  expect_length(fit@condition, 0)
  res <- testSpiDE(fit, condition = "condition", procedure = "all")
  expect_identical(res@contrast, "Responder - Non-responder")
  expect_setequal(unique(res@table$test), c("pooled", "condition"))
})
