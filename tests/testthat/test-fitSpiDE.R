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

test_that("depth = 'nonlinear' puts a 3-df spline of log depth in the design", {
  spe <- buildNiches(.toySPE(n_samples = 8, n_per = 120, n_genes = 6, seed = 3), sigma = 20,
                     verbose = FALSE)
  fit <- fitSpiDE(spe, index = "A", sigma = 20, depth = "nonlinear", verbose = FALSE)
  expect_identical(fit@params$depth, "nonlinear")
  x <- fit@index$A
  expect_true(any(is.finite(x$beta)))
  ct <- as.character(spe$cell_type)
  ik <- spiDE:::.indexCells(ct, as.character(spe$sample_id), "A", 10L)
  lib <- log(Matrix::colSums(SummarizedExperiment::assay(spe, "counts")))
  covk <- spiDE:::.indexCovariates(cbind(loglib = lib[ik]), "nonlinear")
  expect_equal(colnames(covk), paste0("loglib_ns", 1:3))
  expect_equal(unname(colMeans(covk)), rep(0, 3), tolerance = 1e-12)
})

# Absorbed to the precision of the 1e-3 ridge on the patient intercepts (the
# shifted intercepts are penalised slightly differently): ~1e-4 on the toy.
test_that("a per-gene constant offset is absorbed by the patient intercepts", {
  spe <- buildNiches(.toySPE(n_samples = 8, n_per = 120, n_genes = 6, seed = 3), sigma = 20,
                     verbose = FALSE)
  O <- matrix(seq(-1, 1, length.out = nrow(spe)), nrow(spe), ncol(spe),
              dimnames = dimnames(spe))
  SummarizedExperiment::assay(spe, "shift") <- O
  f0 <- fitSpiDE(spe, index = "A", sigma = 20, verbose = FALSE)
  f1 <- fitSpiDE(spe, index = "A", sigma = 20, offset = "shift", verbose = FALSE)
  expect_equal(f1@index$A$beta, f0@index$A$beta, tolerance = 1e-3)
  expect_equal(f1@index$A$v_tile, f0@index$A$v_tile, tolerance = 1e-3)
  s0 <- fitSpiDE(spe, index = "A", sigma = 20, engine = "sandwich", verbose = FALSE)
  s1 <- fitSpiDE(spe, index = "A", sigma = 20, engine = "sandwich", offset = "shift", verbose = FALSE)
  expect_equal(s1@index$A$coef$estimate, s0@index$A$coef$estimate, tolerance = 1e-3)
  expect_equal(s1@index$A$coef$se, s0@index$A$coef$se, tolerance = 1e-3)
  expect_error(fitSpiDE(spe, index = "A", sigma = 20, offset = "nope", verbose = FALSE),
               "must name an assay")
})

test_that("a log-depth offset replaces the depth covariate", {
  spe <- buildNiches(.toySPE(n_samples = 8, n_per = 120, n_genes = 6, seed = 3), sigma = 20,
                     verbose = FALSE)
  ll <- log(pmax(Matrix::colSums(SummarizedExperiment::assay(spe, "counts")), 1))
  SummarizedExperiment::assay(spe, "logdepth") <- matrix(ll, nrow(spe), ncol(spe), byrow = TRUE,
                                                         dimnames = dimnames(spe))
  f <- fitSpiDE(spe, index = "A", sigma = 20, depth = "none", offset = "logdepth", verbose = FALSE)
  r <- testSpiDE(f)
  expect_true(all(is.finite(r@table$t[!is.na(r@table$t)])))
  expect_identical(f@params$offset, "logdepth")
})

test_that("genes restricts the tested genes but not the library size", {
  spe <- buildNiches(.toySPE(n_samples = 8, n_per = 120, n_genes = 6, seed = 3), sigma = 20,
                     verbose = FALSE)
  f <- fitSpiDE(spe, index = "A", sigma = 20, genes = c("G1", "G2"), verbose = FALSE)
  expect_setequal(f@index$A$genes, c("G1", "G2"))
  full <- fitSpiDE(spe, index = "A", sigma = 20, verbose = FALSE)
  expect_equal(f@index$A$beta["G1", , ], full@index$A$beta["G1", , ], tolerance = 1e-8)
})
