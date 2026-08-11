test_that(".spanormComponents extracts and validates the stored fit", {
  spe <- toy_spanorm_spe()
  comp <- spiDE:::.spanormComponents(spe)
  expect_named(comp, c("alpha", "gmean", "psi", "W", "bio"))
  expect_identical(dim(comp$W), c(ncol(spe), 3L))
  expect_identical(comp$bio, c(FALSE, TRUE, TRUE))
  expect_length(comp$psi, nrow(spe))
})

test_that(".spanormComponents errors usefully without a fit or on mismatch", {
  spe <- buildNiches(spiDE:::.toySPE(n_genes = 20), sigma = 30, verbose = FALSE)
  expect_error(spiDE:::.spanormComponents(spe), "SpaNorm::SpaNorm")
  spe2 <- toy_spanorm_spe()
  expect_error(spiDE:::.spanormComponents(spe2[, 1:10]), "does not match")
})

test_that(".stage1Epsilon linearises log counts minus the LS effect", {
  spe <- toy_spanorm_spe()
  comp <- spiDE:::.spanormComponents(spe)
  cells <- 1:50
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  se <- spiDE:::.stage1Epsilon(Y, comp, cells, epsilon = "addback")
  expect_identical(dim(se$eps), c(nrow(Y), 50L))
  expect_identical(dim(se$w), dim(se$eps))
  expect_true(all(se$w > 0 & is.finite(se$eps)))
  # at y = 0 the addback response is eta_bio - 1 exactly (z = -1)
  Wc <- comp$W[cells, , drop = FALSE]
  eta_bio <- comp$gmean +
    tcrossprod(comp$alpha[, comp$bio, drop = FALSE],
               Wc[, comp$bio, drop = FALSE])
  zero <- Y[, cells] == 0
  expect_equal(se$eps[zero], (eta_bio - 1)[zero])
  # "residual" drops the biology term
  sr <- spiDE:::.stage1Epsilon(Y, comp, cells, epsilon = "residual")
  expect_equal(se$eps - sr$eps, eta_bio, tolerance = 1e-12)
})
