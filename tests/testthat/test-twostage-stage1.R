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
