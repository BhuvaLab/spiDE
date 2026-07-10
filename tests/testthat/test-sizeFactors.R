test_that("computeSizeFactors writes a per-sample LS offset", {
  spe <- .toySPE()
  spe <- computeSizeFactors(spe, count = "nCount", area = "Area")

  expect_true("LS" %in% colnames(SummarizedExperiment::colData(spe)))
  expect_length(spe$LS, ncol(spe))
  expect_false(anyNA(spe$LS))

  # LS is constant within a sample (it is a sample-level offset)
  by_sample <- tapply(spe$LS, spe$sample_id, function(x) length(unique(round(x, 10))))
  expect_true(all(by_sample == 1))

  # number of distinct LS values equals number of samples
  expect_equal(length(unique(spe$LS)), length(unique(spe$sample_id)))
})

test_that("computeSizeFactors errors on missing count/area columns", {
  spe <- .toySPE()
  expect_error(computeSizeFactors(spe, count = "nope", area = "Area"), "not found")
  expect_error(computeSizeFactors(spe, count = "nCount", area = "nope"), "not found")
})

test_that(".computeLS centres each cell type to mean zero", {
  spe <- .toySPE()
  sf <- .computeLS(spe, "cell_type", "sample_id", "nCount", "Area")
  # each column (cell type) is centred across samples
  expect_true(all(abs(colMeans(sf, na.rm = TRUE)) < 1e-8))
})
