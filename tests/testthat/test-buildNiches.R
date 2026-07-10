test_that("buildNiches adds one reducedDim per bandwidth with correct shape", {
  spe <- .toySPE()
  spe <- buildNiches(spe, sigma = c(10, 20))

  expect_true(all(c("Niche10", "Niche20") %in%
    SingleCellExperiment::reducedDimNames(spe)))

  n10 <- SingleCellExperiment::reducedDim(spe, "Niche10")
  expect_equal(nrow(n10), ncol(spe))
  expect_setequal(colnames(n10), sort(unique(spe$cell_type)))
})

test_that("niche densities are non-negative and free of NA", {
  spe <- .toySPE()
  spe <- buildNiches(spe, sigma = 20)
  n <- SingleCellExperiment::reducedDim(spe, "Niche20")
  expect_true(all(n >= 0))
  expect_false(anyNA(n))
})

test_that("niche density recovers a planted spatial gradient", {
  # B cells are concentrated at high x, so the B-niche should track x
  spe <- .toySPE()
  spe <- buildNiches(spe, sigma = 15)
  n <- SingleCellExperiment::reducedDim(spe, "Niche15")
  x <- SpatialExperiment::spatialCoords(spe)[, 1]
  expect_gt(cor(n[, "B"], x), 0.4)
})

test_that("rows align with cells regardless of sample interleaving", {
  spe <- .toySPE()
  # shuffle cells so samples are interleaved, not grouped
  set.seed(99)
  spe <- spe[, sample(ncol(spe))]
  spe <- buildNiches(spe, sigma = 15)
  n <- SingleCellExperiment::reducedDim(spe, "Niche15")
  # each cell's niche must be computed within its own sample: recompute for one
  # sample and compare
  sid <- spe$sample_id[1]
  sub <- spe[, spe$sample_id == sid]
  ref <- .effectiveNiche(sub, sigma = 15)
  ref <- .fillMissingDims(ref, colnames(n), "column")
  ref[is.na(ref)] <- 0
  got <- n[spe$sample_id == sid, colnames(ref), drop = FALSE]
  expect_equal(unname(got), unname(ref), tolerance = 1e-8)
})

test_that("SerialParam and MulticoreParam give identical niches", {
  skip_on_os("windows")
  spe <- .toySPE()
  s1 <- buildNiches(spe, sigma = 20, BPPARAM = BiocParallel::SerialParam())
  s2 <- buildNiches(spe, sigma = 20, BPPARAM = BiocParallel::MulticoreParam(2))
  expect_equal(
    SingleCellExperiment::reducedDim(s1, "Niche20"),
    SingleCellExperiment::reducedDim(s2, "Niche20")
  )
})
