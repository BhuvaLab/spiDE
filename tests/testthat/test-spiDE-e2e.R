test_that("spiDE runs end-to-end and returns a populated results table", {
  spe <- .toySPE()
  res <- spiDE(spe,
    condition = "condition", sigma = c(10, 30, 50),
    covariates = "Age", verbose = FALSE
  )
  expect_s4_class(res, "SpiDEResults")
  tab <- results(res)
  expect_s3_class(tab, "data.frame")
  expect_true(nrow(tab) > 0)
  expect_true(all(c(
    "gene", "ct_index", "ct_niche", "bandwidth.max", "coef", "t",
    "DirectionGene", "DirectionIndex", "DirectionNiche",
    "fdr.gene", "fdr.index", "fdr.niche"
  ) %in% colnames(tab)))
  # all reported FDRs are valid probabilities
  expect_true(all(tab$fdr.niche >= 0 & tab$fdr.niche <= 1))
})

test_that("the planted neighbourhood effect is recovered and niche-specific", {
  spe <- .toySPE()
  res <- spiDE(spe,
    condition = "condition", sigma = c(10, 30, 50),
    covariates = "Age", verbose = FALSE
  )
  tab <- results(res)

  # the planted G1 (A index, B niche, up in Responders) is recovered
  ab <- tab[tab$gene == "G1" & tab$ct_index == "A" & tab$ct_niche == "B", ]
  expect_equal(nrow(ab), 1)
  expect_equal(ab$DirectionNiche, "Up")
  expect_gt(abs(ab$t), 5)

  # the signal is niche-specific: for G1 in index A, the B niche is the
  # strongest of all niche associations (the planted driver)
  g1a <- tab[tab$gene == "G1" & tab$ct_index == "A", ]
  expect_equal(g1a$ct_niche[which.max(abs(g1a$t))], "B")
})

test_that("combine = 'cauchy' produces valid results recovering G1/A/B", {
  spe <- .toySPE()
  res <- spiDE(spe,
    condition = "condition", sigma = c(10, 30, 50),
    covariates = "Age", combine = "cauchy", verbose = FALSE
  )
  expect_s4_class(res, "SpiDEResults")
  tab <- results(res)
  expect_true(nrow(tab) > 0)
  expect_true(all(tab$fdr.niche >= 0 & tab$fdr.niche <= 1))

  # the planted G1 (A index, B niche, up) is recovered and niche-specific
  ab <- tab[tab$gene == "G1" & tab$ct_index == "A" & tab$ct_niche == "B", ]
  expect_equal(nrow(ab), 1)
  expect_equal(ab$DirectionNiche, "Up")
  g1a <- tab[tab$gene == "G1" & tab$ct_index == "A", ]
  expect_equal(g1a$ct_niche[which.max(abs(g1a$t))], "B")
})

test_that("spiDE is deterministic and invariant to block.size", {
  # exact block-size invariance is a property of the reproducible CPU path;
  # the GPU path is float32 and its matmul is not bit-identical across
  # differing block (matrix) shapes, so pin the backend for this determinism
  # check (GPU block-size invariance is tested to tolerance in
  # test-gpuInference.R).
  spe <- .toySPE()
  a <- spiDE(spe, condition = "condition", sigma = c(10, 30),
    covariates = "Age", backend = "cpu", verbose = FALSE)
  b <- spiDE(spe, condition = "condition", sigma = c(10, 30),
    covariates = "Age", backend = "cpu", verbose = FALSE, block.size = 3)
  expect_equal(results(a), results(b))
})

test_that("testSpiDE errors if inference is missing and no spe is supplied", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  res <- fitSpiDE(spe, condition = "condition", sigma = 20, verbose = FALSE)
  expect_error(testSpiDE(res), "spe")
})
