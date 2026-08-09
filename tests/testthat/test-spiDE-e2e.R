# The first three tests asserted on IDENTICAL pipeline runs: tests 1 and 2 used
# the default combiner and test 3 passed combine = "cauchy", which IS the
# default. Three full three-bandwidth runs of the same thing cost ~340 s of the
# 511 s this file used to take, against Bioconductor's 600 s budget for the
# whole check. The run is hoisted here and shared; every assertion is unchanged.
#
# The combiner default is asserted explicitly below, so if it ever changes this
# file starts testing something other than Cauchy loudly rather than silently.
spe_e2e <- .toySPE()
res_e2e <- spiDE(spe_e2e, condition = "condition", sigma = c(10, 30, 50),
                 random = "none",
                 covariates = "Age", combine = "cauchy", verbose = FALSE)
tab_e2e <- results(res_e2e)

test_that("cauchy is still the default combiner (guards the shared fit above)", {
  # formals() on an S4 method's @.Data returns a signature object rather than
  # the argument list, so the default has to be read off the printed method.
  tx <- paste(capture.output(print(selectMethod("testSpiDE", "SpiDEResults"))),
              collapse = " ")
  hit <- regmatches(tx, regexpr('combine = c\\([^)]*\\)', tx))
  expect_length(hit, 1L)
  expect_match(hit, '^combine = c\\(\\s*"cauchy"')
})

test_that("spiDE runs end-to-end and returns a populated results table", {
  res <- res_e2e
  expect_s4_class(res, "SpiDEResults")
  tab <- tab_e2e
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
  tab <- tab_e2e

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
  res <- res_e2e
  expect_s4_class(res, "SpiDEResults")
  tab <- tab_e2e
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
  # One bandwidth, not two: block.size partitions GENES inside the inference
  # stage, so its invariance is per-bandwidth by construction and a second
  # sigma doubles the fits without adding coverage. Cross-bandwidth combination
  # is exercised by the shared three-bandwidth run at the top of this file.
  spe <- .toySPE()
  a <- spiDE(spe, condition = "condition", sigma = 30, random = "none",
    covariates = "Age", backend = "cpu", verbose = FALSE)
  b <- spiDE(spe, condition = "condition", sigma = 30, random = "none",
    covariates = "Age", backend = "cpu", verbose = FALSE, block.size = 3)
  expect_equal(results(a), results(b))
})

test_that("testSpiDE errors if inference is missing and no spe is supplied", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  res <- fitSpiDE(spe, condition = "condition", sigma = 20, random = "none", verbose = FALSE)
  expect_error(testSpiDE(res), "spe")
})
