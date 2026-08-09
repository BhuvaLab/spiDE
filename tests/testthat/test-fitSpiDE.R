test_that("fitSpiDE returns a SpiDEResults with one fit per bandwidth", {
  spe <- buildNiches(.toySPE(), sigma = c(10, 20))
  res <- fitSpiDE(spe, condition = "condition", random = "none", verbose = FALSE)

  expect_s4_class(res, "SpiDEResults")
  expect_equal(bandwidths(res), c(10, 20))
  expect_length(fits(res), 2)
  expect_true(all(vapply(fits(res), validObject, logical(1))))
})

test_that("the NB fit produces valid dispersions and log-likelihoods", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  res <- fitSpiDE(spe, condition = "condition", sigma = 20, random = "none", verbose = FALSE)
  f <- fits(res)[[1]]

  expect_equal(dim(f@alpha), c(f@ngenes, ncol(f@W)))
  expect_false(anyNA(f@psi))
  expect_true(all(f@psi > 0))
  expect_true(all(is.finite(f@loglik)))
  expect_true(all(f@loglik <= 0))
})

test_that("fitSpiDE recovers the planted B-niche effect on G1 in A cells", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  spe <- computeSizeFactors(spe, count = "nCount", area = "Area")
  res <- fitSpiDE(spe,
    condition = "condition", sigma = 20, random = "none",
    covariates = c("Age", "LS"), verbose = FALSE
  )
  f <- fits(res)[[1]]
  cm <- f@coefmap
  col <- cm$covariate[cm$type == "ResponseNiche" &
    cm$index == "A" & cm$niche == "B"]
  # G1 should have the largest coefficient at the A:B ResponseNiche column
  expect_equal(names(which.max(f@alpha[, col])), "G1")
})

test_that("fitSpiDE errors when niches are missing", {
  spe <- .toySPE()
  expect_error(fitSpiDE(spe, condition = "condition"), "buildNiches|niche")
})
