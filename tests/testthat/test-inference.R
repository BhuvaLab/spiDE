# helper: a fitted SpiDEFit on the toy data
.toyFit <- function(sigma = 20) {
  spe <- buildNiches(.toySPE(), sigma = sigma)
  spe <- computeSizeFactors(spe, count = "nCount", area = "Area")
  res <- fitSpiDE(spe,
    condition = "condition", sigma = sigma,
    covariates = c("Age", "LS"), verbose = FALSE
  )
  list(fit = fits(res)[[1]], Y = as.matrix(SummarizedExperiment::assay(spe, "counts")))
}

test_that(".waldBrownGene: se and t-stat match the closed form", {
  set.seed(1)
  ncells <- 200
  W <- cbind(1, scale(runif(ncells)), scale(runif(ncells)))
  colnames(W) <- c("Response", "CellTypeA:R:B", "CellTypeC:R:B")
  alpha_g <- c(0.2, 1.0, -0.5)
  wt <- runif(ncells, 0.5, 1.5)
  psi <- 0.3
  cov_niche <- c(FALSE, TRUE, TRUE)
  res <- spiDE:::.waldBrownGene(alpha_g, W, wt, psi, cov_niche,
    index_ct = c("A", "C"), uniq_index = c("A", "C"))

  varcov <- solve(crossprod(W * wt, W))
  expect_equal(unname(res$se), unname(sqrt(psi * diag(varcov))), tolerance = 1e-8)
  expect_equal(unname(res$t_stat),
    unname(alpha_g / sqrt(psi * diag(varcov))), tolerance = 1e-8)
})

test_that(".waldBrownGene: p.pos and p.neg are complementary in direction", {
  set.seed(2)
  ncells <- 300
  W <- cbind(1, scale(runif(ncells)), scale(runif(ncells)))
  colnames(W) <- c("Response", "i1:B", "i2:B")
  # strongly positive niche coefficients -> small p.pos, large p.neg
  res <- spiDE:::.waldBrownGene(c(0, 3, 3), W, rep(1, ncells), 0.2,
    c(FALSE, TRUE, TRUE), c("i1", "i2"), c("i1", "i2"))
  expect_lt(res$p.pos["Gene"], 0.05)
  expect_gt(res$p.neg["Gene"], 0.5)
})

test_that("inference recovers the planted neighbourhood effect", {
  tf <- .toyFit()
  f <- spiDE:::.blockedInference(tf$fit, tf$Y)

  expect_true(all(f@p.combined.pos >= 0 & f@p.combined.pos <= 1))
  expect_true(all(f@p.combined.neg >= 0 & f@p.combined.neg <= 1))
  # G1 is the most significant gene for index A in the up direction
  expect_equal(names(which.min(f@p.combined.pos[, "A"])), "G1")
  # and its A-niche effect is up, not down
  expect_lt(f@p.combined.pos["G1", "A"], f@p.combined.neg["G1", "A"])
})

test_that("Cauchy combination recovers the planted effect and is valid", {
  tf <- .toyFit()
  f <- spiDE:::.blockedInference(tf$fit, tf$Y, combine = "cauchy")

  expect_true(all(f@p.combined.pos >= 0 & f@p.combined.pos <= 1))
  expect_true(all(f@p.combined.neg >= 0 & f@p.combined.neg <= 1))
  # G1 is still the most significant gene for index A in the up direction
  expect_equal(names(which.min(f@p.combined.pos[, "A"])), "G1")
  expect_lt(f@p.combined.pos["G1", "A"], f@p.combined.neg["G1", "A"])

  # Cauchy path is invariant to block.size
  fb <- spiDE:::.blockedInference(tf$fit, tf$Y, combine = "cauchy",
                                  block.size = 3)
  expect_equal(f@p.combined.pos, fb@p.combined.pos)
})

test_that("inference is invariant to block.size, BPPARAM and DelayedArray", {
  tf <- .toyFit()
  base <- spiDE:::.blockedInference(tf$fit, tf$Y)
  blocked <- spiDE:::.blockedInference(tf$fit, tf$Y, block.size = 3)
  expect_equal(base@t_stat, blocked@t_stat)
  expect_equal(base@p.combined.pos, blocked@p.combined.pos)
  expect_equal(base@loglik, blocked@loglik)

  skip_if_not_installed("DelayedArray")
  yd <- DelayedArray::DelayedArray(tf$Y)
  da <- spiDE:::.blockedInference(tf$fit, yd, block.size = 5)
  expect_equal(base@t_stat, da@t_stat)
  expect_equal(base@p.combined.pos, da@p.combined.pos)

  skip_on_os("windows")
  mc <- spiDE:::.blockedInference(tf$fit, tf$Y,
    block.size = 4, BPPARAM = BiocParallel::MulticoreParam(2))
  expect_equal(base@p.combined.pos, mc@p.combined.pos)
})

test_that(".chunkGenes partitions all genes exactly once", {
  expect_equal(spiDE:::.chunkGenes(10, NULL), list(seq_len(10)))
  ch <- spiDE:::.chunkGenes(10, 3)
  expect_equal(sort(unname(unlist(ch))), seq_len(10))
  expect_true(all(lengths(ch) <= 3))
})
