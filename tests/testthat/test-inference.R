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

test_that(".wwFlat / .gramBatch reproduce the per-gene weighted Gram matrix", {
  set.seed(20)
  ncells <- 40
  p <- 5
  b <- 6
  W <- matrix(rnorm(ncells * p), ncells, p)
  wt <- matrix(runif(b * ncells, 0.5, 2), b, ncells)

  ref <- array(0, c(b, p, p))
  for (g in seq_len(b)) ref[g, , ] <- crossprod(W * wt[g, ], W)

  ww <- spiDE:::.wwFlat(W)
  expect_equal(dim(ww), c(ncells, p * p))
  info <- spiDE:::.gramBatch(ww, wt, p)
  expect_equal(dim(info), c(b, p, p))
  expect_equal(info, ref, tolerance = 1e-10)

  # ridge penalty is added to every slice's diagonal
  pen <- runif(p)
  ref_pen <- ref
  for (g in seq_len(b)) ref_pen[g, , ] <- ref[g, , ] + diag(pen)
  info_pen <- spiDE:::.gramBatch(ww, wt, p, penalty_diag = pen)
  expect_equal(info_pen, ref_pen, tolerance = 1e-10)
})

test_that(".subsetBatch and .batchDiag match base-array indexing", {
  set.seed(21)
  b <- 4
  p <- 5
  arr <- array(rnorm(b * p * p), c(b, p, p))
  sel <- c(2, 4)
  expect_equal(spiDE:::.subsetBatch(arr, sel), arr[, sel, sel, drop = FALSE])

  d <- spiDE:::.batchDiag(arr)
  ref <- t(vapply(seq_len(b), function(g) diag(arr[g, , ]), numeric(p)))
  expect_equal(d, ref)
})

test_that(".waldCauchyBlock returns a proper matrix for a single-gene block", {
  # regression test: vapply(uniq_index, FUN, numeric(b)) drops to a plain
  # (unnamed) vector rather than a (b, n_index) matrix when b == 1, which
  # silently corrupted cbind(Gene = p_gene, p_ct) for single-gene blocks.
  set.seed(22)
  ncells <- 30
  p <- 4
  W <- matrix(rnorm(ncells * p), ncells, p)
  colnames(W) <- paste0("c", seq_len(p))
  wtb <- matrix(runif(ncells, 0.5, 2), nrow = 1)
  alpha_block <- matrix(rnorm(p), nrow = 1)
  cov_niche <- c(FALSE, TRUE, TRUE, TRUE)
  index_ct <- c("A", "A", "B")
  uniq_index <- c("A", "B")

  ww <- spiDE:::.wwFlat(W)
  res <- spiDE:::.waldCauchyBlock(alpha_block, W, wtb, scale_block = 0.3,
                                  cov_niche = cov_niche, index_ct = index_ct,
                                  uniq_index = uniq_index, WW_flat = ww)
  expect_equal(dim(res$p.pos), c(1, 3))
  expect_equal(dim(res$p.neg), c(1, 3))
  expect_equal(colnames(res$p.pos), c("Gene", "A", "B"))
  expect_true(all(res$p.pos >= 0 & res$p.pos <= 1))
})

test_that("batched Cauchy path matches the per-gene .waldBrownGene loop (CPU)", {
  tf <- .toyFit()
  f <- spiDE:::.blockedInference(tf$fit, tf$Y, combine = "cauchy",
                                 backend = "cpu")

  # per-gene reference for the same fit
  W_full <- tf$fit@W
  covtype <- as.character(tf$fit@covtype)
  cols_gene <- grepl("Response", covtype)
  Wsub <- W_full[, cols_gene, drop = FALSE]
  cov_niche <- covtype[cols_gene] == "ResponseNiche"
  coefmap_sub <- tf$fit@coefmap[cols_gene, , drop = FALSE]
  index_ct <- coefmap_sub$index[cov_niche]
  uniq_index <- unique(index_ct)
  alpha_sub <- tf$fit@alpha[, cols_gene, drop = FALSE]
  psi <- tf$fit@psi
  mu <- SpaNorm::calculateMu(rep(0, nrow(tf$fit@alpha)), tf$fit@alpha, W_full)
  wt <- 1 / (1 / mu + psi)

  ref <- lapply(seq_len(nrow(alpha_sub)), function(i) {
    spiDE:::.waldBrownGene(alpha_sub[i, ], Wsub, wt[i, ], psi[i], cov_niche,
                           index_ct, uniq_index, combine = "cauchy")
  })
  ref_t <- do.call(rbind, lapply(ref, `[[`, "t_stat"))
  ref_ppos <- do.call(rbind, lapply(ref, `[[`, "p.pos"))

  expect_equal(unname(f@t_stat), unname(ref_t), tolerance = 1e-8)
  expect_equal(unname(f@p.combined.pos), unname(ref_ppos), tolerance = 1e-8)
})
