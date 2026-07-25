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

test_that(".gramBatch reproduces the per-gene weighted Gram matrix", {
  set.seed(20)
  ncells <- 40
  p <- 5
  b <- 6
  W <- matrix(rnorm(ncells * p), ncells, p)
  wt <- matrix(runif(b * ncells, 0.5, 2), b, ncells)

  ref <- array(0, c(b, p, p))
  for (g in seq_len(b)) ref[g, , ] <- crossprod(W * wt[g, ], W)

  info <- spiDE:::.gramBatch(W, wt)
  expect_equal(dim(info), c(b, p, p))
  expect_equal(info, ref, tolerance = 1e-10)

  # ridge penalty is added to every slice's diagonal
  pen <- runif(p)
  ref_pen <- ref
  for (g in seq_len(b)) ref_pen[g, , ] <- ref[g, , ] + diag(pen)
  info_pen <- spiDE:::.gramBatch(W, wt, penalty_diag = pen)
  expect_equal(info_pen, ref_pen, tolerance = 1e-10)

  # sub-batching the genes must not change the result: the same slices come
  # back whether built all at once or a few at a time. This is the property
  # that makes .covBatchSize()'s bound safe to apply at any size.
  chunked <- array(0, c(b, p, p))
  for (ii in list(1:2, 3:4, 5:6)) {
    chunked[ii, , ] <- spiDE:::.gramBatch(W, wt[ii, , drop = FALSE])
  }
  expect_equal(chunked, ref, tolerance = 1e-10)
})

test_that(".covBatchSize shrinks with design width and stays >= 1", {
  # the point of the sub-batch: covariance memory is linear in p, so a wider
  # design yields a smaller batch rather than an unbounded allocation
  narrow <- spiDE:::.covBatchSize(20000, 50, "cpu")
  wide <- spiDE:::.covBatchSize(20000, 5000, "cpu")
  expect_gt(narrow, wide)
  expect_gte(wide, 1L)

  # a design too wide for any budget still yields a usable batch of 1,
  # never 0 or negative
  expect_equal(spiDE:::.covBatchSize(20000, 5e5, "cpu"), 1L)

  # a larger budget allows a larger batch
  big <- withr::with_options(
    list(spiDE.cov.mem.budget = 2e10),
    spiDE:::.covBatchSize(20000, 5000, "cpu"))
  expect_gt(big, wide)
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

  res <- spiDE:::.waldCauchyBlock(alpha_block, W, wtb, scale_block = 0.3,
                                  cov_niche = cov_niche, index_ct = index_ct,
                                  uniq_index = uniq_index, Wgram = W)
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

test_that(".waldCauchyBlock is invariant to the covariance sub-batch size", {
  # Regression test for the Khatri-Rao blowup: the batched Cauchy path used
  # to precompute an (ncells x p^2) cross term once per bandwidth, outside
  # the gene loop, so its memory could not be bounded by any block size --
  # 63 GB for a 602-column random-intercept design, 4.2 TB for a
  # 4906-column random-slope one, making combine = "cauchy" unusable with
  # random effects on both backends. The fix bounds it via a gene sub-batch
  # instead, which is only sound if results do not depend on that size.
  set.seed(23)
  ncells <- 60
  p <- 8
  b <- 7
  W <- matrix(rnorm(ncells * p), ncells, p)
  colnames(W) <- paste0("c", seq_len(p))
  wtb <- matrix(runif(b * ncells, 0.5, 2), b, ncells)
  alpha_block <- matrix(rnorm(b * p), b, p)
  scale_block <- runif(b, 0.2, 0.5)
  cov_niche <- c(FALSE, rep(TRUE, p - 1))
  index_ct <- rep(c("A", "B"), length.out = p - 1)
  uniq_index <- c("A", "B")

  run <- function(cb) {
    spiDE:::.waldCauchyBlock(alpha_block, W, wtb, scale_block,
                             cov_niche = cov_niche, index_ct = index_ct,
                             uniq_index = uniq_index, Wgram = W,
                             cov.batch = cb)
  }
  full <- run(NULL)
  for (cb in c(1, 2, 3, b, b + 5)) {
    expect_equal(run(cb), full, tolerance = 1e-10,
                 info = paste("cov.batch =", cb))
  }
})

test_that("mixed-effects Cauchy inference runs on a design too wide to flatten", {
  # p = 120 makes the old (ncells x p^2) cross term 60x larger than the
  # design itself; the sub-batched path handles it in bounded memory and
  # must still agree with the per-gene reference exactly.
  set.seed(24)
  ncells <- 200
  p <- 120
  b <- 4
  W <- matrix(rnorm(ncells * p), ncells, p)
  colnames(W) <- paste0("c", seq_len(p))
  sel <- 1:6
  wtb <- matrix(runif(b * ncells, 0.5, 2), b, ncells)
  alpha_block <- matrix(rnorm(b * length(sel)), b, length(sel))
  scale_block <- runif(b, 0.2, 0.5)
  penalty <- c(rep(0, 20), rep(2, p - 20))
  cov_niche <- c(FALSE, rep(TRUE, length(sel) - 1))
  index_ct <- rep(c("A", "B"), length.out = length(sel) - 1)
  uniq_index <- c("A", "B")
  Wsub <- W[, sel, drop = FALSE]

  res <- spiDE:::.waldCauchyBlock(alpha_block, Wsub, wtb, scale_block,
                                  cov_niche = cov_niche, index_ct = index_ct,
                                  uniq_index = uniq_index, Wgram = W,
                                  W_full = W, penalty = penalty, sel = sel,
                                  df = 10, cov.batch = 2)

  # per-gene reference: the full penalised inverse, restricted to sel
  ref_se <- t(vapply(seq_len(b), function(g) {
    vc <- solve(crossprod(W * wtb[g, ], W) + diag(penalty))
    sqrt(scale_block[g] * diag(vc)[sel])
  }, numeric(length(sel))))
  expect_equal(unname(res$se), unname(ref_se), tolerance = 1e-8)
  expect_true(all(is.finite(res$t_stat)))
})
