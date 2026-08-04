# Live-fit numerical tests for the gene-set layer. These need a real fit, so
# they run on the weekly long-tests builder rather than in the short suite.

test_that(".interGeneCor reproduces cor() on real Pearson residuals", {
  # The streaming estimator never forms the gene x gene correlation matrix. On
  # a fixture small enough to form it, the two must agree exactly -- if they do
  # not, every gene-set p-value in the package is mis-scaled.
  set.seed(1)
  data(toySpiDE)
  spe <- buildNiches(toySpiDE, sigma = 20)
  res <- spiDE(spe, condition = "condition", sigma = 20, verbose = FALSE)
  f <- fits(res)[[1]]
  Y <- SummarizedExperiment::assay(spe, "counts")

  mu <- as.matrix(SpaNorm::toRMatrix(
    SpaNorm::calculateMu(f@gmean, f@alpha, f@W)))
  r <- (as.matrix(Y) - mu) / sqrt(mu + f@psi * mu^2)
  ref <- mean(stats::cor(t(r))[lower.tri(diag(nrow(r)))])

  got <- spiDE:::.interGeneCor(Y, f, seq_len(nrow(spe)), block.size = 1000L,
                               backend = "cpu")
  expect_equal(got, ref, tolerance = 1e-8)
})

test_that(".interGeneCor is invariant to block size and to worker count", {
  # Blocks contribute additive partial sums, so neither chunking nor
  # parallelism may move the answer. A bug here would be silent: the estimate
  # would simply be a bit wrong, in a way that depends on the machine.
  set.seed(1)
  data(toySpiDE)
  spe <- buildNiches(toySpiDE, sigma = 20)
  res <- spiDE(spe, condition = "condition", sigma = 20, verbose = FALSE)
  f <- fits(res)[[1]]
  Y <- SummarizedExperiment::assay(spe, "counts")
  g <- seq_len(nrow(spe))

  one <- spiDE:::.interGeneCor(Y, f, g, block.size = 1000L, backend = "cpu")
  many <- spiDE:::.interGeneCor(Y, f, g, block.size = 3L, backend = "cpu")
  expect_equal(one, many, tolerance = 1e-10)

  par <- spiDE:::.interGeneCor(
    Y, f, g, block.size = 3L, backend = "cpu",
    BPPARAM = BiocParallel::MulticoreParam(2))
  expect_equal(one, par, tolerance = 1e-10)
})

test_that("the GPU backend agrees with the CPU backend", {
  skip_if_not(SpaNorm::checkGPU(), "torch GPU/MPS not available")
  set.seed(1)
  data(toySpiDE)
  spe <- buildNiches(toySpiDE, sigma = 20)
  res <- spiDE(spe, condition = "condition", sigma = 20, verbose = FALSE)
  f <- fits(res)[[1]]
  Y <- SummarizedExperiment::assay(spe, "counts")
  g <- seq_len(nrow(spe))
  cpu <- spiDE:::.interGeneCor(Y, f, g, block.size = 5L, backend = "cpu")
  gpu <- spiDE:::.interGeneCor(Y, f, g, block.size = 5L, backend = "gpu")
  # MPS is float32-only, so the tolerance defers to the device
  expect_equal(gpu, cpu, tolerance = gpu_tol())
})

test_that("spiGSEA recovers a planted gene-set signal", {
  set.seed(1)
  data(toySpiDE)
  spe <- buildNiches(toySpiDE, sigma = c(20, 40))
  res <- spiDE(spe, condition = "condition", sigma = c(20, 40), verbose = FALSE)
  gn <- rownames(spe)
  # G1 carries the planted Responder x B-niche effect in index cell type A
  sets <- list(planted = gn[1:6], other = gn[7:14])
  out <- spiGSEA(res, spe, sets, min.size = 3, fdr = 1, verbose = FALSE)

  expect_s3_class(out, "data.frame")
  expect_true(all(c("geneset", "ct_index", "ct_niche", "z", "Direction",
                    "fdr.geneset", "fdr.index", "fdr.niche") %in% names(out)))
  expect_true(all(out$size >= 3))
  # the planted set must be the strongest row for the B niche in index A
  ba <- out[out$ct_niche == "B" & out$ct_index == "A", , drop = FALSE]
  if (nrow(ba)) {
    expect_equal(ba$geneset[which.max(abs(ba$z))], "planted")
  }
})

test_that("spiGSEA errors informatively on an un-tested object", {
  set.seed(1)
  data(toySpiDE)
  spe <- buildNiches(toySpiDE, sigma = 20)
  f <- fitSpiDE(spe, "condition", sigma = 20, verbose = FALSE)
  expect_error(spiGSEA(f, spe, list(a = rownames(spe)[1:6]), verbose = FALSE),
               "testSpiDE")
})

test_that("spiGSEA errors when no set survives the size window", {
  set.seed(1)
  data(toySpiDE)
  spe <- buildNiches(toySpiDE, sigma = 20)
  res <- spiDE(spe, condition = "condition", sigma = 20, verbose = FALSE)
  expect_error(
    spiGSEA(res, spe, list(a = rownames(spe)[1:3]), min.size = 100,
            verbose = FALSE),
    "no gene set"
  )
})

test_that("spiGSEA handles a SINGLE gene set", {
  # Regression: with one set, vapply over bandwidths returns a length-k vector
  # instead of a 1 x k matrix, which .cauchyCombine then transposed against a
  # 1 x k weight matrix and errored on. Testing one set is ordinary usage.
  set.seed(1)
  data(toySpiDE)
  spe <- buildNiches(toySpiDE, sigma = c(20, 40))
  res <- spiDE(spe, condition = "condition", sigma = c(20, 40), verbose = FALSE)
  one <- list(only = rownames(spe)[1:6])
  expect_no_error(
    out <- spiGSEA(res, spe, one, min.size = 3, fdr = 1, verbose = FALSE))
  expect_true(all(out$geneset == "only"))
  # and it must agree with the same set tested alongside others
  many <- spiGSEA(res, spe, c(one, list(other = rownames(spe)[7:14])),
                  min.size = 3, fdr = 1, verbose = FALSE)
  a <- out[order(out$ct_index, out$ct_niche), ]
  b <- many[many$geneset == "only", ]
  b <- b[order(b$ct_index, b$ct_niche), ]
  expect_equal(a$z, b$z, tolerance = 1e-10)
})

test_that("testSpiDE() stores the correlation spiGSEA then applies", {
  # The slot is reached the same way as every other: fits(res)[[i]]$rho.
  set.seed(1)
  data(toySpiDE)
  spe <- buildNiches(toySpiDE, sigma = c(20, 40))
  res <- spiDE(spe, condition = "condition", sigma = c(20, 40), verbose = FALSE)
  rho <- vapply(fits(res), function(f) f$rho, numeric(1))
  expect_length(rho, 2L)
  expect_true(all(rho > -1 & rho < 1))
  out <- spiGSEA(res, genesets = list(a = rownames(spe)[1:6]), min.size = 3,
                 fdr = 1, verbose = FALSE)
  expect_equal(unname(attr(out, "rho")), unname(rho), tolerance = 0)
})
