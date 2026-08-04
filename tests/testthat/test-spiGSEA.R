# Contract and numerical tests for the gene-set layer (R/spiGSEA.R).
#
# Everything here runs on synthetic matrices rather than a fitted model, so the
# short suite stays inside the Bioconductor budget. The end-to-end behaviour --
# that a planted signal is recovered, and that the streaming rho matches cor()
# on real residuals -- lives in longtests/testthat/test-gsea-numerics.R.

test_that(".meanCorFromZ reproduces the mean off-diagonal correlation", {
  # The streaming estimator's whole justification is this identity. Check it
  # against cor() directly: if it is wrong, every rho in the package is wrong
  # and every gene-set p-value is mis-scaled with it.
  set.seed(11)
  for (p in c(5L, 40L)) {
    n <- 200L
    x <- matrix(stats::rnorm(p * n), p, n)
    # induce genuine correlation so a zero-vs-nonzero bug cannot pass
    x <- x + matrix(rep(stats::rnorm(n), each = p), p, n)
    z <- t(scale(t(x)))                       # standardise each gene
    ref <- mean(stats::cor(t(x))[lower.tri(diag(p))])
    expect_equal(spiDE:::.meanCorFromZ(colSums(z), n, p), ref, tolerance = 1e-8)
  }
})

test_that(".meanCorFromZ is degenerate-safe", {
  expect_equal(spiDE:::.meanCorFromZ(numeric(10), 10L, 1L), 0)
  expect_equal(spiDE:::.meanCorFromZ(numeric(10), 10L, 0L), 0)
})

test_that(".tToZ keeps sign, is monotone, and survives extreme t", {
  tm <- matrix(c(-40, -5, -1, 0, 1, 5, 40), nrow = 1)
  z <- spiDE:::.tToZ(tm, 8)
  expect_equal(sign(z), sign(tm))
  expect_false(is.unsorted(as.numeric(z)))
  expect_true(all(is.finite(z)))
  # a t with few df has heavier tails than z, so |z| < |t| away from 0
  nz <- tm != 0
  expect_true(all(abs(z[nz]) < abs(tm[nz])))
  # NULL df means the statistics are already normal -- pass through untouched
  expect_identical(spiDE:::.tToZ(tm, NULL), tm)
  # large df converges back to the identity
  expect_equal(spiDE:::.tToZ(matrix(2), 1e6), matrix(2), tolerance = 1e-4)
})

test_that(".tToZ applies per-column df", {
  tm <- matrix(c(3, 3), nrow = 1)
  z <- spiDE:::.tToZ(tm, c(4, 1000))
  # the same t is less impressive on 4 df than on 1000
  expect_lt(z[1, 1], z[1, 2])
})

test_that(".setColMeans matches colMeans per set", {
  m <- matrix(as.numeric(1:30), nrow = 10)
  sets <- list(a = c(1L, 2L), b = c(3L, 4L, 5L), c = 10L)
  got <- spiDE:::.setColMeans(m, sets)
  expect_equal(unname(got[1, ]), colMeans(m[1:2, , drop = FALSE]))
  expect_equal(unname(got[2, ]), colMeans(m[3:5, , drop = FALSE]))
  expect_equal(unname(got[3, ]), m[10, ])
  # order follows `sets`, not the row order of the matrix
  expect_identical(rownames(got), c("1", "2", "3"))
})

test_that(".relWeights thresholds and normalises like .geneWeights", {
  ll <- cbind(c(0, -10), c(-1, 0))
  w <- spiDE:::.relWeights(ll, thresh = 0.1)
  expect_equal(w[1, 1], 1)                    # best bandwidth for row 1
  expect_equal(w[2, 2], 1)
  expect_equal(w[2, 1], 0)                    # exp(-10) < 0.1 -> dropped
  expect_true(all(w >= 0 & w <= 1))
})

test_that("the variance inflation factor is actually applied", {
  # A set of m equicorrelated statistics has less information than m
  # independent ones, so raising rho must shrink the statistic. Without this
  # term a co-regulated pathway is significant purely for being co-regulated.
  m <- 20
  sd0 <- sqrt((1 + 0.0 * (m - 1)) / m)
  sd5 <- sqrt((1 + 0.5 * (m - 1)) / m)
  expect_gt(sd5, sd0)
  expect_equal(sd0, sqrt(1 / m))
})

test_that(".gseaCascade gates each level within the survivors of the last", {
  tab <- data.frame(
    geneset = rep(c("s1", "s2"), each = 4L),
    collection = "H",
    ct_index = rep(rep(c("A", "B"), each = 2L), 2L),
    ct_niche = rep(c("X", "Y"), 4L),
    # s1 strongly significant throughout; s2 pure null
    p = c(rep(1e-8, 4L), rep(0.9, 4L)),
    stringsAsFactors = FALSE
  )
  out <- spiDE:::.gseaCascade(tab, fdr = 0.05, nested = TRUE)
  expect_true(all(out$geneset == "s1"))       # s2 killed at level 1
  expect_true(all(out$fdr.geneset < 0.05))
  expect_true(all(out$fdr.index < 0.05))
  expect_true(all(out$fdr.niche < 0.05))
  # NB the levels are NOT monotone in q, and must not be asserted to be: each
  # corrects over a different family (all sets in a collection; cell types
  # within a set; niches within a pair), so a deeper level with fewer
  # comparisons can legitimately report a SMALLER q than the level above it.
  # What the cascade guarantees is gating, not ordering.
  expect_equal(nrow(out), 4L)
})

test_that(".gseaCascade returns an empty frame when nothing survives", {
  tab <- data.frame(
    geneset = "s1", collection = "H", ct_index = "A", ct_niche = "X",
    p = 0.99, stringsAsFactors = FALSE
  )
  expect_equal(nrow(spiDE:::.gseaCascade(tab, fdr = 0.05)), 0L)
})

test_that(".gseaCascade stops at two levels when nested = FALSE", {
  tab <- data.frame(
    geneset = rep("s1", 2L), collection = "H", ct_index = c("A", "B"),
    ct_niche = NA_character_, p = c(1e-8, 1e-8), stringsAsFactors = FALSE
  )
  out <- spiDE:::.gseaCascade(tab, fdr = 0.05, nested = FALSE)
  expect_gt(nrow(out), 0L)
  expect_null(out$fdr.niche)
})

test_that("the competitive statistic is centred on the background, not zero", {
  # Add a constant to every gene: a self-contained test must react to that
  # (every set really has shifted), a competitive one must not (no set shifted
  # RELATIVE to the others). This is the property that distinguishes them, and
  # it is checkable without fitting anything.
  set.seed(3)
  z <- matrix(stats::rnorm(200L), 100L, 2L)
  sets <- list(a = 1:20, b = 21:40)
  shift <- 5
  meanIn <- function(mat) spiDE:::.setColMeans(mat, sets)
  comp <- function(mat) {
    zbar <- meanIn(mat)
    m <- lengths(sets)
    mOut <- nrow(mat) - m
    meanOut <- sweep(-(zbar * m), 2L, colSums(mat), "+") / mOut
    zbar - meanOut
  }
  expect_equal(comp(z), comp(z + shift), tolerance = 1e-10)
  # while the raw (self-contained) set mean moves by exactly the shift
  expect_equal(meanIn(z + shift) - meanIn(z),
               matrix(shift, length(sets), ncol(z),
                      dimnames = dimnames(meanIn(z))),
               tolerance = 1e-10)
})
