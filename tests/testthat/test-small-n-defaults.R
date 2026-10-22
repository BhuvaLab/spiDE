# The two defaults of 0.99.37, chosen for designs of 5-7 patients per condition: the
# pooled test's proportional df (research/smalldf/README.md) and the condition test's
# heterogeneity filter (research/condtest/README.md), checked against direct formulas.

spe16h <- buildNiches(.toySPE(n_samples = 16, n_per = 120, n_genes = 12, seed = 2), sigma = 30, verbose = FALSE)
fit16h <- fitSpiDE(spe16h, index = "A", sigma = 30, verbose = FALSE)

# a synthetic column family: G genes x S patients, unequal patient precision
dfFamily <- function(G = 150, S = 10, seed = 21) {
  set.seed(seed)
  mean_expr <- exp(seq(-3, 3, length.out = G))
  v <- matrix(stats::rgamma(G * S, 2, 2 / 0.05), G, S)
  b <- matrix(stats::rnorm(G * S, 0.05, sqrt(v + 0.01)), G, S)
  b[1, 1:2] <- NA
  list(b = b, v = v, tau = rep(0.01, G), mean_expr = mean_expr)
}

test_that("the pooled test equals its direct limma computation with the proportional df", {
  s <- dfFamily()
  w <- 1 / (s$v + s$tau); w[!is.finite(s$b)] <- NA
  w <- w / rowMeans(w, na.rm = TRUE)
  f <- limma::lmFit(s$b, matrix(1, ncol(s$b), 1), weights = w)
  f$Amean <- log(s$mean_expr + 1e-3)
  f <- limma::eBayes(f, robust = TRUE, trend = TRUE)
  neff <- rowSums(w, na.rm = TRUE)^2 / rowSums(w^2, na.rm = TRUE)
  m <- rowSums(is.finite(w))
  pro <- spiDE:::.pooledColumnTest(s$b, s$v, s$tau, s$mean_expr)
  ok <- rowSums(is.finite(s$b)) >= 6
  expect_equal(pro$df[ok], (f$df.total * neff / m)[ok], tolerance = 1e-12)
  # limma's df is the upper bound
  expect_true(all(pro$df[ok] <= f$df.total[ok] + 1e-12))
})

test_that("the heterogeneity p is Cochran's Q across patients, the statistic behind tau2", {
  x <- fit16h@index$A
  vp <- sweep(x$v_model, c(2, 3), x$factor, "*")
  hp <- spiDE:::.heterogeneityP(x$beta, vp)
  tau <- spiDE:::.dlTau2(x$beta, vp)
  expect_identical(dim(hp), c(length(x$genes), length(x$niches)))
  for (g in seq_along(x$genes)) for (j in seq_along(x$niches)) {
    b <- x$beta[g, , j]; w <- 1 / vp[g, , j]
    ok <- is.finite(b) & is.finite(w) & w > 0
    if (sum(ok) < 2) { expect_true(is.na(hp[g, j])); next }
    bbar <- sum(w[ok] * b[ok]) / sum(w[ok])
    Q <- sum(w[ok] * (b[ok] - bbar)^2)
    expect_equal(hp[g, j], stats::pchisq(Q, sum(ok) - 1, lower.tail = FALSE), tolerance = 1e-10)
    # tau2 is positive exactly when Q exceeds its null expectation m - 1
    expect_identical(tau[g, j] > 0, Q > sum(ok) - 1)
  }
  # one patient is not enough
  b1 <- x$beta; b1[1, -1, 1] <- NA
  expect_true(is.na(spiDE:::.heterogeneityP(b1, vp)[1, 1]))
  # two patients with variances near zero overflow Q (finite weights, infinite sum): NA, not p = 0
  bo <- array(c(0, 3, 1), c(1, 3, 1)); vo <- array(c(2e-308, 2e-308, 1), c(1, 3, 1))
  expect_true(is.na(spiDE:::.heterogeneityP(bo, vo)[1, 1]))
})

test_that("the slopes engine's default is the pooled-or-heterogeneity family", {
  r <- testSpiDE(fit16h, condition = "condition")
  expect_identical(r@procedure, "pooled_or_heterogeneity")
  expect_identical(r@table, testSpiDE(fit16h, condition = "condition", procedure = "pooled_or_heterogeneity")@table)
  tb <- r@table
  pl <- tb[tb$test == "pooled" & is.finite(tb$p.heterogeneity), ]
  het <- paste(pl$gene, pl$index, pl$niche)[stats::p.adjust(pl$p.heterogeneity, "BH") < r@fdr]
  pp <- tb[tb$test == "pooled" & is.finite(tb$p), ]
  pooled <- paste(pp$gene, pp$index, pp$niche)[pp$q < r@fdr]
  cond <- tb[tb$test == "condition", ]
  key <- paste(cond$gene, cond$index, cond$niche)
  fam <- function(proc) testSpiDE(fit16h, condition = "condition", procedure = proc)@table
  for (proc in c("pooled_or_heterogeneity", "heterogeneity", "filtered")) {
    pass <- switch(proc, pooled_or_heterogeneity = union(pooled, het), heterogeneity = het, filtered = pooled)
    cf <- fam(proc); cf <- cf[cf$test == "condition", ]
    expect_identical(cf$in_family, is.finite(cf$p) & key %in% pass)
    expect_equal(cf$q[cf$in_family], stats::p.adjust(cf$p[cf$in_family], "BH"))
    expect_true(all(is.na(cf$q[!cf$in_family])))
  }
  # the default family contains every triplet the previous default tested
  old <- fam("filtered"); new <- tb
  expect_true(all(new$in_family[new$test == "condition"][old$in_family[old$test == "condition"]]))
  # each condition row carries its triplet's heterogeneity p
  pk <- paste(tb$gene, tb$index, tb$niche)[tb$test == "pooled"]
  expect_identical(cond$p.heterogeneity, tb$p.heterogeneity[tb$test == "pooled"][match(key, pk)])
})

test_that("the previous family stays available and only the family moves", {
  new <- testSpiDE(fit16h, condition = "condition")@table
  old <- testSpiDE(fit16h, condition = "condition", procedure = "filtered")@table
  cl <- new$test == "condition"
  # the condition test itself does not depend on the family
  expect_identical(new[cl, c("gene", "niche", "estimate", "se", "t", "df", "p")],
                   old[cl, c("gene", "niche", "estimate", "se", "t", "df", "p")])
  pl <- new$test == "pooled"
  expect_identical(new$t[pl], old$t[pl])
  # the previous filter is the pooled test at the fdr
  pk <- paste(old$gene, old$niche)[pl & old$q < 0.05]
  expect_true(all(paste(old$gene, old$niche)[cl & old$in_family] %in% pk))
})

test_that("the sandwich engine keeps the pooled filter and refuses the heterogeneity filter", {
  fs <- fitSpiDE(spe16h, condition = "condition", index = "A", sigma = 30, engine = "sandwich", verbose = FALSE)
  r <- testSpiDE(fs)
  expect_identical(r@procedure, "filtered")
  expect_true(all(is.na(r@table$p.heterogeneity)))
  expect_error(testSpiDE(fs, procedure = "heterogeneity"), "slopes engine")
  expect_error(testSpiDE(fs, procedure = "pooled_or_heterogeneity"), "slopes engine")
  # spiDE() refuses the combination before it fits anything
  t0 <- proc.time()[["elapsed"]]
  expect_error(spiDE(spe16h, condition = "condition", sigma = 30, index = "A", engine = "sandwich",
                     procedure = "heterogeneity", verbose = FALSE), "slopes engine")
  expect_lt(proc.time()[["elapsed"]] - t0, 2)
})

test_that("show() summarises both tests of the slopes engine", {
  rp <- testSpiDE(fit16h, condition = "condition")
  expect_output(show(rp), "pooled across patients")
  expect_output(show(rp), "condition-specific")
})

test_that("the defaults are identical serially and in parallel", {
  skip_on_os("windows")
  r1 <- testSpiDE(fit16h, condition = "condition")@table
  r2 <- testSpiDE(fit16h, condition = "condition", BPPARAM = BiocParallel::MulticoreParam(2))@table
  expect_identical(r1, r2)
})

test_that(".cochran and .kishN are the direct formulas", {
  x <- fit16h@index$A
  vp <- sweep(x$v_model, c(2, 3), x$factor, "*")
  cq <- spiDE:::.cochran(x$beta, vp)
  g <- 2L; j <- 1L
  b <- x$beta[g, , j]; w <- 1 / vp[g, , j]
  ok <- is.finite(b) & is.finite(w) & w > 0
  expect_equal(cq$m[g, j], sum(ok))
  expect_equal(cq$sw[g, j], sum(w[ok]), tolerance = 1e-12)
  expect_equal(cq$sw2[g, j], sum(w[ok]^2), tolerance = 1e-12)
  expect_equal(cq$Q[g, j], sum(w[ok] * (b[ok] - sum(w[ok] * b[ok]) / sum(w[ok]))^2), tolerance = 1e-10)
  W <- rbind(c(1, 1, 1, 1), c(1, 2, NA, 3), c(NA, NA, NA, NA))
  expect_equal(spiDE:::.kishN(W)[1:2], c(4, 36 / 14))
})
