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

test_that("both df rules equal their direct limma computation, and only the df move", {
  s <- dfFamily()
  w <- 1 / (s$v + s$tau); w[!is.finite(s$b)] <- NA
  w <- w / rowMeans(w, na.rm = TRUE)
  f <- limma::lmFit(s$b, matrix(1, ncol(s$b), 1), weights = w)
  f$Amean <- log(s$mean_expr + 1e-3)
  f <- limma::eBayes(f, robust = TRUE, trend = TRUE)
  neff <- rowSums(w, na.rm = TRUE)^2 / rowSums(w^2, na.rm = TRUE)
  m <- rowSums(is.finite(w))
  pro <- spiDE:::.pooledColumnTest(s$b, s$v, s$tau, s$mean_expr)
  cap <- spiDE:::.pooledColumnTest(s$b, s$v, s$tau, s$mean_expr, df.rule = "capped")
  ok <- rowSums(is.finite(s$b)) >= 6
  expect_equal(pro$df[ok], (f$df.total * neff / m)[ok], tolerance = 1e-12)
  expect_equal(cap$df[ok], pmin(f$df.total, pmax(neff - 1, 1))[ok], tolerance = 1e-12)
  expect_identical(pro$t, cap$t)
  expect_identical(pro$estimate, cap$estimate)
  # limma's df is the upper bound of both; with unequal weights the proportional df exceeds the cap here
  expect_true(all(pro$df[ok] <= f$df.total[ok] + 1e-12))
  expect_gt(mean(pro$df[ok] > cap$df[ok]), 0.5)
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
})

test_that("the slopes engine's defaults are the heterogeneity family and the proportional df", {
  r <- testSpiDE(fit16h, condition = "condition")
  expect_identical(r@procedure, "heterogeneity")
  expect_identical(r@pooled.df, "proportional")
  expect_identical(r@table, testSpiDE(fit16h, condition = "condition", procedure = "heterogeneity",
                                      pooled.df = "proportional")@table)
  tb <- r@table
  pl <- tb[tb$test == "pooled" & is.finite(tb$p.heterogeneity), ]
  pass <- paste(pl$gene, pl$index, pl$niche)[stats::p.adjust(pl$p.heterogeneity, "BH") < r@fdr]
  cond <- tb[tb$test == "condition", ]
  key <- paste(cond$gene, cond$index, cond$niche)
  expect_identical(cond$in_family, is.finite(cond$p) & key %in% pass)
  expect_equal(cond$q[cond$in_family], stats::p.adjust(cond$p[cond$in_family], "BH"))
  expect_true(all(is.na(cond$q[!cond$in_family])))
  # each condition row carries its triplet's heterogeneity p
  pk <- paste(tb$gene, tb$index, tb$niche)[tb$test == "pooled"]
  expect_identical(cond$p.heterogeneity, tb$p.heterogeneity[tb$test == "pooled"][match(key, pk)])
})

test_that("the previous defaults stay available and only the family and the pooled df move", {
  new <- testSpiDE(fit16h, condition = "condition")@table
  old <- testSpiDE(fit16h, condition = "condition", procedure = "filtered", pooled.df = "capped")@table
  cl <- new$test == "condition"
  # the condition test itself does not depend on either choice
  expect_identical(new[cl, c("gene", "niche", "estimate", "se", "t", "df", "p")],
                   old[cl, c("gene", "niche", "estimate", "se", "t", "df", "p")])
  pl <- new$test == "pooled"
  expect_identical(new$t[pl], old$t[pl])
  # the previous filter is the pooled test at the fdr, under the capped df
  pk <- paste(old$gene, old$niche)[pl & old$q < 0.05]
  expect_true(all(paste(old$gene, old$niche)[cl & old$in_family] %in% pk))
})

test_that("the sandwich engine keeps the pooled filter and refuses the heterogeneity filter", {
  fs <- fitSpiDE(spe16h, condition = "condition", index = "A", sigma = 30, engine = "sandwich", verbose = FALSE)
  r <- testSpiDE(fs)
  expect_identical(r@procedure, "filtered")
  expect_identical(r@pooled.df, character())
  expect_true(all(is.na(r@table$p.heterogeneity)))
  expect_error(testSpiDE(fs, procedure = "heterogeneity"), "slopes engine")
  expect_identical(testSpiDE(fs, pooled.df = "capped")@table, r@table)
})

test_that("spiGSEA takes the pooled df of the results it is given; results saved before 0.99.37 were capped", {
  sets <- list(first = paste0("G", 1:5), second = paste0("G", 6:12))
  rp <- testSpiDE(fit16h, condition = "condition")
  rc <- testSpiDE(fit16h, condition = "condition", pooled.df = "capped")
  gp <- spiGSEA(rp, sets, test = "pooled", min.size = 3)
  gc <- spiGSEA(rc, sets, test = "pooled", min.size = 3)
  expect_identical(gp, spiGSEA(rc, sets, test = "pooled", min.size = 3, pooled.df = "proportional"))
  expect_identical(gc, spiGSEA(rp, sets, test = "pooled", min.size = 3, pooled.df = "capped"))
  expect_identical(gp$t, gc$t)
  old <- rp
  attr(old, "pooled.df") <- NULL
  expect_false(methods::.hasSlot(old, "pooled.df"))
  expect_identical(spiDE:::.pooledDfOf(old), "capped")
  expect_identical(spiGSEA(old, sets, test = "pooled", min.size = 3), gc)
  expect_output(show(rp), "proportional df")
})

test_that("the defaults are identical serially and in parallel", {
  skip_on_os("windows")
  r1 <- testSpiDE(fit16h, condition = "condition")@table
  r2 <- testSpiDE(fit16h, condition = "condition", BPPARAM = BiocParallel::MulticoreParam(2))@table
  expect_identical(r1, r2)
})
