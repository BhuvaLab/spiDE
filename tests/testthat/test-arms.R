# The test-stage arms of feature/shared-tau2 (research/sharedtau/README.md):
# the heterogeneity floor, the uncapped pooled df and the moderated condition
# variance, checked against direct formulas; the default arm is the shipped test.

spe16a <- buildNiches(.toySPE(n_samples = 16, n_per = 120, n_genes = 12, seed = 2), sigma = 30, verbose = FALSE)
fit16a <- fitSpiDE(spe16a, index = "A", sigma = 30, verbose = FALSE)

# a synthetic column family: G genes x S patients, heterogeneity rising with expression
synthFamily <- function(G = 150, S = 20, seed = 11) {
  set.seed(seed)
  mean_expr <- exp(seq(-3, 3, length.out = G))
  tau <- 0.02 * (1 + seq_len(G) / G)
  v <- matrix(stats::rgamma(G * S, 4, 4 / 0.05), G, S)
  b <- matrix(stats::rnorm(G * S, 0.1, sqrt(v + tau)), G, S)
  list(b = b, v = v, mean_expr = mean_expr)
}

test_that("the default arm is the shipped test", {
  expect_identical(spiDE:::.armSpec(), list(heterogeneity = "own", pooled_df = "capped", condition_variance = "hc2"))
  r0 <- testSpiDE(fit16a, condition = "condition")
  r1 <- testSpiDE(fit16a, condition = "condition", .arm = spiDE:::.armSpec())
  expect_identical(r0@table, r1@table)
})

test_that("the DL moments give the truncated DL tau2", {
  s <- synthFamily()
  b <- array(s$b, c(nrow(s$b), ncol(s$b), 1)); v <- array(s$v, dim(b))
  mo <- spiDE:::.dlMoments(b, v)
  w <- 1 / s$v
  bb <- rowSums(w * s$b) / rowSums(w)
  Q <- rowSums(w * (s$b - bb)^2)
  cc <- rowSums(w) - rowSums(w^2) / rowSums(w)
  expect_equal(as.numeric(mo$t), (Q - (ncol(s$b) - 1)) / cc, tolerance = 1e-12)
  expect_equal(as.numeric(spiDE:::.dlTau2(b, v)), pmax((Q - (ncol(s$b) - 1)) / cc, 0), tolerance = 1e-12)
})

test_that("the heterogeneity floor: pooled ratio under min.rows, loess reproduces a linear trend, never below own", {
  s <- synthFamily()
  b <- array(s$b, c(nrow(s$b), ncol(s$b), 1)); v <- array(s$v, dim(b))
  mo <- spiDE:::.dlMoments(b, v)
  own <- spiDE:::.dlTau2(b, v)
  # fewer usable rows than min.rows: the family value is sum(Q - (m - 1)) / sum(c)
  fl <- spiDE:::.heterogeneityFloor(b, v, s$mean_expr, min.rows = 1e6)
  pooled <- max(sum(mo$Qm) / sum(mo$c), 0)
  expect_equal(as.numeric(fl), pmax(as.numeric(own), pooled), tolerance = 1e-10)
  # the floor is never below the gene's own tau2, and finite everywhere
  fl2 <- spiDE:::.heterogeneityFloor(b, v, s$mean_expr)
  expect_true(all(is.finite(fl2)) && all(fl2 >= own - 1e-15))
  # identical genes give a constant floor
  bc <- array(rep(s$b[1, ], each = nrow(s$b)), dim(b)); vc <- array(rep(s$v[1, ], each = nrow(s$v)), dim(b))
  flc <- spiDE:::.heterogeneityFloor(bc, vc, s$mean_expr)
  expect_equal(diff(range(flc)), 0, tolerance = 1e-10)
})

test_that("the family trend reproduces a linear trend, clamps outside the range, pools small families", {
  x <- seq(-3, 3, length.out = 200); cw <- stats::runif(200, 0.5, 5)
  t <- 0.01 + 0.002 * x
  expect_equal(spiDE:::.familyTrend(t, cw, x, seq_along(x)), t, tolerance = 1e-10)
  xo <- c(x, 10, NA); to <- c(t, NA, NA); co <- c(cw, NA, NA)
  fo <- spiDE:::.familyTrend(to, co, xo, seq_along(x))
  expect_equal(fo[201], t[200], tolerance = 1e-10)
  expect_equal(fo[202], sum(cw * t) / sum(cw), tolerance = 1e-12)
  expect_equal(spiDE:::.familyTrend(t, cw, x, 1:10), rep(sum(cw[1:10] * t[1:10]) / sum(cw[1:10]), 200), tolerance = 1e-12)
})

test_that("the heterogeneity floor is finite outside the expression range and for genes without patients", {
  s <- synthFamily()
  b <- array(s$b, c(nrow(s$b), ncol(s$b), 1)); v <- array(s$v, dim(b))
  b[1, , 1] <- NA                                   # a gene with no usable patient
  me <- s$mean_expr; me[2] <- 1e6; me[3] <- NA      # outside the fitted range; unknown expression
  fl <- spiDE:::.heterogeneityFloor(b, v, me)
  expect_true(all(is.finite(fl)))
  expect_true(fl[1, 1] >= 0)
})

test_that("the heterogeneity floor and every arm's pooled test and filter never look at the condition", {
  pt <- fit16a@patients
  set.seed(4)
  fit_perm <- fit16a
  fit_perm@patients$perm <- factor(stats::setNames(sample(as.character(pt$condition)), pt$patient)[pt$patient])
  arms <- list(spiDE:::.armSpec("floor"), spiDE:::.armSpec("floor", "uncapped", "moderated"),
               spiDE:::.armSpec("equal", "uncapped"), spiDE:::.armSpec(condition_variance = "moderated"))
  for (a in arms) {
    r1 <- testSpiDE(fit16a, condition = "condition", .arm = a)@table
    r2 <- testSpiDE(fit_perm, condition = "perm", .arm = a)@table
    expect_identical(r1[r1$test == "pooled", c("gene", "niche", "p", "q")], r2[r2$test == "pooled", c("gene", "niche", "p", "q")])
  }
})

test_that("the uncapped pooled test refers limma's t to df.total", {
  s <- synthFamily()
  tau <- rep(0.02, nrow(s$b))
  capped <- spiDE:::.pooledColumnTest(s$b, s$v, tau, s$mean_expr)
  unc <- spiDE:::.pooledColumnTest(s$b, s$v, tau, s$mean_expr, df_rule = "uncapped")
  expect_identical(capped$t, unc$t)
  w <- 1 / (s$v + tau); w <- w / rowMeans(w)
  fit <- limma::eBayes(within(limma::lmFit(s$b, matrix(1, ncol(s$b), 1), weights = w),
                              Amean <- log(s$mean_expr + 1e-3)), robust = TRUE, trend = TRUE)
  expect_equal(unc$df, fit$df.total, tolerance = 1e-12)
  expect_true(all(capped$df <= unc$df + 1e-12))
})

test_that("the condition test's pieces rebuild the HC2 test, and moderation is a direct squeezeVar", {
  s <- synthFamily()
  trt <- rep(0:1, length.out = ncol(s$b))
  tau <- rep(0.02, nrow(s$b))
  M <- spiDE:::.conditionPieces(s$b, s$v, tau, trt)
  hc2 <- spiDE:::.robustConditionTest(s$b, s$v, tau, trt)
  expect_equal(hc2$se, sqrt(M[, "V"]))
  expect_equal(hc2$df, unname(M[, "nu"]))
  # A22 is the model variance of the condition coefficient under the unnormalised weights
  g <- 5
  X <- cbind(1, trt); W <- 1 / (s$v[g, ] + tau[g])
  expect_equal(unname(M[g, "A22"]), solve(crossprod(X * sqrt(W)))[2, 2], tolerance = 1e-10)
  mod <- spiDE:::.robustConditionTest(s$b, s$v, tau, trt, moderate = TRUE, mean_expr = s$mean_expr)
  ok <- which(is.finite(M[, "V"]))
  sq <- limma::squeezeVar(M[ok, "V"] / M[ok, "A22"], M[ok, "nu"], covariate = log(s$mean_expr[ok] + 1e-3), robust = TRUE)
  expect_equal(mod$se[ok], sqrt(sq$var.post * M[ok, "A22"]), tolerance = 1e-12)
  expect_equal(mod$df[ok], M[ok, "nu"] + rep_len(sq$df.prior, length(ok)), tolerance = 1e-12)
  expect_identical(mod$estimate, hc2$estimate)
})

test_that("moderation leaves untestable rows NA and small families unmoderated", {
  s <- synthFamily()
  trt <- rep(0:1, length.out = ncol(s$b))
  tau <- rep(0.02, nrow(s$b))
  b <- s$b; b[1:3, ] <- NA
  mod <- spiDE:::.robustConditionTest(b, s$v, tau, trt, moderate = TRUE, mean_expr = s$mean_expr)
  expect_true(all(is.na(mod$p[1:3])) && all(is.finite(mod$p[-(1:3)])))
  small <- spiDE:::.robustConditionTest(s$b[1:10, ], s$v[1:10, ], tau[1:10], trt, moderate = TRUE,
                                        mean_expr = s$mean_expr[1:10])
  expect_false(attr(small, "moderated"))
  attr(small, "moderated") <- NULL
  expect_identical(small, spiDE:::.robustConditionTest(s$b[1:10, ], s$v[1:10, ], tau[1:10], trt))
  expect_true(attr(mod, "moderated"))
})

test_that("every arm is identical serially and in parallel", {
  skip_on_os("windows")
  for (a in list(spiDE:::.armSpec("floor", "uncapped", "moderated"), spiDE:::.armSpec(condition_variance = "moderated"))) {
    r1 <- testSpiDE(fit16a, condition = "condition", .arm = a)@table
    r2 <- testSpiDE(fit16a, condition = "condition", .arm = a, BPPARAM = BiocParallel::MulticoreParam(2))@table
    expect_identical(r1$p, r2$p)
  }
})
