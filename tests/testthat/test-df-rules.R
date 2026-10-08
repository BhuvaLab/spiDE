# The pooled test's df rules of feature/small-sample-df (research/smalldf/README.md),
# checked against direct limma computations; the default rule is the shipped test.

spe16d <- buildNiches(.toySPE(n_samples = 16, n_per = 120, n_genes = 12, seed = 2), sigma = 30, verbose = FALSE)
fit16d <- fitSpiDE(spe16d, index = "A", sigma = 30, verbose = FALSE)

# a synthetic column family: G genes x S patients, unequal patient precision
dfFamily <- function(G = 150, S = 10, seed = 21) {
  set.seed(seed)
  mean_expr <- exp(seq(-3, 3, length.out = G))
  v <- matrix(stats::rgamma(G * S, 2, 2 / 0.05), G, S)
  b <- matrix(stats::rnorm(G * S, 0.05, sqrt(v + 0.01)), G, S)
  b[1, 1:2] <- NA
  list(b = b, v = v, tau = rep(0.01, G), mean_expr = mean_expr)
}
directFit <- function(s) {
  w <- 1 / (s$v + s$tau); w[!is.finite(s$b)] <- NA
  w <- w / rowMeans(w, na.rm = TRUE)
  fit <- limma::lmFit(s$b, matrix(1, ncol(s$b), 1), weights = w)
  fit$Amean <- log(s$mean_expr + 1e-3)
  ww <- w; ww[!is.finite(ww) | !is.finite(s$b)] <- NA
  list(fit = limma::eBayes(fit, robust = TRUE, trend = TRUE),
       neff = rowSums(ww, na.rm = TRUE)^2 / rowSums(ww^2, na.rm = TRUE), m = rowSums(is.finite(ww)))
}

test_that("the default arm is the shipped test", {
  expect_identical(spiDE:::.armSpec(), list(pooled_df = "capped"))
  r0 <- testSpiDE(fit16d, condition = "condition")@table
  r1 <- testSpiDE(fit16d, condition = "condition", .arm = spiDE:::.armSpec())@table
  expect_identical(r0, r1)
})

test_that("every df rule equals its direct limma computation, and only the df move", {
  s <- dfFamily(); d <- directFit(s)
  get <- function(rule) spiDE:::.pooledColumnTest(s$b, s$v, s$tau, s$mean_expr, df_rule = rule)
  cap <- get("capped"); pro <- get("proportional"); kpp <- get("kish_plus_prior"); lim <- get("limma")
  ok <- rowSums(is.finite(s$b)) >= 6
  expect_equal(cap$df[ok], pmin(d$fit$df.total, pmax(d$neff - 1, 1))[ok], tolerance = 1e-12)
  expect_equal(pro$df[ok], (d$fit$df.total * d$neff / d$m)[ok], tolerance = 1e-12)
  expect_equal(kpp$df[ok], pmin(d$fit$df.total, pmax(d$neff - 1, 1) + d$fit$df.prior)[ok], tolerance = 1e-12)
  expect_equal(lim$df[ok], d$fit$df.total[ok], tolerance = 1e-12)
  for (r in list(pro, kpp, lim)) { expect_identical(r$t, cap$t); expect_identical(r$estimate, cap$estimate) }
  # the capped df is the smallest and limma's the largest
  expect_true(all(cap$df[ok] <= kpp$df[ok] + 1e-12) && all(kpp$df[ok] <= lim$df[ok] + 1e-12))
  expect_true(all(pro$df[ok] <= lim$df[ok] + 1e-12))
})

test_that("every arm's pooled test and filter never look at the condition", {
  pt <- fit16d@patients
  set.seed(4)
  fp <- fit16d
  fp@patients$perm <- factor(stats::setNames(sample(as.character(pt$condition)), pt$patient)[pt$patient])
  for (r in c("proportional", "kish_plus_prior", "limma")) {
    a <- spiDE:::.armSpec(r)
    r1 <- testSpiDE(fit16d, condition = "condition", .arm = a)@table
    r2 <- testSpiDE(fp, condition = "perm", .arm = a)@table
    expect_identical(r1[r1$test == "pooled", c("gene", "niche", "p", "q")], r2[r2$test == "pooled", c("gene", "niche", "p", "q")])
    # the condition test is unchanged by the pooled df rule (only its filter can move)
    r0 <- testSpiDE(fit16d, condition = "condition")@table
    expect_identical(r1$p[r1$test == "condition"], r0$p[r0$test == "condition"])
  }
})

test_that("every arm is identical serially and in parallel", {
  skip_on_os("windows")
  for (r in c("proportional", "kish_plus_prior")) {
    a <- spiDE:::.armSpec(r)
    r1 <- testSpiDE(fit16d, condition = "condition", .arm = a)@table
    r2 <- testSpiDE(fit16d, condition = "condition", .arm = a, BPPARAM = BiocParallel::MulticoreParam(2))@table
    expect_identical(r1$p, r2$p)
  }
})
