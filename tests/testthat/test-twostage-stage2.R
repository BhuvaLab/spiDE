test_that("SpiDEResults carries a diagnostics list slot", {
  r <- new("SpiDEResults", fits = list(), sigma = 30, condition = "condition",
           mode = "condition", index = "A", niche = "B",
           covariates = character(), coldata = S4Vectors::DataFrame(),
           gene.weights = NULL, p.cauchy.pos = NULL, p.cauchy.neg = NULL,
           results = data.frame(), fdr = 0.05, call = NULL,
           diagnostics = list(r2 = data.frame()))
  expect_identical(names(r@diagnostics), "r2")
  # default construction still works and defaults to an empty list
  r0 <- new("SpiDEResults", fits = list(), sigma = 30, condition = "c",
            mode = "condition", index = "A", niche = "B",
            covariates = character(), coldata = S4Vectors::DataFrame(),
            gene.weights = NULL, p.cauchy.pos = NULL, p.cauchy.neg = NULL,
            results = data.frame(), fdr = 0.05, call = NULL)
  expect_identical(r0@diagnostics, list())
})

test_that(".poolPatientSlopes precision-pools cores within a patient", {
  B <- array(NA_real_, c(1, 1, 3),
             dimnames = list("g", "n", c("s1", "s2", "s3")))
  V <- B
  B[1, 1, ] <- c(1, 3, 10); V[1, 1, ] <- c(1, 0.5, 2)
  s2p <- c(s1 = "p1", s2 = "p1", s3 = "p2")
  pooled <- spiDE:::.poolPatientSlopes(B, V, s2p)
  expect_identical(dimnames(pooled$beta)[[3]], c("p1", "p2"))
  expect_equal(pooled$beta["g", "n", "p1"], (1 / 1 + 3 / 0.5) / (1 + 2))
  expect_equal(pooled$var["g", "n", "p1"], 1 / 3)
  expect_equal(pooled$beta["g", "n", "p2"], 10)   # single core passes through
  # an NA core is ignored, not contagious
  B[1, 1, 2] <- NA
  pooled2 <- spiDE:::.poolPatientSlopes(B, V, s2p)
  expect_equal(pooled2$beta["g", "n", "p1"], 1)
})

test_that(".tau2DL recovers heterogeneity and floors at zero", {
  set.seed(11)
  S <- 40; G <- 200
  X <- cbind(1, rep(0:1, each = S / 2))
  v <- matrix(runif(G * S, 0.05, 0.1), G, S)
  # homogeneous slopes: tau2 ~ 0
  B0 <- matrix(rnorm(G * S, sd = sqrt(v)), G, S)
  expect_lt(spiDE:::.tau2DL(B0, v, X), 0.02)
  # true tau2 = 0.5 added between patients
  u <- matrix(rnorm(G * S, sd = sqrt(0.5)), G, S)
  B1 <- B0 + u
  t2 <- spiDE:::.tau2DL(B1, v, X)
  expect_gt(t2, 0.25); expect_lt(t2, 1)
})
