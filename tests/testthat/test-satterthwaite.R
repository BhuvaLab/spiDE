test_that(".ptByCol handles NULL / scalar / vector / matrix df", {
  tm <- matrix(c(0, 1, -1, 2), nrow = 2)         # 2 genes x 2 cols
  # NULL -> normal
  expect_equal(spiDE:::.ptByCol(tm, NULL), stats::pnorm(tm))
  # scalar -> same df everywhere
  expect_equal(spiDE:::.ptByCol(tm, 5), stats::pt(tm, df = 5))
  # per-column vector -> col 1 uses df=5, col 2 uses df=50
  ref <- cbind(stats::pt(tm[, 1], df = 5), stats::pt(tm[, 2], df = 50))
  expect_equal(spiDE:::.ptByCol(tm, c(5, 50)), ref)
  # single-gene vector input with per-column df
  expect_equal(spiDE:::.ptByCol(c(0, 2), c(5, 50)),
               c(stats::pt(0, 5), stats::pt(2, 50)))
  # lower.tail forwarded
  expect_equal(spiDE:::.ptByCol(tm, 5, lower.tail = FALSE),
               stats::pt(tm, df = 5, lower.tail = FALSE))
})

test_that(".satterthwaiteDF matches lmerTest df on a Gaussian random-intercept LMM", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("lme4")
  set.seed(1)
  n_g <- 24L; n_per <- 8L; n <- n_g * n_per
  g <- factor(rep(seq_len(n_g), each = n_per))
  x <- rnorm(n)                                  # varies WITHIN group
  u <- rnorm(n_g, 0, 0.8)
  y <- 1 + 0.5 * x + u[as.integer(g)] + rnorm(n, 0, 1.0)
  m  <- lmerTest::lmer(y ~ x + (1 | g), REML = TRUE)
  df_lmer <- summary(m)$coefficients["x", "df"]
  vc   <- as.data.frame(lme4::VarCorr(m))
  tau2 <- vc$vcov[vc$grp == "g"]
  sig2 <- vc$vcov[vc$grp == "Residual"]
  # spiDE standardised working-model inputs: weights = precision, phi = 1
  Z <- stats::model.matrix(~ 0 + g)
  W <- cbind(`(Intercept)` = 1, x = x, Z)
  reg <- c(NA, NA, rep("g", ncol(Z)))
  wbar <- rep(1 / sig2, n)
  pen  <- c(0, 0, rep(1 / tau2, ncol(Z)))
  A <- crossprod(W * sqrt(wbar))
  minv <- SpaNorm::invert_mat(A + diag(pen))
  df_s <- spiDE:::.satterthwaiteDF(A, minv, pen, reg,
                                   stats::setNames(tau2, "g"),
                                   tested = 2L, ncells = n, tested_names = "x")
  expect_equal(unname(df_s), unname(df_lmer), tolerance = 0.05)
})

test_that(".satterthwaiteDF gives large df for a within-group contrast, small for between", {
  set.seed(2)
  n_g <- 20L; n_per <- 10L; n <- n_g * n_per
  g <- factor(rep(seq_len(n_g), each = n_per))
  xin  <- rnorm(n)                                   # within-group
  xbtw <- rnorm(n_g)[as.integer(g)]                  # constant within group
  Z <- stats::model.matrix(~ 0 + g)
  W <- cbind(`(Intercept)` = 1, xin = xin, xbtw = xbtw, Z)
  reg <- c(NA, NA, NA, rep("g", ncol(Z)))
  tau2 <- 0.5; pen <- c(0, 0, 0, rep(1 / tau2, ncol(Z)))
  A <- crossprod(W); minv <- SpaNorm::invert_mat(A + diag(pen))
  df <- spiDE:::.satterthwaiteDF(A, minv, pen, reg, c(g = tau2),
                                 tested = c(2L, 3L), ncells = n,
                                 tested_names = c("xin", "xbtw"))
  expect_gt(df[["xin"]], 5 * df[["xbtw"]])           # within >> between
  expect_lt(df[["xbtw"]], n_g)                        # between ~ O(n_g), not O(n)
})
