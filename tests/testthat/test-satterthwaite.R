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

test_that(".satterthwaiteDF matches lmerTest with a random slope", {
  skip_if_not_installed("lmerTest"); skip_if_not_installed("lme4")
  set.seed(3)
  n_g <- 30L; n_per <- 12L; n <- n_g * n_per
  g <- factor(rep(seq_len(n_g), each = n_per))
  x <- rnorm(n)
  u0 <- rnorm(n_g, 0, 0.7); u1 <- rnorm(n_g, 0, 0.5)
  y <- 1 + 0.4 * x + u0[as.integer(g)] + u1[as.integer(g)] * x + rnorm(n, 0, 1)
  m <- lmerTest::lmer(y ~ x + (1 + x || g), REML = TRUE)   # independent int+slope
  df_lmer <- summary(m)$coefficients["x", "df"]
  vc <- as.data.frame(lme4::VarCorr(m))
  t_int <- vc$vcov[vc$grp == "g" & is.na(vc$var2) & vc$var1 == "(Intercept)"]
  t_slp <- vc$vcov[grepl("g", vc$grp) & vc$var1 == "x"]
  sig2  <- vc$vcov[vc$grp == "Residual"]
  Zi <- stats::model.matrix(~ 0 + g); Zs <- Zi * x
  W  <- cbind(`(Intercept)` = 1, x = x, Zi, Zs)
  reg <- c(NA, NA, rep("SampleInt", ncol(Zi)), rep("SampleSlope", ncol(Zs)))
  wbar <- rep(1 / sig2, n)
  pen  <- c(0, 0, rep(1 / t_int, ncol(Zi)), rep(1 / t_slp, ncol(Zs)))
  A <- crossprod(W * sqrt(wbar)); minv <- SpaNorm::invert_mat(A + diag(pen))
  df_s <- spiDE:::.satterthwaiteDF(A, minv, pen, reg,
            c(SampleInt = t_int, SampleSlope = t_slp),
            tested = 2L, ncells = n, tested_names = "x")
  expect_equal(unname(df_s), unname(df_lmer), tolerance = 0.08)
})

test_that(".varParamCov matches a brute-force n x n REML information (2 RE groups)", {
  set.seed(3)
  n_g <- 12L; n_per <- 8L; n <- n_g * n_per
  g <- factor(rep(seq_len(n_g), each = n_per)); x <- rnorm(n)
  Zi <- stats::model.matrix(~ 0 + g); Zs <- Zi * x
  W <- cbind(`(Intercept)` = 1, x = x, Zi, Zs)
  reg <- c(NA, NA, rep("SampleInt", ncol(Zi)), rep("SampleSlope", ncol(Zs)))
  t_int <- 0.5; t_slp <- 0.3; sig2 <- 1
  wbar <- rep(1 / sig2, n)
  pen <- c(0, 0, rep(1 / t_int, ncol(Zi)), rep(1 / t_slp, ncol(Zs)))
  A <- crossprod(W * sqrt(wbar)); minv <- SpaNorm::invert_mat(A + diag(pen))
  vp <- spiDE:::.varParamCov(A, minv, pen, reg, n)
  # brute-force REML expected information I_ab = 0.5 tr(P V_a P V_b), phi = 1,
  # P = W-bar - W-bar C minv C' W-bar; param order (phi, groups...)
  Wd <- diag(wbar)
  P <- Wd - Wd %*% W %*% minv %*% t(W) %*% Wd
  Vlist <- c(list(diag(1 / wbar)),
             lapply(vp$groups, function(gr) {
               Z <- W[, which(reg == gr), drop = FALSE]; Z %*% t(Z)
             }))
  K <- length(Vlist); Ibrute <- matrix(0, K, K)
  for (a in seq_len(K)) for (b in seq_len(K))
    Ibrute[a, b] <- 0.5 * sum(diag(P %*% Vlist[[a]] %*% P %*% Vlist[[b]]))
  expect_equal(unname(vp$cov), unname(solve(Ibrute)), tolerance = 1e-6)
})

test_that("df.method='satterthwaite' yields a per-column @df; Response df ~ S-2", {
  spe <- buildNiches(spiDE:::.toyClustered(n_samples = 16, sd_patient = 0.7),
                     sigma = 30)
  fb <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                 df.method = "between", verbose = FALSE)
  fs <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                 df.method = "satterthwaite", verbose = FALSE)
  db <- fits(fb)[[1]]@df
  ds <- fits(fs)[[1]]@df
  # between: scalar S-2 = 14; satterthwaite: one df per tested column
  expect_length(db, 1L)
  ct <- as.character(fits(fs)[[1]]@covtype)
  n_tested <- sum(grepl("Response", ct))
  expect_length(ds, n_tested)
  # Response column df ~ S-2 (the 0.99x regression anchor)
  cm <- fits(fs)[[1]]@coefmap
  resp_name <- cm$covariate[ct == "Response"]
  expect_equal(unname(ds[resp_name]), 14, tolerance = 0.20)   # 16 samples -> ~14
  # niche-interaction dfs are larger than the Response df (within-sample info)
  rn_names <- cm$covariate[ct == "ResponseNiche"]
  expect_gt(stats::median(ds[rn_names]), ds[resp_name])
})

test_that("df.method='between' is byte-identical to the pre-existing scalar df", {
  spe <- buildNiches(spiDE:::.toyClustered(n_samples = 12, sd_patient = 0.7),
                     sigma = 30)
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  fb <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                 df.method = "between", verbose = FALSE)
  expect_equal(fits(fb)[[1]]@df, 10)                        # 12 samples -> S-2
  tb <- spiDE:::.blockedInference(fits(fb)[[1]], Y)
  expect_false(any(is.na(tb@p.combined.pos)))
})

test_that("df.method='satterthwaite' works for random='slope'", {
  spe <- buildNiches(spiDE:::.toyClustered(n_samples = 16, sd_patient = 0.7),
                     sigma = 30)
  fs <- fitSpiDE(spe, "condition", sigma = 30, random = "slope",
                 df.method = "satterthwaite", verbose = FALSE)
  ds <- fits(fs)[[1]]@df
  ct <- as.character(fits(fs)[[1]]@covtype)
  expect_length(ds, sum(grepl("Response", ct)))
  expect_true(all(is.finite(ds)) && all(ds >= 1))
})
