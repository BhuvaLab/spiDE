# Per-gene convergence after fitNB. SpaNorm's multi-gene IRLS shares one cell
# weight vector across genes and stops on the aggregate log-likelihood, so
# bright genes sit 1-4 SE from their own optimum with a dispersion ~1.6x too
# large. These tests pin the polish stage: it never decreases the penalised
# log-likelihood, it reaches a stationary point, it leaves an already-converged
# gene alone, and the Schur-complement solve equals the dense one.

set.seed(11)

# a small design with an indicator block, in the shape .polishFit() sees
toy_design <- function(n = 240, p_x = 6, n_grp = 8) {
  X <- cbind(1, matrix(rnorm(n * (p_x - 1)), n, p_x - 1))
  colnames(X) <- c("(Intercept)", paste0("x", seq_len(p_x - 1)))
  gidx <- rep(seq_len(n_grp), length.out = n)
  Z <- stats::model.matrix(~ 0 + factor(gidx))
  colnames(Z) <- paste0("SampleCellType", seq_len(n_grp))
  W <- cbind(X, Z)
  list(W = W, nested = grepl("^SampleCellType", colnames(W)),
       pen = c(rep(0, p_x), rep(1.7, n_grp)))
}

test_that(".newtonSolver's Schur solve equals the dense solve", {
  d <- toy_design()
  w <- runif(nrow(d$W), 0.2, 3)
  s <- rnorm(ncol(d$W))

  info <- crossprod(d$W * sqrt(w))
  diag(info) <- diag(info) + d$pen
  expect_equal(as.numeric(spiDE:::.newtonSolver(d$W, d$pen, d$nested)$solve(w, s)),
               as.numeric(solve(info, s)), tolerance = 1e-8)
})

test_that(".newtonSolver's xcov equals the X-block of the dense covariance", {
  d <- toy_design()
  w <- runif(nrow(d$W), 0.2, 3)
  info <- crossprod(d$W * sqrt(w))
  diag(info) <- diag(info) + d$pen
  expect_equal(spiDE:::.newtonSolver(d$W, d$pen, d$nested)$xcov(w),
               solve(info)[!d$nested, !d$nested], tolerance = 1e-8,
               ignore_attr = TRUE)
})

test_that(".newtonSolver falls back to the dense path with no nested block", {
  d <- toy_design()
  nested <- rep(FALSE, ncol(d$W))
  w <- runif(nrow(d$W), 0.2, 3)
  s <- rnorm(ncol(d$W))
  info <- crossprod(d$W * sqrt(w))
  diag(info) <- diag(info) + d$pen
  expect_equal(as.numeric(spiDE:::.newtonSolver(d$W, d$pen, nested)$solve(w, s)),
               as.numeric(solve(info, s)), tolerance = 1e-8)
})

test_that(".polishGene reaches a stationary point and raises the log-likelihood", {
  d <- toy_design()
  a_true <- c(1.2, 0.4, -0.3, 0.2, 0, 0.1, rnorm(8, 0, 0.3))
  mu <- exp(d$W %*% a_true)
  y <- rnbinom(length(mu), mu = as.numeric(mu), size = 1 / 0.4)
  a0 <- c(0.2, rep(0, ncol(d$W) - 1))          # a deliberately poor start
  solver <- spiDE:::.newtonSolver(d$W, d$pen, d$nested)

  ll0 <- spiDE:::.nbPenLoglik(y, as.numeric(exp(d$W %*% a0)), 0.4, a0, d$pen)
  out <- spiDE:::.polishGene(y, d$W, a0, 0.4, d$pen, solver,
                             maxit = 50L, tol = 1e-10)

  expect_gt(out$loglik, ll0)
  # penalised score at the polished point, at the polished psi
  mu1 <- as.numeric(exp(d$W %*% out$alpha))
  s1 <- as.numeric(crossprod(d$W, (y - mu1) / (1 + out$psi * mu1))) -
    d$pen * out$alpha
  expect_lt(max(abs(s1)), 1e-4 * max(abs(y)))
  expect_true(is.finite(out$psi) && out$psi > 0)
  expect_gte(out$iterations, 1L)
  expect_false(out$singular)
})

test_that(".polishGene leaves an already-converged gene alone", {
  d <- toy_design()
  a_true <- c(1.0, 0.3, -0.2, 0.1, 0, 0, rnorm(8, 0, 0.2))
  mu <- exp(d$W %*% a_true)
  y <- rnbinom(length(mu), mu = as.numeric(mu), size = 1 / 0.5)
  solver <- spiDE:::.newtonSolver(d$W, d$pen, d$nested)
  # converge once, then polish again from the converged point
  first <- spiDE:::.polishGene(y, d$W, a_true, 0.5, d$pen, solver, 50L, 1e-12)
  again <- spiDE:::.polishGene(y, d$W, first$alpha, first$psi, d$pen, solver,
                               50L, 1e-12)
  expect_equal(again$alpha, first$alpha, tolerance = 1e-3)
  expect_equal(again$psi, first$psi, tolerance = 1e-3)
})

test_that(".polishGene restarts from cell-type means when the start is degenerate", {
  d <- toy_design()
  mu <- exp(d$W %*% c(1.0, rep(0.1, 5), rnorm(8, 0, 0.2)))
  y <- rnbinom(length(mu), mu = as.numeric(mu), size = 1 / 0.4)
  a_bad <- c(-25, rep(0, ncol(d$W) - 1))       # min log mu < -10
  solver <- spiDE:::.newtonSolver(d$W, d$pen, d$nested)
  out <- spiDE:::.polishGene(y, d$W, a_bad, 0.4, d$pen, solver, 50L, 1e-10)

  expect_true(out$restarted)
  expect_gt(min(as.numeric(d$W %*% out$alpha)), -10)
  expect_true(is.finite(out$loglik))
})

test_that(".polishFit is invariant to gene blocking", {
  d <- toy_design()
  ng <- 5
  A0 <- matrix(0, ng, ncol(d$W), dimnames = list(paste0("G", seq_len(ng)),
                                                 colnames(d$W)))
  A0[, 1] <- 0.5
  Y <- t(vapply(seq_len(ng), function(g) {
    mu <- exp(d$W %*% c(1 + 0.2 * g, rep(0.15, 5), rnorm(8, 0, 0.2)))
    rnbinom(nrow(d$W), mu = as.numeric(mu), size = 1 / 0.4)
  }, numeric(nrow(d$W))))
  dimnames(Y) <- list(rownames(A0), NULL)
  psi0 <- rep(0.4, ng)
  re_group <- ifelse(d$nested, "SampleCellTypeInt", NA_character_)

  a <- spiDE:::.polishFit(Y, d$W, A0, psi0, d$pen, re_group,
                          block.size = NULL)
  b <- spiDE:::.polishFit(Y, d$W, A0, psi0, d$pen, re_group,
                          block.size = 2L)
  expect_equal(a$alpha, b$alpha)
  expect_equal(a$psi, b$psi)
  expect_equal(nrow(a$polish), ng)
  expect_setequal(colnames(a$polish),
                  c("iterations", "psi_fitnb", "restarted", "capped", "singular"))
  expect_equal(a$polish$psi_fitnb, psi0)
})

test_that(".polishFit flags a degenerate gene and keeps finite values", {
  d <- toy_design()
  # a gene that is zero everywhere: every working weight is ~0, so the
  # information matrix is singular at the start
  Y <- matrix(0L, 1, nrow(d$W), dimnames = list("G1", NULL))
  A0 <- matrix(-40, 1, ncol(d$W), dimnames = list("G1", colnames(d$W)))
  re_group <- ifelse(d$nested, "SampleCellTypeInt", NA_character_)
  out <- spiDE:::.polishFit(Y, d$W, A0, 0.4, d$pen, re_group)
  expect_true(out$polish$singular || out$polish$restarted)
  expect_true(all(is.finite(out$alpha)))
  expect_true(is.finite(out$psi))
})
