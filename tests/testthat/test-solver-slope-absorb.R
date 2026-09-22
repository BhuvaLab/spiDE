# Absorbing a random block that is BLOCK-diagonal by sample, not merely diagonal.
#
# .newtonSolver()'s Schur absorption was written for the nested (sample x cell
# type) intercept: 0/1 indicators partitioning the cells, so C = Z' diag(w) Z is
# diagonal and C^-1 is a reciprocal. Random SLOPES are indicator x covariate, so
# a sample's slope columns are not orthogonal to each other or to that sample's
# intercept -- C is block-diagonal by sample with a dense block per sample, and
# the scalar path cannot absorb it.
#
# It is still absorbable, and by the same identity: every random column belongs
# to exactly one sample, so cells of sample s load on no other sample's columns.
# The only change is that C^-1 is a per-sample dense solve rather than a
# reciprocal.
#
# The oracle throughout is the DENSE solve on the full design. Absorption is an
# exact re-arrangement, so "exact" is the standard, not "close".

# A design whose random block is block-diagonal by sample and NOT a 0/1
# partition: per-sample intercept plus per-sample slopes on K continuous bases.
.slopeFixture <- function(n = 120, px = 3, S = 4, K = 2, seed = 11) {
  set.seed(seed)
  X <- cbind(1, matrix(stats::rnorm(n * (px - 1)), n, px - 1))
  smp <- factor(c(seq_len(S), sample(seq_len(S), n - S, replace = TRUE)))
  Zint <- stats::model.matrix(~ 0 + smp)
  bases <- matrix(stats::rnorm(n * K), n, K)          # the CellType:niche bases
  Zsl <- do.call(cbind, lapply(seq_len(K), function(k) Zint * bases[, k]))
  Z <- cbind(Zint, Zsl)
  # block id per column: every random column belongs to exactly one sample
  blk <- c(as.integer(smp[!duplicated(smp)][order(unique(as.integer(smp)))]),
           rep(NA_integer_, 0))
  blk <- c(seq_len(S), rep(seq_len(S), times = K))
  list(W = cbind(X, Z), px = px, S = S, K = K,
       group = c(rep(NA_integer_, px), blk),
       pen = c(rep(0, px), rep(0.4, ncol(Z))),
       w = stats::runif(n, 0.2, 2),
       score = stats::rnorm(px + ncol(Z)))
}

# the dense oracle: (X'WX + diag(pen))^-1 s, no absorption anywhere
.denseSolve <- function(W, pen, w, s) {
  info <- crossprod(W * sqrt(w))
  diag(info) <- diag(info) + pen
  as.numeric(solve(info, s))
}

test_that(".newtonSolver absorbs a per-sample block-diagonal random block exactly", {
  f <- .slopeFixture()
  # sanity: this really is NOT the case the scalar path handles -- a sample's
  # slope column is not orthogonal to that sample's intercept
  zi <- which(!is.na(f$group))
  C <- crossprod(f$W[, zi, drop = FALSE] * sqrt(f$w))
  expect_gt(max(abs(C[upper.tri(C)])), 1e-6)

  sol <- spiDE:::.newtonSolver(f$W, f$pen, f$group)
  got <- sol$solve(f$w, f$score)
  expect_equal(got, .denseSolve(f$W, f$pen, f$w, f$score), tolerance = 1e-8)
})

test_that(".newtonSolver's xcov on a block-diagonal random block is the fixed-effect covariance", {
  f <- .slopeFixture()
  sol <- spiDE:::.newtonSolver(f$W, f$pen, f$group)
  got <- sol$xcov(f$w)
  info <- crossprod(f$W * sqrt(f$w)); diag(info) <- diag(info) + f$pen
  # the Schur complement's inverse IS the fixed-effect block of the full inverse
  expect_equal(got, solve(info)[seq_len(f$px), seq_len(f$px)], tolerance = 1e-8)
})

test_that("a logical `nested` still absorbs the indicator block exactly", {
  # the existing call sites pass a logical vector; that path must not move
  set.seed(5)
  n <- 90; px <- 3; G <- 5
  X <- cbind(1, matrix(stats::rnorm(n * (px - 1)), n, px - 1))
  g <- c(seq_len(G), sample(seq_len(G), n - G, replace = TRUE))
  Z <- matrix(0, n, G); Z[cbind(seq_len(n), g)] <- 1
  W <- cbind(X, Z); pen <- c(rep(0, px), rep(0.3, G))
  nested <- c(rep(FALSE, px), rep(TRUE, G))
  w <- stats::runif(n, 0.2, 2); s <- stats::rnorm(px + G)

  sol <- spiDE:::.newtonSolver(W, pen, nested)
  expect_equal(sol$solve(w, s), .denseSolve(W, pen, w, s), tolerance = 1e-8)
})

test_that("a block grouping that mixes indicators and slopes absorbs exactly", {
  # the production shape: SampleInt + SampleSlope + SampleCellTypeInt, all of
  # which belong to one sample, so all of them go in that sample's block
  set.seed(7)
  n <- 150; px <- 3; S <- 3; K <- 2; nct <- 2
  X <- cbind(1, matrix(stats::rnorm(n * (px - 1)), n, px - 1))
  smp <- factor(sample(seq_len(S), n, replace = TRUE))
  ct <- factor(sample(seq_len(nct), n, replace = TRUE))
  Zint <- stats::model.matrix(~ 0 + smp)
  bases <- matrix(stats::rnorm(n * K), n, K)
  Zsl <- do.call(cbind, lapply(seq_len(K), function(k) Zint * bases[, k]))
  grp <- interaction(smp, ct, drop = TRUE)
  Zct <- stats::model.matrix(~ 0 + grp)
  ct_sample <- as.integer(sub("\\..*$", "", levels(grp)))

  W <- cbind(X, Zint, Zsl, Zct)
  pen <- c(rep(0, px), rep(0.5, ncol(Zint) + ncol(Zsl) + ncol(Zct)))
  group <- c(rep(NA_integer_, px), seq_len(S), rep(seq_len(S), times = K), ct_sample)
  w <- stats::runif(n, 0.2, 2); s <- stats::rnorm(ncol(W))

  sol <- spiDE:::.newtonSolver(W, pen, group)
  expect_equal(sol$solve(w, s), .denseSolve(W, pen, w, s), tolerance = 1e-8)
})

test_that("the batched solver refuses a multi-column block rather than guessing", {
  # .absorbBatch()/.newtonSolverBatch() carry their own copy of the absorption
  # and only implement the diagonal (1x1 block) case -- C^-1 is a reciprocal
  # there, and a per-block Cholesky has no batched equivalent written yet. A
  # grouping they cannot honour must stop with a message that says so, because
  # the failure mode otherwise is a wrong Newton step, not an error.
  f <- .slopeFixture()
  expect_error(spiDE:::.newtonSolverBatch(f$W, f$pen, f$group),
               "1x1|one column|block", ignore.case = TRUE)
})

test_that("the batched solver still accepts the logical indicator case", {
  set.seed(5)
  n <- 60; px <- 3; G <- 4
  X <- cbind(1, matrix(stats::rnorm(n * (px - 1)), n, px - 1))
  g <- c(seq_len(G), sample(seq_len(G), n - G, replace = TRUE))
  Z <- matrix(0, n, G); Z[cbind(seq_len(n), g)] <- 1
  W <- cbind(X, Z); pen <- c(rep(0, px), rep(0.3, G))
  nested <- c(rep(FALSE, px), rep(TRUE, G))
  expect_no_error(spiDE:::.newtonSolverBatch(W, pen, nested))
})
