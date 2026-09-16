# The two kernels .polishBatch()'s Newton is built on -- the mean and the
# penalised NB log-likelihood -- on either backend.
#
# Everything else in the loop (the line search, the dispersion search, the
# convergence test) is a consumer of these two, so they go to tensors first.
# base R is the oracle: these are the expressions already in .polishBatch().

.nbFixture <- function(n = 50, p = 4, b = 5, seed = 17) {
  set.seed(seed)
  W <- cbind(1, matrix(stats::rnorm(n * (p - 1)), n, p - 1))
  A <- matrix(stats::rnorm(b * p, 0, 0.3), b, p)
  A[, 1] <- A[, 1] + 1.5
  mu <- exp(A %*% t(W))
  Y <- matrix(stats::rnbinom(b * n, mu = as.numeric(mu), size = 2), b, n)
  list(W = W, A = A, Y = Y, psi = stats::runif(b, 0.1, 1.5),
       pen = stats::runif(p, 0, 0.4))
}

# the expressions verbatim from .polishBatch()
.muRef <- function(f) pmax(exp(f$A %*% t(f$W)), spiDE:::.MU_FLOOR)
.llRef <- function(f, M) {
  rowSums(stats::dnbinom(f$Y, size = 1 / f$psi, mu = M, log = TRUE)) -
    0.5 * as.numeric((f$A^2) %*% f$pen)
}

test_that(".muBatch reproduces the base-R mean, including the floor", {
  f <- .nbFixture()
  expect_equal(spiDE:::.muBatch(f$A, f$W), .muRef(f), tolerance = 1e-12)

  # the floor has to bite, or it is not being tested. Drive the INTERCEPT down
  # and leave the rest at zero: setting every coefficient to -50 gives
  # eta = -50 * rowSums(W), which is large and POSITIVE wherever rowSums(W) is
  # negative, so most cells would not be floored at all.
  f2 <- f; f2$A[1, ] <- c(-50, rep(0, ncol(f$A) - 1L))
  got <- spiDE:::.muBatch(f2$A, f2$W)
  expect_true(all(got[1, ] == spiDE:::.MU_FLOOR))
  expect_equal(got, .muRef(f2), tolerance = 1e-12)
})

test_that(".nbLoglikBatch reproduces dnbinom's penalised row sums", {
  f <- .nbFixture()
  M <- .muRef(f)
  expect_equal(spiDE:::.nbLoglikBatch(f$Y, M, f$psi, f$A, f$pen), .llRef(f, M),
               tolerance = 1e-10)
})

test_that("the NB kernels agree between the base-R and torch branches", {
  skip_if_not_installed("torch")
  f <- .nbFixture()
  M <- .muRef(f)
  ll <- .llRef(f, M)

  tt <- function(x) torch::torch_tensor(x, dtype = torch::torch_float64())
  Wt <- tt(f$W); At <- tt(f$A); Yt <- tt(f$Y)
  Mt <- spiDE:::.muBatch(At, Wt)
  expect_true(SpaNorm::is_torch_tensor(Mt))
  expect_equal(as.matrix(SpaNorm::toRMatrix(Mt)), M, tolerance = 1e-10)

  llt <- spiDE:::.nbLoglikBatch(Yt, Mt, f$psi, At, f$pen)
  expect_equal(as.numeric(SpaNorm::toRMatrix(llt)), ll, tolerance = 1e-9)
})

test_that(".nbLoglikBatch handles an all-zero gene and a large count", {
  # the two ends .polishBatch()'s own tests exercise per gene
  f <- .nbFixture()
  f$Y[1, ] <- 0L
  f$Y[2, ] <- f$Y[2, ] * 100L
  M <- .muRef(f)
  got <- spiDE:::.nbLoglikBatch(f$Y, M, f$psi, f$A, f$pen)
  expect_equal(got, .llRef(f, M), tolerance = 1e-10)
  expect_true(all(is.finite(got)))
})
