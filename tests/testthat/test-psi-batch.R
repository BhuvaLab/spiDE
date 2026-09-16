# The profile dispersion search, batched and fixed-iteration.
#
# This is the ONE deliberate numerical divergence in the batched fitter, and
# the justification is accuracy rather than speed: optimize()'s default
# tolerance is .Machine$double.eps^0.25, about 1.2e-4 ABSOLUTE on a
# log-interval of width log(1e3) - log(1e-3) = 13.8, so the psi the per-gene
# engine reports is only determined to ~1e-4 anyway. Fifty bisection steps on
# the score take that interval to ~1e-14.
#
# So the gate is not "the same answer as optimize()". It is "at least as good
# an answer", measured on the objective both are maximising, plus agreement
# within optimize()'s own tolerance. Bisection is also fixed-iteration, which
# means no divergent control flow across a batch -- the reason the plan chose
# it over reproducing Brent.

.psiFixture <- function(n = 200, b = 6, seed = 5) {
  set.seed(seed)
  mu <- matrix(exp(stats::rnorm(b * n, 2, 0.3)), b, n)
  Y <- matrix(0L, b, n)
  # genes 1..b-1 overdispersed at a range of psi; gene b is exactly Poisson,
  # whose NB optimum sits below the search range and must come back at_bound
  psis <- c(0.05, 0.2, 0.8, 2, 10)[seq_len(b - 1L)]
  for (i in seq_len(b - 1L)) {
    Y[i, ] <- stats::rnbinom(n, mu = mu[i, ], size = 1 / psis[i])
  }
  Y[b, ] <- stats::rpois(n, mu[b, ])
  list(Y = Y, Mu = mu, range = c(1e-3, 1e3))
}

# what the per-gene engine does today, verbatim
.psiRef <- function(f) {
  lo <- log(f$range[1]); hi <- log(f$range[2])
  est <- numeric(nrow(f$Y)); bnd <- logical(nrow(f$Y))
  for (i in seq_len(nrow(f$Y))) {
    y <- f$Y[i, ]; mu <- f$Mu[i, ]
    o <- stats::optimize(function(lp) {
      -sum(stats::dnbinom(y, size = 1 / exp(lp), mu = mu, log = TRUE))
    }, c(lo, hi))
    bnd[i] <- (o$minimum - lo) < 1e-3 * (hi - lo) ||
      (hi - o$minimum) < 1e-3 * (hi - lo)
    est[i] <- exp(o$minimum)
  }
  list(psi = est, at_bound = bnd)
}

.nbll <- function(y, mu, psi) sum(stats::dnbinom(y, size = 1 / psi, mu = mu, log = TRUE))

test_that(".psiProfileBatch agrees with optimize() within optimize's own tolerance", {
  f <- .psiFixture()
  ref <- .psiRef(f)
  got <- spiDE:::.psiProfileBatch(f$Y, f$Mu, f$range)
  expect_named(got, c("psi", "at_bound"), ignore.order = TRUE)
  # compare in the space the search runs in, against optimize's tolerance there
  free <- !ref$at_bound
  expect_true(any(free))
  expect_lt(max(abs(log(got$psi[free]) - log(ref$psi[free]))), 2e-4)
})

test_that(".psiProfileBatch never finds a worse optimum than optimize()", {
  # the gate. Both maximise the same function; the bisection is the more
  # accurate search, so it must not lose.
  f <- .psiFixture()
  ref <- .psiRef(f)
  got <- spiDE:::.psiProfileBatch(f$Y, f$Mu, f$range)
  for (i in seq_len(nrow(f$Y))) {
    ll_new <- .nbll(f$Y[i, ], f$Mu[i, ], got$psi[i])
    ll_old <- .nbll(f$Y[i, ], f$Mu[i, ], ref$psi[i])
    expect_gte(ll_new, ll_old - 1e-8 * abs(ll_old))
  }
})

test_that(".psiProfileBatch reproduces the at_bound rule", {
  f <- .psiFixture()
  ref <- .psiRef(f)
  got <- spiDE:::.psiProfileBatch(f$Y, f$Mu, f$range)
  expect_identical(got$at_bound, ref$at_bound)
  # the Poisson gene is the one that should be flagged
  expect_true(got$at_bound[nrow(f$Y)])
})

test_that(".psiProfileBatch is fixed-iteration, so more steps only refine", {
  # no divergent control flow across the batch: the answer at 60 steps is the
  # answer at 50, refined, never a different root
  f <- .psiFixture()
  a <- spiDE:::.psiProfileBatch(f$Y, f$Mu, f$range, maxit = 50L)
  b <- spiDE:::.psiProfileBatch(f$Y, f$Mu, f$range, maxit = 60L)
  expect_equal(a$psi, b$psi, tolerance = 1e-8)
  expect_identical(a$at_bound, b$at_bound)
})

test_that(".psiProfileBatch agrees between the base-R and torch branches", {
  skip_if_not_installed("torch")
  f <- .psiFixture()
  base <- spiDE:::.psiProfileBatch(f$Y, f$Mu, f$range)
  tt <- function(x) torch::torch_tensor(x, dtype = torch::torch_float64())
  tor <- spiDE:::.psiProfileBatch(tt(f$Y), tt(f$Mu), f$range)
  expect_equal(as.numeric(SpaNorm::toRMatrix(tor$psi)), base$psi, tolerance = 1e-9)
  expect_identical(tor$at_bound, base$at_bound)
})

test_that(".reprofilePsi uses the same kernel as the batched engine", {
  # .reprofilePsi() runs at the END of the tau2 loop, in the default path
  # (psi.method = "profile"), and overwrites @psi. Left on optimize() it would
  # discard the bisection's more accurate answer at the last step, so folding
  # it onto the same kernel is what makes the change mean anything for
  # production rather than only for the engine's internals.
  set.seed(31)
  n <- 150; p <- 3; b <- 5
  W <- cbind(1, matrix(stats::rnorm(n * (p - 1)), n, p - 1))
  alpha <- matrix(stats::rnorm(b * p, 0, 0.2), b, p); alpha[, 1] <- 2
  mu <- exp(alpha %*% t(W))
  Y <- matrix(0L, b, n)
  for (i in seq_len(b - 1L)) Y[i, ] <- stats::rnbinom(n, mu = mu[i, ], size = 3)
  Y[b, ] <- stats::rpois(n, mu[b, ])          # at_bound: keeps its incoming psi
  psi_in <- rep(0.35, b)

  got <- spiDE:::.reprofilePsi(Y, W, alpha, psi_in)
  expect_length(got, b)
  expect_true(all(is.finite(got)))

  # the same rule as before: an at-bound gene keeps the dispersion it came with
  ref_kernel <- spiDE:::.psiProfileBatch(Y, pmax(exp(alpha %*% t(W)), spiDE:::.MU_FLOOR))
  expect_equal(got[ref_kernel$at_bound], psi_in[ref_kernel$at_bound])

  # THE assertion: one optimiser, not two. A free gene's re-profiled dispersion
  # must be the kernel's answer exactly, not a second search that merely agrees
  # with it to a few decimals. An "at least as good as optimize()" gate would
  # pass on the old implementation too -- optimize() is as good as itself -- so
  # it would not drive this change.
  free <- !ref_kernel$at_bound
  expect_true(any(free))
  expect_equal(got[free], ref_kernel$psi[free], tolerance = 1e-12)
})
