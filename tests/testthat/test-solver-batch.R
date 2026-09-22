# A batched factor/solve for the nested design, which is what Phase 2e needs to
# put .polishBatch()'s Newton on a device.
#
# newton() currently keeps a LIST of per-gene factorisations and calls
# solver$solve(fac[[k]], S[k, ]) one gene at a time. A list of R objects cannot
# go to a device and a per-gene solve is a kernel launch per gene per
# iteration, so the state has to become one batched object.
#
# .newtonSolver() is the oracle throughout: same Schur absorption, same
# right-hand side, same back-substitution, one gene at a time.

.solverFixture <- function(n = 80, px = 4, G = 5, b = 4, seed = 3,
                           nested = TRUE) {
  set.seed(seed)
  X <- cbind(1, matrix(stats::rnorm(n * (px - 1)), n, px - 1))
  if (!nested) {
    return(list(W = X, px = px, G = 0L,
                nested = rep(FALSE, px), pen = rep(0, px),
                wt = matrix(stats::runif(b * n, 0.2, 2), b, n)))
  }
  grp <- c(seq_len(G), sample(seq_len(G), n - G, replace = TRUE))
  Z <- matrix(0, n, G)
  Z[cbind(seq_len(n), grp)] <- 1
  # pen_x = 0 is the production shape: the ridge is on the random-effect
  # columns only, which is also what makes a singular dense block reachable
  list(W = cbind(X, Z), px = px, G = G,
       nested = c(rep(FALSE, px), rep(TRUE, G)),
       pen = c(rep(0, px), rep(0.3, G)),
       wt = matrix(stats::runif(b * n, 0.2, 2), b, n))
}

# a per-gene score to solve against
.scores <- function(f) {
  set.seed(99)
  matrix(stats::rnorm(nrow(f$wt) * ncol(f$W)), nrow(f$wt), ncol(f$W))
}

test_that(".absorbBatch can return the parts, not only the Schur complement", {
  # the solve needs B and cvec as well as S, and they must be the same B and
  # cvec .newtonSolver() built
  f <- .solverFixture()
  got <- spiDE:::.absorbBatch(f$W, f$pen, f$nested, f$wt, parts = TRUE)
  expect_named(got, c("S", "B", "cvec", "xi", "zi"), ignore.order = TRUE)
  expect_equal(dim(got$S), c(nrow(f$wt), f$px, f$px))
  expect_equal(dim(got$B), c(nrow(f$wt), f$G, f$px))
  expect_equal(dim(got$cvec), c(nrow(f$wt), f$G))

  sol <- spiDE:::.newtonSolver(f$W, f$pen, f$nested)
  for (g in seq_len(nrow(f$wt))) {
    ref <- sol$factor(f$wt[g, ])
    expect_equal(got$S[g, , ], ref$S, tolerance = 1e-10)
    # .newtonSolver()'s B is px x G; the batched stack carries its transpose.
    # unname(): rowsum() puts the group levels on B's columns and cvec's names,
    # which a (batch, G, px) stack cannot carry and nothing downstream reads --
    # both index positionally.
    expect_equal(t(got$B[g, , ]), unname(ref$B), tolerance = 1e-10)
    expect_equal(got$cvec[g, ], unname(ref$cvec), tolerance = 1e-10)
  }
  # and the default return is unchanged
  expect_equal(spiDE:::.absorbBatch(f$W, f$pen, f$nested, f$wt), got$S)
})

test_that(".newtonSolverBatch's step matches .newtonSolver's, gene by gene", {
  f <- .solverFixture()
  Sc <- .scores(f)
  sb <- spiDE:::.newtonSolverBatch(f$W, f$pen, f$nested)
  st <- sb$factor(f$wt)
  D <- sb$solve(st, Sc)
  expect_equal(dim(D), dim(Sc))
  expect_true(all(st$ok))

  sol <- spiDE:::.newtonSolver(f$W, f$pen, f$nested)
  for (g in seq_len(nrow(f$wt))) {
    expect_equal(D[g, ], sol$solve(f$wt[g, ], Sc[g, ]), tolerance = 1e-9)
  }
})

test_that(".newtonSolverBatch works with no nested block at all", {
  f <- .solverFixture(nested = FALSE)
  Sc <- .scores(f)
  sb <- spiDE:::.newtonSolverBatch(f$W, f$pen, f$nested)
  D <- sb$solve(sb$factor(f$wt), Sc)
  sol <- spiDE:::.newtonSolver(f$W, f$pen, f$nested)
  for (g in seq_len(nrow(f$wt))) {
    expect_equal(D[g, ], sol$solve(f$wt[g, ], Sc[g, ]), tolerance = 1e-9)
  }
})

test_that(".newtonSolverBatch's xcov matches .newtonSolver's", {
  f <- .solverFixture()
  sb <- spiDE:::.newtonSolverBatch(f$W, f$pen, f$nested)
  V <- sb$xcov(sb$factor(f$wt))
  sol <- spiDE:::.newtonSolver(f$W, f$pen, f$nested)
  for (g in seq_len(nrow(f$wt))) {
    expect_equal(V[g, , ], sol$xcov(f$wt[g, ]), tolerance = 1e-9)
  }
})

test_that("a singular gene does not poison its batch", {
  # THE trap this whole batched path has to avoid, and the one
  # .waldCauchyBlock() fell into (inference.R: one singular gene kills the
  # sub-batch's Cholesky and the only recourse is telling the user to shrink
  # cov.batch). With pen_x = 0, a gene carrying no weight has a zero dense
  # block: singular, and its neighbours must be unaffected.
  f <- .solverFixture()
  f$wt[2, ] <- 0
  Sc <- .scores(f)
  sb <- spiDE:::.newtonSolverBatch(f$W, f$pen, f$nested)
  st <- sb$factor(f$wt)
  expect_false(st$ok[2])
  expect_true(all(st$ok[-2]))

  D <- sb$solve(st, Sc)
  expect_true(all(is.na(D[2, ])))
  sol <- spiDE:::.newtonSolver(f$W, f$pen, f$nested)
  for (g in seq_len(nrow(f$wt))[-2]) {
    expect_equal(D[g, ], sol$solve(f$wt[g, ], Sc[g, ]), tolerance = 1e-9)
  }
})

test_that(".newtonSolverBatch agrees between the base-R and torch branches", {
  skip_if_not_installed("torch")
  f <- .solverFixture()
  Sc <- .scores(f)
  base_D <- spiDE:::.newtonSolverBatch(f$W, f$pen, f$nested) |>
    (\(sb) sb$solve(sb$factor(f$wt), Sc))()

  Wt <- torch::torch_tensor(f$W, dtype = torch::torch_float64())
  wtt <- torch::torch_tensor(f$wt, dtype = torch::torch_float64())
  Sct <- torch::torch_tensor(Sc, dtype = torch::torch_float64())
  sbt <- spiDE:::.newtonSolverBatch(Wt, f$pen, f$nested)
  tor_D <- sbt$solve(sbt$factor(wtt), Sct)
  expect_true(SpaNorm::is_torch_tensor(tor_D))
  expect_equal(as.array(tor_D), base_D, tolerance = 1e-9)
})
