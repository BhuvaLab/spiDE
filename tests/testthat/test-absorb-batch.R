# The nested (sample x cell type) Schur absorption, batched and
# backend-agnostic. .newtonSolver()'s per-gene `parts()` is the oracle: it is
# the implementation every converged fit and the CPU inference path already
# use, and test-polish.R pins it against the dense inverse.
#
# The point of the batched form is the GPU inference path, which today skips
# the absorption entirely and builds a dense p x p gram -- 1,107 columns where
# 398 would do on the cohort's design, 7.7x the flops. Everything here is
# testable on CPU torch tensors; only device placement needs an accelerator.

# a design with a genuine nested indicator block: px dense columns, then G
# 0/1 columns that partition the cells
.absorbFixture <- function(n = 60, px = 4, G = 5, b = 3, seed = 7,
                           empty_group = FALSE) {
  set.seed(seed)
  X <- cbind(1, matrix(stats::rnorm(n * (px - 1)), n, px - 1))
  grp <- if (empty_group) {
    # group 2 gets no cells: rowsum() drops it, and the caller indexes cvec
    # positionally, so a dropped group silently misaligns every group after it
    sample(setdiff(seq_len(G), 2L), n, replace = TRUE)
  } else {
    c(seq_len(G), sample(seq_len(G), n - G, replace = TRUE))
  }
  Z <- matrix(0, n, G)
  Z[cbind(seq_len(n), grp)] <- 1
  W <- cbind(X, Z)
  list(W = W, X = X, grp = grp, G = G, px = px,
       nested = c(rep(FALSE, px), rep(TRUE, G)),
       pen = c(rep(0, px), rep(0.3, G)),
       wt = matrix(stats::runif(b * n, 0.2, 2), b, n))
}

.refS <- function(f) {
  sol <- spiDE:::.newtonSolver(f$W, f$pen, f$nested)
  out <- array(0, c(nrow(f$wt), f$px, f$px))
  for (g in seq_len(nrow(f$wt))) out[g, , ] <- sol$factor(f$wt[g, ])$S
  out
}

# S = A - B C^-1 B' written out, independent of .newtonSolver(). Needed for the
# empty-group case: .newtonSolver()'s rowsum() DROPS a group with no cells, so
# cvec comes back short and pen_z recycles against it. That is unreachable from
# the package's own designs -- .buildRandomEffects() builds the nested block
# with interaction(drop = TRUE), so every column holds at least one cell -- but
# .absorbBatch() gets the robustness free from .segmentSum() and should not
# quietly lose it.
.refS_direct <- function(f) {
  Z <- f$W[, f$nested, drop = FALSE]
  out <- array(0, c(nrow(f$wt), f$px, f$px))
  for (g in seq_len(nrow(f$wt))) {
    w <- f$wt[g, ]
    A <- crossprod(f$X * w, f$X) + diag(f$pen[!f$nested], f$px)
    Cd <- colSums(Z * w) + f$pen[f$nested]
    B <- crossprod(f$X * w, Z)
    out[g, , ] <- A - B %*% (t(B) / Cd)
  }
  out
}

test_that(".absorbBatch reproduces .newtonSolver()'s Schur complement per gene", {
  f <- .absorbFixture()
  got <- spiDE:::.absorbBatch(f$W, f$pen, f$nested, f$wt)
  expect_equal(dim(got), c(nrow(f$wt), f$px, f$px))
  expect_equal(got, .refS(f), tolerance = 1e-10)
  # and against the definition, so the two references corroborate each other
  expect_equal(got, .refS_direct(f), tolerance = 1e-10)
})

test_that(".absorbBatch keeps a group that holds no cells", {
  # compared against the written-out definition, NOT .newtonSolver(): the
  # oracle misaligns here, and that is a property of rowsum() rather than of
  # the absorption (see .refS_direct above)
  f <- .absorbFixture(empty_group = TRUE)
  expect_equal(spiDE:::.absorbBatch(f$W, f$pen, f$nested, f$wt),
               .refS_direct(f), tolerance = 1e-10)
})

test_that(".absorbBatch is invariant to the cell tile", {
  # the tile bounds the (batch, tile, px) intermediate; it is performance, not
  # semantics, so every tiling must give the same stack
  f <- .absorbFixture()
  whole <- spiDE:::.absorbBatch(f$W, f$pen, f$nested, f$wt)
  for (tile in c(1L, 7L, 59L, 60L, 1000L)) {
    expect_equal(spiDE:::.absorbBatch(f$W, f$pen, f$nested, f$wt,
                                      cell.tile = tile),
                 whole, tolerance = 1e-12,
                 info = sprintf("cell.tile = %d", tile))
  }
})

test_that(".absorbBatch agrees between the base-R and torch branches", {
  skip_if_not_installed("torch")
  f <- .absorbFixture()
  base_S <- spiDE:::.absorbBatch(f$W, f$pen, f$nested, f$wt)
  Wt <- torch::torch_tensor(f$W, dtype = torch::torch_float64())
  wtt <- torch::torch_tensor(f$wt, dtype = torch::torch_float64())
  tor_S <- spiDE:::.absorbBatch(Wt, f$pen, f$nested, wtt)
  expect_true(SpaNorm::is_torch_tensor(tor_S))
  expect_equal(as.array(tor_S), base_S, tolerance = 1e-10)
})

test_that(".absorbBatch refuses columns that are not a partition", {
  f <- .absorbFixture()
  # cell 1 is in group 1 by construction, so give it group 2 as well
  f$W[1, f$px + 2] <- 1
  expect_error(spiDE:::.absorbBatch(f$W, f$pen, f$nested, f$wt),
               "partition")
})
