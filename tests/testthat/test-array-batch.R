# The handful of array operations newton() does on (genes x cells) and
# (genes x columns) blocks, behind one signature each so the loop can be
# written once and run on a matrix or a tensor.
#
# Every base-R branch is plain indexing, deliberately: converting newton() to
# these must not move the CPU path at all, and test-polish-batch.R's strict
# parity tests are what proves it.

test_that(".setRows writes rows on either backend", {
  set.seed(2)
  X <- matrix(stats::rnorm(20), 5, 4)
  V <- matrix(1:8 / 10, 2, 4)
  ref <- X; ref[c(2L, 4L), ] <- V
  expect_equal(spiDE:::.setRows(X, c(2L, 4L), V), ref)

  skip_if_not_installed("torch")
  tt <- function(x) torch::torch_tensor(x, dtype = torch::torch_float64())
  got <- spiDE:::.setRows(tt(X), c(2L, 4L), tt(V))
  expect_equal(as.matrix(SpaNorm::toRMatrix(got)), ref, tolerance = 1e-12)
})

test_that(".mulRows scales each row by its own scalar", {
  set.seed(3)
  M <- matrix(stats::rnorm(15), 3, 5)
  v <- c(2, -1, 0.5)
  ref <- v * M                     # base R recycles down columns
  expect_equal(spiDE:::.mulRows(v, M), ref)

  skip_if_not_installed("torch")
  tt <- function(x) torch::torch_tensor(x, dtype = torch::torch_float64())
  got <- spiDE:::.mulRows(v, tt(M))
  expect_equal(as.matrix(SpaNorm::toRMatrix(got)), ref, tolerance = 1e-12)
})

test_that(".scaleCols scales each column by its own scalar", {
  set.seed(4)
  A <- matrix(stats::rnorm(12), 3, 4)
  s <- c(1, 0, 2, -0.5)
  ref <- sweep(A, 2L, s, `*`)
  expect_equal(spiDE:::.scaleCols(A, s), ref)

  skip_if_not_installed("torch")
  tt <- function(x) torch::torch_tensor(x, dtype = torch::torch_float64())
  got <- spiDE:::.scaleCols(tt(A), s)
  expect_equal(as.matrix(SpaNorm::toRMatrix(got)), ref, tolerance = 1e-12)
})

test_that(".asHost brings a length-b vector back for the bookkeeping", {
  # the control flow -- which genes are active, which accepted, how many
  # halvings -- stays on the host; only the (genes x cells) work is on device
  v <- c(1.5, -2, 0)
  expect_identical(spiDE:::.asHost(v), v)

  skip_if_not_installed("torch")
  got <- spiDE:::.asHost(torch::torch_tensor(v, dtype = torch::torch_float64()))
  expect_equal(got, v, tolerance = 1e-12)
  expect_true(is.numeric(got))
})

test_that(".matmulB multiplies on either backend", {
  set.seed(6)
  R <- matrix(stats::rnorm(12), 3, 4)
  W <- matrix(stats::rnorm(8), 4, 2)
  expect_equal(spiDE:::.matmulB(R, W), R %*% W)

  skip_if_not_installed("torch")
  tt <- function(x) torch::torch_tensor(x, dtype = torch::torch_float64())
  got <- spiDE:::.matmulB(tt(R), tt(W))
  expect_equal(as.matrix(SpaNorm::toRMatrix(got)), R %*% W, tolerance = 1e-12)
})
