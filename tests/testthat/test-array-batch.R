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

test_that(".rowsFinite flags rows carrying NA, NaN or Inf", {
  X <- matrix(1:12 / 2, 3, 4)
  X[2, 3] <- NA; X[3, 1] <- Inf
  expect_identical(spiDE:::.rowsFinite(X), c(TRUE, FALSE, FALSE))

  skip_if_not_installed("torch")
  got <- spiDE:::.rowsFinite(torch::torch_tensor(X, dtype = torch::torch_float64()))
  expect_identical(got, c(TRUE, FALSE, FALSE))
})

test_that(".rowMeansB returns host numbers on either backend", {
  set.seed(8)
  X <- matrix(stats::rnorm(20), 4, 5)
  expect_equal(spiDE:::.rowMeansB(X), rowMeans(X))

  skip_if_not_installed("torch")
  got <- spiDE:::.rowMeansB(torch::torch_tensor(X, dtype = torch::torch_float64()))
  expect_equal(got, rowMeans(X), tolerance = 1e-12)
  expect_true(is.numeric(got))
})

test_that(".colsOf subsets columns on either backend", {
  set.seed(9)
  X <- matrix(stats::rnorm(20), 4, 5)
  expect_equal(spiDE:::.colsOf(X, c(2L, 5L)), X[, c(2L, 5L), drop = FALSE])

  skip_if_not_installed("torch")
  got <- spiDE:::.colsOf(torch::torch_tensor(X, dtype = torch::torch_float64()),
                         c(2L, 5L))
  expect_equal(as.matrix(SpaNorm::toRMatrix(got)), X[, c(2L, 5L)], tolerance = 1e-12)
})

test_that(".maskedRowMin takes the minimum over the masked cells only", {
  # degenerate() asks for the smallest linear predictor among cells with a
  # POSITIVE count; a gene with no positive count has no such cell and must
  # come back Inf rather than empty or NA
  X <- rbind(c(1, -5, 3), c(2, 2, 2), c(-9, 0, 4))
  mask <- rbind(c(TRUE, FALSE, TRUE), c(FALSE, FALSE, FALSE), c(TRUE, TRUE, TRUE))
  expect_equal(spiDE:::.maskedRowMin(X, mask), c(1, Inf, -9))

  skip_if_not_installed("torch")
  tt <- function(x) torch::torch_tensor(x, dtype = torch::torch_float64())
  got <- spiDE:::.maskedRowMin(tt(X), tt(mask * 1)$to(dtype = torch::torch_bool()))
  expect_equal(got, c(1, Inf, -9), tolerance = 1e-12)
})

test_that(".asLike moves a host matrix onto the reference's backend", {
  M <- matrix(1:6 / 3, 2, 3)
  expect_identical(spiDE:::.asLike(M, matrix(0, 1, 1)), M)

  skip_if_not_installed("torch")
  ref <- torch::torch_tensor(matrix(0, 1, 1), dtype = torch::torch_float64())
  got <- spiDE:::.asLike(M, ref)
  expect_true(SpaNorm::is_torch_tensor(got))
  expect_equal(as.matrix(SpaNorm::toRMatrix(got)), M, tolerance = 1e-12)
})

test_that(".rowMaxB and .asHostMat close the boundary back to the host", {
  set.seed(12)
  X <- matrix(stats::rnorm(20), 4, 5)
  expect_equal(spiDE:::.rowMaxB(X), apply(X, 1L, max))
  expect_identical(spiDE:::.asHostMat(X), X)

  skip_if_not_installed("torch")
  Xt <- torch::torch_tensor(X, dtype = torch::torch_float64())
  expect_equal(spiDE:::.rowMaxB(Xt), apply(X, 1L, max), tolerance = 1e-12)
  got <- spiDE:::.asHostMat(Xt)
  expect_true(is.matrix(got))
  expect_equal(got, X, tolerance = 1e-12)
})
