# The polish stage's backend argument, and the one thing it refuses.
#
# A 447-column penalised Newton solve in single precision is not defensible, so
# a float32 device is an ERROR rather than a warning-and-proceed.
# SpaNorm::getBackendDtype() is float64 on CUDA but float32 on MPS, which is
# the case this exists for: the failure would otherwise be silent and look like
# a convergence problem.

test_that("the device path refuses a single-precision device", {
  # The refusal lives in SpaNorm::polishNB() since the polish moved there.
  # Reach it with no device at all: report a GPU, make the transfer the
  # identity and set the backend dtype, so only the precision differs.
  set.seed(3)
  n <- 60
  W <- cbind(1, scale(rnorm(n)))
  Y <- matrix(rnbinom(2 * n, mu = 3, size = 2), 2)
  A0 <- matrix(c(1, 0), 2, 2, byrow = TRUE)
  run <- function(dtype) {
    testthat::with_mocked_bindings(
      SpaNorm::polishNB(Y, W, A0, c(0.5, 0.5), backend = "gpu"),
      checkGPU = function(...) TRUE,
      toGPUMatrix = function(x, ...) x,
      getBackendDtype = function(...) dtype,
      .package = "SpaNorm"
    )
  }
  expect_error(run("Float"), "single precision")
  expect_error(run("float32"), "single precision")
  expect_no_error(run("Double"))
  expect_no_error(run("float64"))
  # the message has to say what to do, not only what is wrong
  err <- tryCatch(run("Float"), error = function(e) conditionMessage(e))
  expect_match(err, "backend")
})

test_that("asking for a GPU without one falls back to the CPU path exactly", {
  # backend = "gpu" on a machine with no accelerator must give the CPU answer,
  # not an error: the device path is an accelerator, never a requirement
  skip_if(SpaNorm::checkGPU(), "an accelerator is present; this tests the fallback")
  spe <- buildNiches(.toySPE(n_genes = 5), sigma = 20)
  res <- fitSpiDE(spe, "condition", sigma = 20, random = "intercept",
                  re.maxit = 2L, verbose = FALSE)
  a <- polishSpiDE(res, spe, tau2 = FALSE, verbose = FALSE, backend = "cpu")
  b <- polishSpiDE(res, spe, tau2 = FALSE, verbose = FALSE, backend = "gpu")
  fa <- fits(a)[[1]]; fb <- fits(b)[[1]]
  expect_equal(fb@alpha, fa@alpha, tolerance = 1e-12)
  expect_equal(fb@psi, fa@psi, tolerance = 1e-12)
})
