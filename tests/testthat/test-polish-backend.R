# The polish stage's backend argument, and the one thing it refuses.
#
# A 447-column penalised Newton solve in single precision is not defensible, so
# a float32 device is an ERROR rather than a warning-and-proceed.
# SpaNorm::getBackendDtype() is float64 on CUDA but float32 on MPS, which is
# the case this exists for: the failure would otherwise be silent and look like
# a convergence problem.

test_that(".requireFloat64 refuses a single-precision device", {
  expect_error(spiDE:::.requireFloat64("Float"), "single precision")
  expect_error(spiDE:::.requireFloat64("float32"), "single precision")
  expect_silent(spiDE:::.requireFloat64("Double"))
  expect_silent(spiDE:::.requireFloat64("float64"))
  # the message has to say what to do, not only what is wrong
  err <- tryCatch(spiDE:::.requireFloat64("Float"), error = function(e) conditionMessage(e))
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
