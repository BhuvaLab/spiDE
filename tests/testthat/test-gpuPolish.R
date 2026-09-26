# The polish stage on an accelerator. Mirrors test-gpuInference.R: every test
# here skips without one. The engine is SpaNorm::polishNB() (SpaNorm's
# test-polishEngine.R and test-polishNB.R test its CPU side); spiDE's CPU-side
# checks are test-polish-batch.R and test-polish-backend.R.
#
# MPS is skipped on top of that, because the polish REFUSES a single-precision
# device rather than warning and proceeding -- test-polish-backend.R pins the
# refusal itself, which needs no device at all.

# backend = "cpu" on the FIT is deliberate. These tests are about the POLISH
# stage's backend, and fitSpiDE() defaults to "auto", so on a GPU node the fit
# would go to the device too -- where SpaNorm's internal fitNBGivenPsi() dies on this
# fixture ("NAs are not allowed in subscripted assignments", H100, job
# 28552404: a nested mixed fit over six genes). That is the fit stage's
# fragility, upstream of anything here, and letting it into the fixture would
# mean these tests never reach the polish at all.
.mixedFit <- function(sigma = 20) {
  spe <- buildNiches(.toySPE(n_genes = 6), sigma = sigma)
  res <- fitSpiDE(spe, "condition", sigma = sigma, random = "intercept",
                  re.celltype = TRUE, re.maxit = 2L, verbose = FALSE,
                  backend = "cpu")
  list(spe = spe, res = res)
}

test_that("the device path isolates to the device, not to the policy", {
  # polishSpiDE(backend = "gpu") changes TWO things against backend = "cpu":
  # the device, and the factorisation policy (the device requires a shared
  # factorisation; the CPU keeps the per-gene list because it is 1.7x faster
  # there). Comparing them directly confounds the two. This compares a shared
  # factorisation with itself, on host tensors against device tensors, so the
  # only difference is placement.
  #
  # The host side reaches the shared factorisation through SpaNorm::polishNB()
  # by reporting a GPU and making the device transfer the identity, so it runs
  # on base matrices. It used to run on CPU torch tensors, which torch 0.17
  # builds by aliasing R memory without keeping the R object alive: a
  # tensor-path result on a CPU device depends on garbage-collection timing.
  skip_if_no_gpu()
  skip_if(SpaNorm::getBackendDevice() == "mps", "float64 is refused on MPS")
  d <- .mixedFit()
  f <- fits(d$res)[[1]]
  Y <- as.matrix(SummarizedExperiment::assay(d$spe, "counts"))[rownames(f@alpha), ]
  pen <- f@penalty; if (length(pen) == 1L) pen <- rep(pen, ncol(f@W))
  nested <- !is.na(f@re_group) & f@re_group == "SampleCellTypeInt"
  ct_cols <- as.character(f@covtype) == "CellType"
  polish <- function() {
    SpaNorm::polishNB(Y, f@W, f@alpha, f@psi, lambda.a = pen, absorb = nested,
                      start.cols = ct_cols, engine = "batch", backend = "gpu")
  }

  host <- testthat::with_mocked_bindings(
    polish(),
    checkGPU = function(...) TRUE,
    toGPUMatrix = function(x, ...) x,
    .package = "SpaNorm"
  )
  dev <- polish()

  # FIRST: something must actually have been polished. Without this the test
  # passes when BOTH sides polish nothing -- which is exactly what happened
  # (H100, job 28555677): psi = 0 made every likelihood NaN on the tensor path,
  # every gene was dropped as non-finite, and "dev agrees with host" was
  # perfectly true and perfectly useless.
  expect_true(all(host$polish$polished))
  expect_true(all(is.finite(host$loglik)))

  expect_equal(dev$alpha, host$alpha, tolerance = gpu_tol())
  expect_equal(dev$psi, host$psi, tolerance = gpu_tol())
  expect_identical(dev$polish$polished, host$polish$polished)
  expect_identical(dev$polish$singular, host$polish$singular)
})

test_that("polishSpiDE on the GPU matches the CPU answer end to end", {
  skip_if_no_gpu()
  skip_if(SpaNorm::getBackendDevice() == "mps", "float64 is refused on MPS")
  d <- .mixedFit()
  cpu <- polishSpiDE(d$res, d$spe, tau2 = FALSE, verbose = FALSE, backend = "cpu")
  gpu <- polishSpiDE(d$res, d$spe, tau2 = FALSE, verbose = FALSE, backend = "gpu")
  fc <- fits(cpu)[[1]]; fg <- fits(gpu)[[1]]
  # this one DOES cross the policy as well as the device, which is the honest
  # user-facing question: does asking for a GPU give the same answer?
  expect_equal(fg@alpha, fc@alpha, tolerance = 1e-5)
  expect_equal(fg@psi, fc@psi, tolerance = 1e-5)
  expect_identical(fg@polish$polished, fc@polish$polished)
  # and it must not be WORSE on the objective, the gate used throughout
  expect_true(all(fg@loglik >= fc@loglik - 1e-7 * abs(fc@loglik)))
})

test_that("the device path refuses the per-gene engine", {
  skip_if_no_gpu()
  skip_if(SpaNorm::getBackendDevice() == "mps", "float64 is refused on MPS")
  d <- .mixedFit()
  expect_error(
    polishSpiDE(d$res, d$spe, tau2 = FALSE, verbose = FALSE,
                backend = "gpu", engine = "gene"),
    "engine")
})
