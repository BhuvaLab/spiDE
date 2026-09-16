# The polish stage on an accelerator. Mirrors test-gpuInference.R: every test
# here skips without one, and the CPU-side logic is tested unconditionally in
# test-polish-batch.R and test-polish-backend.R.
#
# MPS is skipped on top of that, because the polish REFUSES a single-precision
# device rather than warning and proceeding -- test-polish-backend.R pins the
# refusal itself, which needs no device at all.

.mixedFit <- function(sigma = 20) {
  spe <- buildNiches(.toySPE(n_genes = 6), sigma = sigma)
  res <- fitSpiDE(spe, "condition", sigma = sigma, random = "intercept",
                  re.celltype = TRUE, re.maxit = 2L, verbose = FALSE)
  list(spe = spe, res = res)
}

test_that("the device path isolates to the device, not to the policy", {
  # polishSpiDE(backend = "gpu") changes TWO things against backend = "cpu":
  # the device, and the factorisation policy (the device requires a shared
  # factorisation; the CPU keeps the per-gene list because it is 1.7x faster
  # there). Comparing them directly confounds the two. This compares a shared
  # factorisation with itself, on host tensors against device tensors, so the
  # only difference is placement.
  skip_if_no_gpu()
  skip_if(SpaNorm::getBackendDevice() == "mps", "float64 is refused on MPS")
  d <- .mixedFit()
  f <- fits(d$res)[[1]]
  Y <- as.matrix(SummarizedExperiment::assay(d$spe, "counts"))
  pen <- f@penalty; if (length(pen) == 1L) pen <- rep(pen, ncol(f@W))
  nested <- !is.na(f@re_group) & f@re_group == "SampleCellTypeInt"
  ct_cols <- as.character(f@covtype) == "CellType"
  solver <- spiDE:::.newtonSolver(f@W, pen, nested)

  tt <- function(x) torch::torch_tensor(x, dtype = torch::torch_float64())
  host <- spiDE:::.polishBatch(tt(Y), tt(f@W), tt(f@alpha), f@psi, pen, solver,
                               ct_cols = ct_cols, shared.factor = TRUE,
                               nested = nested)
  dev <- spiDE:::.polishBatch(
    SpaNorm::toGPUMatrix(Y, backend = "gpu"),
    SpaNorm::toGPUMatrix(f@W, backend = "gpu"),
    SpaNorm::toGPUMatrix(f@alpha, backend = "gpu"),
    f@psi, pen, solver, ct_cols = ct_cols, shared.factor = TRUE,
    nested = nested)

  expect_equal(dev$alpha, host$alpha, tolerance = gpu_tol())
  expect_equal(dev$psi, host$psi, tolerance = gpu_tol())
  expect_identical(dev$polished, host$polished)
  expect_identical(dev$singular, host$singular)
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
