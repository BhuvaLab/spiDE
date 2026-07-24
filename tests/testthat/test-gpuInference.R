# GPU-accelerated inference: parity with the CPU path, memory-bounded blocking,
# the BiocParallel/GPU contention guard, and fit-stage GPU parity. GPU-path
# tests skip when no accelerator is available (mirroring SpaNorm's own
# test-gpuFunctions.R); the CPU-only helpers are tested unconditionally in
# test-inference.R.

# a fitted SpiDEFit + counts on the toy data (fixed effects)
.toyFit <- function(sigma = 20) {
  spe <- buildNiches(.toySPE(), sigma = sigma)
  spe <- computeSizeFactors(spe, count = "nCount", area = "Area")
  res <- fitSpiDE(spe, condition = "condition", sigma = sigma,
                  covariates = c("Age", "LS"), verbose = FALSE)
  list(fit = fits(res)[[1]],
       Y = as.matrix(SummarizedExperiment::assay(spe, "counts")))
}

test_that(".inferenceBlockSize returns NULL on the CPU path", {
  expect_null(spiDE:::.inferenceBlockSize(1000, 5000, 10, "cpu"))
})

test_that(".inferenceBlockSize splits under a tiny budget, single block under a large one", {
  skip_if_no_gpu()
  SpaNorm::resetGPUCache()
  on.exit(SpaNorm::resetGPUCache(), add = TRUE)

  small <- spiDE:::.inferenceBlockSize(1000, 5000, 10, "gpu",
                                       gpu.mem.budget = 5e6)
  expect_true(is.numeric(small) && small >= 1L && small < 1000L)
  expect_null(spiDE:::.inferenceBlockSize(1000, 5000, 10, "gpu",
                                          gpu.mem.budget = 1e12))
})

test_that("GPU inference matches the CPU path (fixed effects, cauchy)", {
  skip_if_no_gpu()
  tf <- .toyFit()
  cpu <- spiDE:::.blockedInference(tf$fit, tf$Y, combine = "cauchy",
                                   backend = "cpu")
  gpu <- spiDE:::.blockedInference(tf$fit, tf$Y, combine = "cauchy",
                                   backend = "gpu")
  expect_equal(gpu@t_stat, cpu@t_stat, tolerance = gpu_tol())
  expect_equal(gpu@se, cpu@se, tolerance = gpu_tol())
  expect_equal(gpu@p.combined.pos, cpu@p.combined.pos, tolerance = gpu_tol())
  expect_equal(gpu@p.combined.neg, cpu@p.combined.neg, tolerance = gpu_tol())
})

test_that("GPU inference matches the CPU path (fixed effects, brown NB-math)", {
  skip_if_no_gpu()
  skip_if_not_installed("poolr")
  tf <- .toyFit()
  cpu <- spiDE:::.blockedInference(tf$fit, tf$Y, combine = "brown",
                                   backend = "cpu")
  gpu <- spiDE:::.blockedInference(tf$fit, tf$Y, combine = "brown",
                                   backend = "gpu")
  expect_equal(gpu@t_stat, cpu@t_stat, tolerance = gpu_tol())
  expect_equal(gpu@p.combined.pos, cpu@p.combined.pos, tolerance = gpu_tol())
})

test_that("GPU inference matches the CPU path (mixed effects)", {
  skip_if_no_gpu()
  spe <- buildNiches(.toySPE(), sigma = 20)
  spe <- computeSizeFactors(spe, count = "nCount", area = "Area")
  set.seed(1)
  res <- fitSpiDE(spe, condition = "condition", sigma = 20,
                  random = "intercept", verbose = FALSE)
  fm <- fits(res)[[1]]
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))

  cpu <- spiDE:::.blockedInference(fm, Y, combine = "cauchy", backend = "cpu")
  gpu <- spiDE:::.blockedInference(fm, Y, combine = "cauchy", backend = "gpu")
  expect_equal(gpu@t_stat, cpu@t_stat, tolerance = gpu_tol())
  expect_equal(gpu@se, cpu@se, tolerance = gpu_tol())
  expect_equal(gpu@p.combined.pos, cpu@p.combined.pos, tolerance = gpu_tol())
})

test_that("GPU inference is invariant to a memory-bounded block size", {
  skip_if_no_gpu()
  tf <- .toyFit()
  single <- spiDE:::.blockedInference(tf$fit, tf$Y, combine = "cauchy",
                                      backend = "gpu")
  # a tiny budget forces multiple blocks; results must be unchanged
  blocked <- spiDE:::.blockedInference(tf$fit, tf$Y, combine = "cauchy",
                                       backend = "gpu", gpu.mem.budget = 1e5)
  expect_equal(blocked@t_stat, single@t_stat, tolerance = gpu_tol())
  expect_equal(blocked@p.combined.pos, single@p.combined.pos,
               tolerance = gpu_tol())
})

test_that("GPU backend forces serial BPPARAM (with a warning) and still matches", {
  skip_if_no_gpu()
  skip_on_os("windows")
  tf <- .toyFit()
  serial <- spiDE:::.blockedInference(tf$fit, tf$Y, combine = "cauchy",
                                      backend = "gpu")
  expect_warning(
    mc <- spiDE:::.blockedInference(tf$fit, tf$Y, combine = "cauchy",
                                    backend = "gpu",
                                    BPPARAM = BiocParallel::MulticoreParam(2)),
    "SerialParam"
  )
  expect_equal(mc@p.combined.pos, serial@p.combined.pos, tolerance = gpu_tol())
})

test_that("fitSpiDE GPU backend matches the CPU fit on the toy data", {
  skip_if_no_gpu()
  spe <- buildNiches(.toySPE(), sigma = 20)
  spe <- computeSizeFactors(spe, count = "nCount", area = "Area")
  cpu <- fitSpiDE(spe, condition = "condition", sigma = 20,
                  covariates = c("Age", "LS"), backend = "cpu",
                  verbose = FALSE)
  gpu <- fitSpiDE(spe, condition = "condition", sigma = 20,
                  covariates = c("Age", "LS"), backend = "gpu",
                  verbose = FALSE)
  expect_equal(fits(gpu)[[1]]@alpha, fits(cpu)[[1]]@alpha, tolerance = gpu_tol())
  expect_equal(fits(gpu)[[1]]@psi, fits(cpu)[[1]]@psi, tolerance = gpu_tol())
})
