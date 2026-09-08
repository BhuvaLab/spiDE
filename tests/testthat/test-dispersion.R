# The dispersion rule of the convergence step and the standard-error scale of
# the inference, measured on the synthetic benchmark and the real cohort
# (research: fdr-ordering/FINDINGS.md, 2026-09-08): the psi estimate is
# irrelevant to the null and the moderated one is a third cheaper; the
# quasi-likelihood scale is the only one that holds the nominal level at
# every sample size. Both are the defaults. The QL machinery itself lives in
# SpaNorm (qlDispersion, nbUnitDeviance, nbDevianceMoments) with its oracle
# tests against edgeR; spiDE only wires it.

toy_fit <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      spe <- spiDE:::.toySPE()
      spe <- buildNiches(spe, sigma = 30, verbose = FALSE)
      cache <<- list(spe = spe, fit = fitSpiDE(spe, "condition", sigma = 30, verbose = FALSE))
    }
    cache
  }
})

test_that("the default convergence step keeps fitNB's moderated dispersion at the converged mean", {
  tf <- toy_fit()
  fm <- fits(tf$fit)[[1]]
  expect_equal(unname(fm@psi), unname(fm@polish$psi_fitnb))
  expect_true(all(fm@polish$polished | fm@polish$iterations == 0L))
  fp <- fits(fitSpiDE(tf$spe, "condition", sigma = 30, polish.psi = "profile", verbose = FALSE))[[1]]
  expect_false(isTRUE(all.equal(unname(fp@psi), unname(fm@psi))))
})

test_that("the default inference scales the standard errors by the quasi-likelihood dispersion", {
  tf <- toy_fit()
  rd <- fits(testSpiDE(tf$fit, spe = tf$spe, fdr = 1))[[1]]
  rq <- fits(testSpiDE(tf$fit, spe = tf$spe, fdr = 1, dispersion = "ql"))[[1]]
  rp <- fits(testSpiDE(tf$fit, spe = tf$spe, fdr = 1, dispersion = "pearson"))[[1]]
  expect_equal(rd@t_stat, rq@t_stat)
  expect_true(all(is.finite(rq@t_stat)))
  expect_false(isTRUE(all.equal(rq@t_stat, rp@t_stat)))
  ratio <- rq@t_stat / rp@t_stat
  expect_true(all(abs(log(ratio[is.finite(ratio) & rp@t_stat != 0])) < log(3)))
  # the wrapper is an S4 generic: its defaults live in the method definition
  m <- paste(deparse(methods::selectMethod("spiDE", "SpatialExperiment")), collapse = " ")
  expect_match(m, 'dispersion = c\\("ql", "pearson"\\)')
  expect_match(m, 'polish.psi = c\\("moderated", "profile"\\)')
})

test_that("spiDE carries no private copy of the QL machinery", {
  ns <- asNamespace("spiDE")
  expect_false(exists(".qlDispersion", envir = ns, inherits = FALSE))
  expect_false(exists(".nbDevianceMoments", envir = ns, inherits = FALSE))
  expect_false(exists(".nbUnitDeviance", envir = ns, inherits = FALSE))
})

test_that("a fixed-effects, unconverged fit keeps its legacy scale under the default, and says so", {
  tf <- toy_fit()
  f0 <- fitSpiDE(tf$spe, "condition", sigma = 30, random = "none", converge = FALSE, verbose = FALSE)
  expect_message(r0 <- testSpiDE(f0, spe = tf$spe, fdr = 1), "legacy")
  rl <- suppressMessages(testSpiDE(f0, spe = tf$spe, fdr = 1, dispersion = "pearson"))
  expect_equal(fits(r0)[[1]]@t_stat, fits(rl)[[1]]@t_stat)
})

test_that("the QL scale gives the same statistics on the GPU backend", {
  skip_if_not(SpaNorm::checkGPU(), "no accelerator")
  tf <- toy_fit()
  rc <- fits(testSpiDE(tf$fit, spe = tf$spe, fdr = 1, backend = "cpu"))[[1]]
  rg <- fits(testSpiDE(tf$fit, spe = tf$spe, fdr = 1, backend = "gpu"))[[1]]
  expect_equal(rg@t_stat, rc@t_stat, tolerance = 1e-6)
})
