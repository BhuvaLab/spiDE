# The polish is a stage of the pipeline (fit -> polish -> test -> gsea), not
# an option of the fit: it converges each gene, sets its dispersion, and for a
# mixed fit re-estimates the variance components from the converged fit. The
# shared fit's Schall loop reads its own unconverged coefficients, which on the
# clustered fixture reports a between-sample variance of 10 against a planted
# 0.49 (research fdr-ordering/FINDINGS.md, 2026-09-08).

clustered <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      spe <- buildNiches(spiDE:::.toyClustered(n_samples = 8, n_per = 40, n_genes = 12,
                                                sd_patient = 0.7), sigma = 30, verbose = FALSE)
      cache <<- list(spe = spe,
                     fit = fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                                    df.method = "satterthwaite", verbose = FALSE))
    }
    cache
  }
})

test_that("fitSpiDE() fits the shared model and nothing more", {
  fm <- methods::selectMethod("fitSpiDE", "SpatialExperiment")
  expect_false(any(c("converge", "converge.maxit", "converge.tol", "polish.psi") %in%
                     names(formals(fm))))
  f <- fits(clustered()$fit)[[1]]
  expect_null(f@polish)
})

test_that("polishSpiDE() re-estimates the variance components from the converged fit", {
  cl <- clustered()
  pol <- polishSpiDE(cl$fit, cl$spe, verbose = FALSE)
  f <- fits(pol)[[1]]
  after <- f@tau2[["SampleInt"]]
  # planted sd_patient 0.7 -> 0.49. On this small fixture the shared fit's
  # loop is already close; on the full-size one it reports 10, and the long
  # test test-mixed-numerics.R holds the polished value to the window there.
  expect_gt(after, 0.25)
  expect_lt(after, 0.9)
  expect_equal(unname(f@penalty[f@re_group %in% "SampleInt"]), rep(1 / after, sum(f@re_group %in% "SampleInt")))
  # the reference df follows the new components
  expect_true(is.numeric(f@df) && !is.null(names(f@df)))
  expect_false(isTRUE(all.equal(f@df, fits(cl$fit)[[1]]@df)))
  # and the stage can be told to leave them alone
  keep <- fits(polishSpiDE(cl$fit, cl$spe, tau2 = FALSE, verbose = FALSE))[[1]]
  expect_equal(keep@tau2, fits(cl$fit)[[1]]@tau2)
})

test_that("the polish's dispersion rule defaults to the profile value", {
  cl <- clustered()
  pr <- fits(polishSpiDE(cl$fit, cl$spe, verbose = FALSE))[[1]]
  md <- fits(polishSpiDE(cl$fit, cl$spe, psi = "moderated", verbose = FALSE))[[1]]
  expect_equal(unname(md@psi), unname(md@polish$psi_fitnb))
  expect_false(isTRUE(all.equal(unname(pr@psi), unname(pr@polish$psi_fitnb))))
  expect_lt(median(pr@psi), median(md@psi))   # the shared fit's dispersion is inflated here
})

test_that("spiDE(polish = FALSE) is the fit and the test without the stage", {
  spe <- buildNiches(spiDE:::.toySPE(), sigma = 30, verbose = FALSE)
  a <- spiDE(spe, "condition", sigma = 30, polish = FALSE, fdr = 1, verbose = FALSE)
  b <- testSpiDE(fitSpiDE(spe, "condition", sigma = 30, verbose = FALSE), spe = spe, fdr = 1)
  expect_equal(results(a), results(b))
  expect_null(fits(a)[[1]]@polish)
  d <- spiDE(spe, "condition", sigma = 30, fdr = 1, verbose = FALSE)
  e <- testSpiDE(polishSpiDE(fitSpiDE(spe, "condition", sigma = 30, verbose = FALSE), spe, verbose = FALSE), spe = spe, fdr = 1)
  expect_equal(results(d), results(e))
  expect_false(is.null(fits(d)[[1]]@polish))
})
