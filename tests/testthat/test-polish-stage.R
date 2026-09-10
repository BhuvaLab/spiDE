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

test_that("the diagnostics of a re-estimating polish record fitNB's dispersion and the re-polish work", {
  cl <- clustered()
  f <- fits(polishSpiDE(cl$fit, cl$spe, verbose = FALSE))[[1]]
  # psi_fitnb is the shared fit's dispersion, not the previous pass's polished
  # value (the 0.99.19 loop handed each re-polish the last pass's psi as psi0,
  # so the driver's "polished / fitNB" ratio read ~1)
  expect_equal(unname(f@polish$psi_fitnb), unname(fits(cl$fit)[[1]]@psi))
  expect_true("repolish.iterations" %in% colnames(f@polish))
  expect_true(all(f@polish$repolish.iterations >= 0L))
  expect_true(any(f@polish$repolish.iterations > 0L))
})

# The variance-component loop. Schall's update is a fixed-point iteration that
# converges linearly, and on the 0.99.19 cohort runs the nested component was
# still 6-12% above its extrapolated limit when the cap of three stopped it.
# That is not a rounding matter: on the clustered fixture, a tenfold error in
# a near-zero nested component moves individual t-statistics by up to 0.13
# (median 0.012), a 5% error in the per-sample one by up to 0.05.

test_that(".tau2Iterate reaches a geometric fixed point in a few steps and the plain map in many", {
  target <- c(SampleInt = 0.24, SampleCellTypeInt = 0.05)
  contract <- function(tau2) {
    x <- unlist(tau2)
    as.list(target + 0.6 * (x - target))
  }
  acc <- spiDE:::.tau2Iterate(list(SampleInt = 1, SampleCellTypeInt = 1), contract,
                              maxit = 20L, tol = 1e-2, accelerate = TRUE)
  plain <- spiDE:::.tau2Iterate(list(SampleInt = 1, SampleCellTypeInt = 1), contract,
                                maxit = 20L, tol = 1e-2, accelerate = FALSE)
  # an exactly geometric sequence is extrapolated to its limit by Aitken
  expect_lt(max(abs(log(unlist(acc$tau2)) - log(target))), 1e-6)
  expect_lte(acc$iterations, 5L)
  expect_gt(plain$iterations, 8L)
  expect_lt(max(abs(log(unlist(plain$tau2)) - log(target))), 2e-2)
})

test_that(".tau2Iterate leaves an oscillating sequence to the plain map", {
  target <- c(SampleInt = 0.24, SampleCellTypeInt = 0.05)
  flip <- function(tau2) as.list(target - 0.5 * (unlist(tau2) - target))
  acc <- spiDE:::.tau2Iterate(list(SampleInt = 0.4, SampleCellTypeInt = 0.08), flip,
                              maxit = 30L, tol = 1e-3, accelerate = TRUE)
  expect_true(all(is.finite(unlist(acc$tau2))))
  expect_lt(max(abs(log(unlist(acc$tau2)) - log(target))), 2e-3)
})

test_that("the polish converges both variance components, not just the first", {
  # a nested (sample x cell type) component planted at 0.3^2 = 0.09 beside the
  # per-sample 0.49, so both fixed points are interior. The reference is the
  # plain Schall map run to a tight tolerance; the default loop must land
  # within 1e-2 on the log scale of it on BOTH components, where the capped
  # three-step loop stopped short on the nested one.
  spe <- buildNiches(spiDE:::.toyClustered(n_samples = 8, n_per = 40, n_genes = 12,
                                            sd_patient = 0.7, sd_nested = 0.3),
                     sigma = 30, verbose = FALSE)
  fit <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept", verbose = FALSE)
  ref <- fits(polishSpiDE(fit, spe, tau2.maxit = 40L, tau2.tol = 1e-4,
                          tau2.accelerate = FALSE, verbose = FALSE))[[1]]@tau2
  msgs <- capture_messages(def <- polishSpiDE(fit, spe, verbose = TRUE))
  got <- fits(def)[[1]]@tau2
  expect_lt(abs(log(got[["SampleInt"]]) - log(ref[["SampleInt"]])), 1e-2)
  expect_lt(abs(log(got[["SampleCellTypeInt"]]) - log(ref[["SampleCellTypeInt"]])), 1e-2)
  # and it got there in fewer variance-component steps than the reference
  expect_lt(sum(grepl("tau2 from the converged fit", msgs)), 12L)
  # the penalty and the components agree at the end
  f <- fits(def)[[1]]
  for (g in names(got)) {
    expect_equal(unname(unique(f@penalty[f@re_group %in% g])), 1 / got[[g]])
  }
})

# Three points from the numerical-robustness review of the loop (2026-09-10).

test_that(".tau2Iterate stops with a diagnosis on a non-finite component", {
  # tau2 is shared across every gene of the bandwidth, so a bad value must
  # stop with a cause, not surface as "missing value where TRUE/FALSE needed"
  bad <- function(tau2) list(SampleInt = 0.3, SampleCellTypeInt = NaN)
  expect_error(spiDE:::.tau2Iterate(list(SampleInt = 1, SampleCellTypeInt = 1), bad,
                                    maxit = 3L),
               "non-finite")
})

test_that("a singular information at the end of the loop warns that the df was not refreshed", {
  cl <- clustered()
  testthat::local_mocked_bindings(invert_mat = function(...) stop("singular"),
                                  .package = "SpaNorm")
  w <- capture_warnings(pol <- polishSpiDE(cl$fit, cl$spe, verbose = FALSE))
  expect_match(paste(w, collapse = " "), "reference df")
  f <- fits(pol)[[1]]
  expect_true(is.numeric(f@df) && !is.null(names(f@df)) && all(is.finite(f@df)))
})

test_that("the diagnostics flag a re-polish that hit its cap", {
  cl <- clustered()
  f <- fits(polishSpiDE(cl$fit, cl$spe, verbose = FALSE))[[1]]
  expect_true(all(c("repolish.capped", "repolish.singular") %in% colnames(f@polish)))
  expect_false(any(f@polish$repolish.capped))
  expect_false(any(f@polish$repolish.singular))
  # tol = 0 makes every Newton call run to its cap, warm passes included
  g <- fits(polishSpiDE(cl$fit, cl$spe, maxit = 2L, tol = 0, verbose = FALSE))[[1]]
  expect_true(all(g@polish$capped))
  expect_true(all(g@polish$repolish.capped))
})
