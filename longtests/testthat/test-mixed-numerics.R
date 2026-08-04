# Numerical behaviour of the mixed-effects fit, with LIVE fits.
#
# These were in tests/ and cost ~2,200 s between them -- each alone exceeded
# Bioconductor's 10-minute budget for `R CMD check --no-build-vignettes`. They
# are here rather than precomputed because every one of them asserts something
# about what fitSpiDE() *computes*, not about the shape of what it returns.
# Freezing these into fixtures would let a regression in the fitting code pass
# unnoticed, which is the one failure they exist to catch. The shape assertions
# that CAN be frozen stayed behind in tests/, against inst/extdata/testfits.
#
# The Long Tests builder runs weekly with a 6-hour budget.

test_that("df.method='satterthwaite' anchors the between-sample contrast at S-2", {
  spe <- buildNiches(
    spiDE:::.toyClustered(n_samples = 16, n_per = 30, n_genes = 20,
                          sd_patient = 0.7), sigma = 30)
  fs <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                 df.method = "satterthwaite", verbose = FALSE)
  ff <- fits(fs)[[1]]
  ds <- ff@df
  ct <- as.character(ff@covtype)
  cm <- ff@coefmap

  # Under cell-means coding the between-sample condition contrast is carried by
  # the CellType:condition columns, one per cell type, not by a single
  # "Response" main effect. The claim is unchanged: a between-sample contrast
  # gets df ~ S - 2, a within-sample one gets more.
  resp_name <- cm$covariate[ct == "ResponseCellType"]
  expect_gt(length(resp_name), 0L)
  expect_equal(unname(stats::median(ds[resp_name])), 14, tolerance = 0.20)

  rn_names <- cm$covariate[ct == "ResponseNiche"]
  expect_gt(stats::median(ds[rn_names]), stats::median(ds[resp_name]))
})

test_that("df.method='between' reproduces the scalar S-2 reference exactly", {
  spe <- buildNiches(
    spiDE:::.toyClustered(n_samples = 12, n_per = 30, n_genes = 20,
                          sd_patient = 0.7), sigma = 30)
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  fb <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                 df.method = "between", verbose = FALSE)
  expect_equal(fits(fb)[[1]]@df, 10)                    # 12 samples -> S - 2
  tb <- spiDE:::.blockedInference(fits(fb)[[1]], Y)
  expect_false(any(is.na(tb@p.combined.pos)))
})

test_that("satterthwaite is the df.method default on an unqualified mixed fit", {
  # guards the 0.99.7 behaviour change. Must be a live fit: the point is what
  # fitSpiDE() does when df.method is NOT supplied, which a stored object
  # cannot demonstrate.
  spe <- buildNiches(
    spiDE:::.toyClustered(n_samples = 12, n_per = 30, n_genes = 20,
                          sd_patient = 0.7), sigma = 30)
  fd <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                 verbose = FALSE)
  ff <- fits(fd)[[1]]
  expect_gt(length(ff@df), 1L)
  expect_equal(length(ff@df), sum(grepl("Response", as.character(ff@covtype))))
  expect_true(all(is.finite(ff@df)))
})

test_that("the mixed fit recovers the planted between-sample variance", {
  spe <- buildNiches(spiDE:::.toyClustered(n_genes = 15, sd_patient = 0.7),
                     sigma = 30)
  fit <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                  df.method = "between", verbose = FALSE)
  tau2 <- fits(fit)[[1]]@tau2[["SampleInt"]]
  # planted variance is 0.7^2 = 0.49; Schall/PQL should land close
  expect_gt(tau2, 0.3)
  expect_lt(tau2, 0.75)
})

test_that("mixed effects deflate the anti-conservative response inference", {
  # n_per is deliberately left at the default: cells-per-sample IS the
  # pseudo-replication this test measures, so shrinking it would weaken the
  # phenomenon rather than just the runtime.
  spe <- buildNiches(spiDE:::.toyClustered(n_genes = 15, sd_patient = 0.7),
                     sigma = 30)
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))

  resp_p <- function(random) {
    f <- fitSpiDE(spe, "condition", sigma = 30, random = random,
                  df.method = "between", verbose = FALSE)
    ff <- spiDE:::.blockedInference(fits(f)[[1]], Y)
    rc <- ff@coefmap$covariate[as.character(ff@covtype) == "ResponseCellType"]
    tt <- ff@t_stat[, rc]
    df <- if (is.null(ff@df)) Inf else ff@df
    2 * stats::pt(-abs(tt), df)
  }
  p_fixed <- resp_p("none")
  p_mixed <- resp_p("intercept")

  # all genes are null for the response effect: the fixed fit is badly
  # anti-conservative, the mixed fit is far closer to nominal
  expect_gt(mean(p_fixed < 0.05), 0.30)
  expect_lt(mean(p_mixed < 0.05), 0.20)
  expect_lt(mean(p_mixed < 0.05), mean(p_fixed < 0.05))
})

test_that("df.method='satterthwaite' works for random='slope'", {
  spe <- buildNiches(
    spiDE:::.toyClustered(n_samples = 16, n_per = 30, n_genes = 20,
                          sd_patient = 0.7), sigma = 30)
  fs <- fitSpiDE(spe, "condition", sigma = 30, random = "slope",
                 df.method = "satterthwaite", verbose = FALSE)
  ds <- fits(fs)[[1]]@df
  ct <- as.character(fits(fs)[[1]]@covtype)
  expect_length(ds, sum(grepl("Response", ct)))
  expect_true(all(is.finite(ds)) && all(ds >= 1))
})

test_that("re.prop below the per-stratum floor leaves the fit unchanged", {
  # the toy strata (~27 cells) are below the 100-cell floor, so subsampling must
  # not change tau2 or the coefficients. Live fits: freezing them would let a
  # regression in the subsampling path pass unnoticed.
  spe <- buildNiches(spiDE:::.toySPE(), sigma = 20)
  fa <- fitSpiDE(spe, "condition", sigma = 20, random = "intercept",
                 re.prop = 0.1, df.method = "between", verbose = FALSE)
  fb <- fitSpiDE(spe, "condition", sigma = 20, random = "intercept",
                 re.prop = 1, df.method = "between", verbose = FALSE)
  expect_equal(fits(fa)[[1]]@tau2, fits(fb)[[1]]@tau2)
  expect_equal(fits(fa)[[1]]@alpha, fits(fb)[[1]]@alpha)
})
