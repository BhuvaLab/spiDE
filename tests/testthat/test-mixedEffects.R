# Tests for the mixed-effects (random-effects via ridge) correction for
# cell-level pseudo-replication. See R/fitSpiDE.R (.fitNBmixed) and R/inference.R.

test_that("random='none' is unchanged: no RE slots, results identical", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  a <- fitSpiDE(spe, "condition", sigma = 20, verbose = FALSE)
  fa <- fits(a)[[1]]
  expect_null(fa@re_group)
  expect_null(fa@tau2)
  expect_null(fa@penalty)
  expect_null(fa@df)
  # a random fit adds the Random covtype and populates the slots
  b <- fitSpiDE(spe, "condition", sigma = 20, random = "slope", verbose = FALSE)
  fb <- fits(b)[[1]]
  expect_true("Random" %in% levels(fb@covtype))
  expect_true(all(c("SampleInt", "SampleSlope") %in% fb@re_group))
  expect_equal(length(fb@penalty), ncol(fb@W))
  expect_true(is.finite(fb@df))
})

test_that("random-effect design targets the response-related effects", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  des <- spiDE:::.buildNicheDesign(spe, "condition", 20, random = "slope")
  reg <- des$re_group
  n_samples <- length(unique(spe$sample_id))
  # one random intercept per sample (counterpart of the response main effect)
  expect_equal(sum(reg == "SampleInt", na.rm = TRUE), n_samples)
  # random slopes = samples x (CellType:niche base terms) i.e. the non-response
  # bases of the ResponseNiche terms
  n_niche_base <- sum(des$covtype == "Niche")
  expect_equal(sum(reg == "SampleSlope", na.rm = TRUE), n_samples * n_niche_base)
})

test_that("checkSample rejects within-sample condition and sample-level covariates", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  # a covariate constant within sample is confounded with the random intercept
  spe$Age <- ave(rnorm(ncol(spe)), spe$sample_id)
  expect_error(
    fitSpiDE(spe, "condition", sigma = 20, covariates = "Age",
             random = "intercept", verbose = FALSE),
    "constant within sample"
  )
  # a condition that varies within a sample is not a patient-level label
  spe$bad <- rep(c("x", "y"), length.out = ncol(spe))
  expect_error(
    fitSpiDE(spe, "bad", sigma = 20, random = "intercept", verbose = FALSE),
    "varies within sample"
  )
})

test_that(".stratifiedCellIdx honours the proportion, floor and cap per stratum", {
  set.seed(1)
  # three samples x two cell types; stratum sizes span the floor and the prop
  cell_type <- rep(c("A", "B"), times = c(1500, 40))
  sample <- rep(c("S1", "S2", "S3"), length.out = length(cell_type))
  cell_type <- sample(cell_type) # shuffle so strata interleave
  idx <- spiDE:::.stratifiedCellIdx(cell_type, sample, prop = 0.1,
                                    min.cells = 100L)
  expect_type(idx, "logical")
  expect_length(idx, length(cell_type))

  strata <- interaction(sample, cell_type, drop = TRUE)
  for (s in levels(strata)) {
    n <- sum(strata == s)
    kept <- sum(idx[strata == s])
    expect_equal(kept, min(n, max(ceiling(0.1 * n), 100L)))
  }
  # every stratum keeps at least one cell (samples stay represented)
  expect_true(all(tapply(idx, strata, any)))

  # prop >= 1 keeps everything
  expect_true(all(spiDE:::.stratifiedCellIdx(cell_type, sample, prop = 1)))
})

test_that("re.prop is validated and small strata are taken whole", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  expect_error(
    fitSpiDE(spe, "condition", sigma = 20, random = "intercept",
             re.prop = 1.5, verbose = FALSE),
    "re.prop"
  )
  # the toy strata (~27 cells) are below the 100-cell floor, so re.prop must not
  # change the fit relative to using all cells
  fa <- fitSpiDE(spe, "condition", sigma = 20, random = "intercept",
                 re.prop = 0.1, verbose = FALSE)
  fb <- fitSpiDE(spe, "condition", sigma = 20, random = "intercept",
                 re.prop = 1, verbose = FALSE)
  expect_equal(fits(fa)[[1]]@tau2, fits(fb)[[1]]@tau2)
  expect_equal(fits(fa)[[1]]@alpha, fits(fb)[[1]]@alpha)
})

test_that("the mixed fit recovers the planted between-sample variance", {
  spe <- buildNiches(.toyClustered(sd_patient = 0.7), sigma = 30)
  fit <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                  verbose = FALSE)
  tau2 <- fits(fit)[[1]]@tau2[["SampleInt"]]
  # planted variance is 0.7^2 = 0.49; Schall/PQL should land close
  expect_gt(tau2, 0.3)
  expect_lt(tau2, 0.75)
})

test_that("mixed effects deflate the anti-conservative response inference", {
  spe <- buildNiches(.toyClustered(sd_patient = 0.7), sigma = 30)
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))

  resp_p <- function(random) {
    f <- fitSpiDE(spe, "condition", sigma = 30, random = random, verbose = FALSE)
    ff <- spiDE:::.blockedInference(fits(f)[[1]], Y)
    rc <- ff@coefmap$covariate[as.character(ff@covtype) == "Response"]
    t <- ff@t_stat[, rc]
    df <- if (is.null(ff@df)) Inf else ff@df
    2 * stats::pt(-abs(t), df)
  }
  p_fixed <- resp_p("none")
  p_mixed <- resp_p("intercept")

  # all genes are null for the response effect: the fixed fit is badly
  # anti-conservative, the mixed fit is far closer to nominal
  expect_gt(mean(p_fixed < 0.05), 0.30)
  expect_lt(mean(p_mixed < 0.05), 0.20)
  expect_lt(mean(p_mixed < 0.05), mean(p_fixed < 0.05))
})
