# polishSpiDE(): the per-gene convergence stage as a post-hoc adjustment on an
# existing fit. It must be the SAME stage fitSpiDE(converge = TRUE) runs.

spe_p <- buildNiches(.toySPE(), sigma = 20)

test_that("polishSpiDE(fitSpiDE(converge = FALSE)) is fitSpiDE(converge = TRUE)", {
  # fitNB subsamples cells for dispersion without seeding itself, so seed both
  # routes identically and the only difference left is WHERE the polish runs
  set.seed(11)
  a <- fitSpiDE(spe_p, "condition", sigma = 20, random = "intercept",
                converge = TRUE, verbose = FALSE)
  set.seed(11)
  b <- fitSpiDE(spe_p, "condition", sigma = 20, random = "intercept",
                converge = FALSE, verbose = FALSE)
  expect_null(fits(b)[[1]]@polish)
  b <- polishSpiDE(b, spe_p, verbose = FALSE)
  fa <- fits(a)[[1]]; fb <- fits(b)[[1]]
  expect_equal(fb@alpha, fa@alpha, tolerance = 1e-10)
  expect_equal(fb@psi, fa@psi, tolerance = 1e-10)
  expect_equal(fb@polish$iterations, fa@polish$iterations)
  expect_equal(rownames(fb@polish), rownames(fb@alpha))
  expect_equal(fb@loglik, fa@loglik, tolerance = 1e-8)
})

test_that("polishing clears inference, and testSpiDE recomputes it", {
  set.seed(11)
  r <- spiDE(spe_p, "condition", sigma = 20, random = "none", converge = FALSE,
             fdr = 1, verbose = FALSE)
  expect_false(is.null(fits(r)[[1]]@t_stat))
  expect_true(nrow(results(r)) > 0)
  p <- polishSpiDE(r, spe_p, verbose = FALSE)
  expect_null(fits(p)[[1]]@t_stat)
  expect_null(fits(p)[[1]]@se)
  expect_equal(nrow(results(p)), 0L)
  expect_null(p@gene.weights)
  p <- testSpiDE(p, spe = spe_p, fdr = 1)
  expect_true(nrow(results(p)) > 0)
  # the planted effect is at least as sharp after polishing (it was 1.68 ->
  # 10.19 on this fixture in condition mode with covariates; here no covariate)
  ab <- results(r); ab <- ab[ab$gene == "G1" & ab$ct_index == "A" & ab$ct_niche == "B", "t"]
  ab2 <- results(p); ab2 <- ab2[ab2$gene == "G1" & ab2$ct_index == "A" & ab2$ct_niche == "B", "t"]
  expect_gt(abs(ab2), abs(ab) * 0.9)
})

test_that("a fixed-effects fit polishes against the lambda.a it was made with", {
  set.seed(3)
  f0 <- fitSpiDE(spe_p, "condition", sigma = 20, random = "none", lambda.a = 0.5,
                 converge = FALSE, verbose = FALSE)
  set.seed(3)
  f1 <- fitSpiDE(spe_p, "condition", sigma = 20, random = "none", lambda.a = 0.5,
                 converge = TRUE, verbose = FALSE)
  p_right <- polishSpiDE(f0, spe_p, lambda.a = 0.5, verbose = FALSE)
  p_wrong <- polishSpiDE(f0, spe_p, lambda.a = 0, verbose = FALSE)
  expect_equal(fits(p_right)[[1]]@alpha, fits(f1)[[1]]@alpha, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(fits(p_wrong)[[1]]@alpha, fits(f1)[[1]]@alpha,
                                tolerance = 1e-6)))
})

test_that("re-polishing an already converged fit changes nothing material", {
  set.seed(5)
  f1 <- fitSpiDE(spe_p, "condition", sigma = 20, random = "none",
                 converge = TRUE, verbose = FALSE)
  f2 <- polishSpiDE(f1, spe_p, verbose = FALSE)
  # The stopping rule is a RELATIVE log-likelihood gain of 1e-8, which leaves
  # coefficients ~1e-4 off the optimum; a second polish closes that and no
  # more. Unidentified intercepts (a cell type with no counts for the gene sit
  # at log-mu ~ -40 to -60) are compared on the same relative scale.
  expect_equal(fits(f2)[[1]]@alpha, fits(f1)[[1]]@alpha, tolerance = 1e-3)
  expect_equal(fits(f2)[[1]]@psi, fits(f1)[[1]]@psi, tolerance = 1e-3)
  # and the objective did not get worse
  expect_true(all(fits(f2)[[1]]@loglik >= fits(f1)[[1]]@loglik - 1e-6 * abs(fits(f1)[[1]]@loglik)))
})

test_that("an object with no GLM fit is refused clearly", {
  set.seed(5)
  f <- fitSpiDE(spe_p, "condition", sigma = 20, random = "none", converge = FALSE,
                verbose = FALSE)
  f@fits <- list()
  expect_error(polishSpiDE(f, spe_p, verbose = FALSE), "nothing to polish")
})

test_that("a fit serialised before the polish slot existed can be polished", {
  p <- system.file("extdata", "testfits", "fit_toyspe_none.rds", package = "spiDE")
  skip_if(!nzchar(p), "fixture not installed")
  old <- readRDS(p)
  skip_if("polish" %in% names(attributes(old)), "fixture was regenerated")
  # wrap the bare SpiDEFit the way fitSpiDE() returns it
  res <- new("SpiDEResults", fits = list(old), sigma = old@sigma, condition = "condition",
             mode = old@mode, index = character(), niche = character(),
             covariates = character(), coldata = SummarizedExperiment::colData(spe_p),
             gene.weights = NULL, p.cauchy.pos = NULL, p.cauchy.neg = NULL,
             results = data.frame(), fdr = NA_real_, call = quote(x))
  spe_old <- buildNiches(.toySPE(), sigma = old@sigma)
  out <- polishSpiDE(res, spe_old, verbose = FALSE)
  expect_false(is.null(fits(out)[[1]]@polish))
  expect_equal(nrow(fits(out)[[1]]@polish), old@ngenes)
})
