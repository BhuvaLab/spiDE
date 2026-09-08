# polishSpiDE(): the per-gene convergence stage as a post-hoc adjustment on an
# existing fit: the stage between the fit and the test.

spe_p <- buildNiches(.toySPE(), sigma = 20)

test_that("polishing clears inference, and testSpiDE recomputes it", {
  set.seed(11)
  r <- spiDE(spe_p, "condition", sigma = 20, random = "none", polish = FALSE,
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
                 verbose = FALSE)
  f1 <- polishSpiDE(f0, spe_p, lambda.a = 0.5, verbose = FALSE)
  p_right <- polishSpiDE(f0, spe_p, lambda.a = 0.5, verbose = FALSE)
  p_wrong <- polishSpiDE(f0, spe_p, lambda.a = 0, verbose = FALSE)
  expect_equal(fits(p_right)[[1]]@alpha, fits(f1)[[1]]@alpha, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(fits(p_wrong)[[1]]@alpha, fits(f1)[[1]]@alpha,
                                tolerance = 1e-6)))
})

test_that("re-polishing an already converged fit changes nothing material", {
  set.seed(5)
  f1 <- polishSpiDE(fitSpiDE(spe_p, "condition", sigma = 20, random = "none", verbose = FALSE),
                    spe_p, verbose = FALSE)
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
  f <- fitSpiDE(spe_p, "condition", sigma = 20, random = "none", verbose = FALSE)
  methods::slot(f, "fits", check = FALSE) <- list()
  # an empty fits list fails validity now (one fit per bandwidth), and the
  # stage's own guard sits behind it; either message is the clear refusal
  expect_error(polishSpiDE(f, spe_p, verbose = FALSE), "nothing to polish|does not match")
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
