# The nested (sample x cell type) intercept and the per-gene convergence stage,
# checked NUMERICALLY on a fixture with a planted between-sample composition
# confound and zero within-sample niche slope. Live mixed fits at 2,400 cells
# cost ~2 min each, which is why this lives in longtests/ (see
# data-raw/make_test_fixtures.R for the contract/numerics split this project
# uses).
#
# Fixture parameters are not arbitrary. At the default 80 cells per sample the
# niche covariate's per-sample MEAN is placement noise (between-sample sd 0.07
# against a within-sample 0.25) and no estimator can see a between-sample
# effect the covariate does not carry -- three redesigns established that. At
# 200 cells and bandwidth 50 the ratio is 0.7; at composition = 6 the old
# design calls the confound at t ~ 3.4 and the nested design reports ~1.3.
#
#   Rscript -e 'testthat::test_file("longtests/testthat/test-nested-intercept.R")'

test_that("the nested intercept removes a between-sample composition confound", {
  spe <- buildNiches(spiDE:::.toySPE(n_samples = 12, n_per = 200, n_genes = 12,
                                     composition = 6, seed = 3),
                     sigma = 50, verbose = FALSE)

  no_nest <- spiDE(spe, condition = "condition", sigma = 50,
                   random = "intercept", re.celltype = FALSE,
                   re.maxit = 2L, fdr = 1, verbose = FALSE)
  nested <- spiDE(spe, condition = "condition", sigma = 50,
                  random = "intercept", re.celltype = TRUE,
                  re.maxit = 2L, fdr = 1, verbose = FALSE)

  pick <- function(res, gene, idx, nch) {
    tab <- results(res)
    r <- tab[tab$gene == gene & tab$ct_index == idx & tab$ct_niche == nch, ]
    expect_equal(nrow(r), 1L)
    r
  }

  # the confound: called without the nested block, not with it
  a <- pick(no_nest, "G2", "A", "B")
  b <- pick(nested, "G2", "A", "B")
  expect_gt(abs(a$t), 2.5)
  expect_lt(a$fdr.niche, 0.05)
  expect_lt(abs(b$t), abs(a$t) / 2)
  expect_gt(b$fdr.niche, 0.05)

  # the genuine within-sample effect survives both
  ga <- pick(no_nest, "G1", "A", "B")
  gb <- pick(nested, "G1", "A", "B")
  expect_gt(abs(ga$t), 3)
  expect_gt(abs(gb$t), 3)
  expect_lt(abs(abs(gb$t) - abs(ga$t)) / abs(ga$t), 0.5)

  # and neither design inflates the null triplets
  for (res in list(no_nest, nested)) {
    tab <- results(res)
    nul <- tab[!(tab$gene %in% c("G1", "G2") & tab$ct_index == "A" &
                   tab$ct_niche == "B"), ]
    expect_lt(sd(nul$t), 1.25)
  }
})

test_that("converging each gene raises every gene's penalised log-likelihood", {
  spe <- buildNiches(spiDE:::.toySPE(n_genes = 15), sigma = 30)
  f0 <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                 re.maxit = 2L, converge = FALSE, verbose = FALSE)
  f1 <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                 re.maxit = 2L, converge = TRUE, verbose = FALSE)
  a0 <- fits(f0)[[1]]
  a1 <- fits(f1)[[1]]
  Y <- SummarizedExperiment::assay(spe, "counts")

  # both fits carry the same penalty vector (same design, same tau2 path), so
  # the comparison is of the same objective at two points
  ll <- function(fit, g) {
    mu <- pmax(as.numeric(exp(fit@W %*% fit@alpha[g, ])), spiDE:::.MU_FLOOR)
    spiDE:::.nbPenLoglik(Y[g, ], mu, fit@psi[g], fit@alpha[g, ], fit@penalty)
  }
  base <- vapply(seq_len(a1@ngenes), function(g) ll(a0, g), numeric(1))
  gains <- vapply(seq_len(a1@ngenes), function(g) ll(a1, g), numeric(1)) - base
  expect_true(all(gains > -1e-6 * abs(base)))
  expect_gt(median(gains), 0)
  # the dispersion does not inflate; on the real cohort it falls to ~0.63x
  expect_lt(median(a1@psi / a0@psi), 1.05)
})

test_that("polishSpiDE() on a converge = FALSE fit is the converge = TRUE fit", {
  spe <- buildNiches(spiDE:::.toySPE(n_genes = 15), sigma = 30)
  set.seed(7)
  f1 <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                 re.maxit = 2L, converge = TRUE, verbose = FALSE)
  set.seed(7)
  f0 <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                 re.maxit = 2L, converge = FALSE, verbose = FALSE)
  f0 <- polishSpiDE(f0, spe, verbose = FALSE)
  expect_equal(fits(f0)[[1]]@alpha, fits(f1)[[1]]@alpha, tolerance = 1e-8)
  expect_equal(fits(f0)[[1]]@psi, fits(f1)[[1]]@psi, tolerance = 1e-8)
})
