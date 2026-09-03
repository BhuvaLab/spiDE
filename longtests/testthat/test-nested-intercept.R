# The nested (sample x cell type) intercept and the per-gene convergence stage,
# checked NUMERICALLY on a fixture with a planted between-sample composition
# confound and zero within-sample niche slope. Live mixed fits cost ~70 s each,
# which is why this lives in longtests/ (see data-raw/make_test_fixtures.R for
# the contract/numerics split this project uses).
#
#   Rscript -e 'testthat::test_file("longtests/testthat/test-nested-intercept.R")'

test_that("the nested intercept removes a between-sample composition confound", {
  spe <- buildNiches(spiDE:::.toySPE(composition = 2.5, n_genes = 12),
                     sigma = 30)

  no_nest <- spiDE(spe, condition = "condition", sigma = 30,
                   random = "intercept", re.celltype = FALSE,
                   re.maxit = 2L, fdr = 1, verbose = FALSE)
  nested <- spiDE(spe, condition = "condition", sigma = 30,
                  random = "intercept", re.celltype = TRUE,
                  re.maxit = 2L, fdr = 1, verbose = FALSE)

  pick <- function(res, gene, idx, nch) {
    tab <- results(res)
    r <- tab[tab$gene == gene & tab$ct_index == idx & tab$ct_niche == nch, ]
    expect_equal(nrow(r), 1L)
    r
  }

  # the confound: called without the nested block, attenuated with it
  a <- pick(no_nest, "G2", "A", "B")
  b <- pick(nested, "G2", "A", "B")
  expect_gt(abs(a$t), 3)
  expect_lt(abs(b$t), abs(a$t) / 2)

  # the genuine within-sample effect survives both
  ga <- pick(no_nest, "G1", "A", "B")
  gb <- pick(nested, "G1", "A", "B")
  expect_gt(abs(ga$t), 3)
  expect_gt(abs(gb$t), 3)
  expect_lt(abs(abs(gb$t) - abs(ga$t)) / abs(ga$t), 0.5)
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
    mu <- as.numeric(exp(fit@W %*% fit@alpha[g, ]))
    spiDE:::.nbPenLoglik(Y[g, ], mu, fit@psi[g], fit@alpha[g, ], fit@penalty)
  }
  base <- vapply(seq_len(a1@ngenes), function(g) ll(a0, g), numeric(1))
  gains <- vapply(seq_len(a1@ngenes), function(g) ll(a1, g), numeric(1)) - base
  expect_true(all(gains > -1e-6 * abs(base)))
  expect_gt(median(gains), 0)
  # the dispersion does not inflate; on the real cohort it falls to ~0.63x
  expect_lt(median(a1@psi / a0@psi), 1.05)
})
