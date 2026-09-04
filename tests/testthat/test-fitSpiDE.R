test_that("fitSpiDE returns a SpiDEResults with one fit per bandwidth", {
  spe <- buildNiches(.toySPE(), sigma = c(10, 20))
  res <- fitSpiDE(spe, condition = "condition", random = "none", verbose = FALSE)

  expect_s4_class(res, "SpiDEResults")
  expect_equal(bandwidths(res), c(10, 20))
  expect_length(fits(res), 2)
  expect_true(all(vapply(fits(res), validObject, logical(1))))
})

test_that("the NB fit produces valid dispersions and log-likelihoods", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  res <- fitSpiDE(spe, condition = "condition", sigma = 20, random = "none", verbose = FALSE)
  f <- fits(res)[[1]]

  expect_equal(dim(f@alpha), c(f@ngenes, ncol(f@W)))
  expect_false(anyNA(f@psi))
  expect_true(all(f@psi > 0))
  expect_true(all(is.finite(f@loglik)))
  expect_true(all(f@loglik <= 0))
})

test_that("fitSpiDE recovers the planted B-niche effect on G1 in A cells", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  spe <- computeSizeFactors(spe, count = "nCount", area = "Area")
  res <- fitSpiDE(spe,
    condition = "condition", sigma = 20, random = "none",
    covariates = c("Age", "LS"), verbose = FALSE
  )
  f <- fits(res)[[1]]
  cm <- f@coefmap
  col <- cm$covariate[cm$type == "ResponseNiche" &
    cm$index == "A" & cm$niche == "B"]

  # The planted effect is recovered as a STATISTIC, not as a raw coefficient.
  # This used to assert which.max(alpha[, col]) == "G1", which is the wrong
  # quantity: a near-empty gene (G10 here, mean count 0.26) can carry a larger
  # point estimate than the planted one on an estimate its own standard error
  # swamps. Converging each gene made that explicit -- G1's dispersion falls
  # from 3.09 to 0.28 once its dynamic range is actually fitted, which is the
  # documented .toySPE() pathology -- so its coefficient shrinks while its t
  # statistic goes from 1.68 (fdr .09, not significant) to 10.19 (fdr 2e-15).
  # Assert the t.
  tab <- results(testSpiDE(res, spe = spe, fdr = 1))
  ab <- tab[tab$ct_index == "A" & tab$ct_niche == "B", ]
  expect_equal(ab$gene[which.max(abs(ab$t))], "G1")
  expect_gt(abs(ab$t[ab$gene == "G1"]), 5)
  expect_gt(f@alpha["G1", col], 0)
})

test_that("fitSpiDE errors when niches are missing", {
  spe <- .toySPE()
  expect_error(fitSpiDE(spe, condition = "condition"), "buildNiches|niche")
})

test_that(".toySPE(composition = 0) is unchanged and composition plants a between-sample effect", {
  a <- spiDE:::.toySPE()
  b <- spiDE:::.toySPE(composition = 0)
  expect_identical(SummarizedExperiment::assay(a, "counts"),
                   SummarizedExperiment::assay(b, "counts"))

  cs <- spiDE:::.toySPE(composition = 2.5)
  cd <- SummarizedExperiment::colData(cs)
  y <- SummarizedExperiment::assay(cs, "counts")["G2", ]
  isA <- cd$cell_type == "A"

  # G2 in A cells differs BETWEEN samples ...
  m <- tapply(y[isA], droplevels(factor(cd$sample_id[isA])), mean)
  expect_gt(max(m) / min(m), 1.5)
  # ... and B-cell prevalence differs between samples too, in the same order
  pb <- tapply(cd$cell_type == "B", factor(cd$sample_id), mean)
  expect_gt(abs(cor(as.numeric(m), as.numeric(pb[names(m)]),
                    method = "spearman")), 0.6)
})

test_that("a covariate with non-finite values is refused by name", {
  # model.matrix() drops those rows, so the design no longer matches the
  # random-effect block and the run died with "number of rows of matrices must
  # match" -- an error naming neither the covariate nor the cause.
  spe <- buildNiches(.toySPE(), sigma = 20)
  SummarizedExperiment::colData(spe)$bad <- log(c(0, runif(ncol(spe) - 1)))
  expect_error(
    fitSpiDE(spe, "condition", sigma = 20, covariates = "bad",
             random = "intercept", verbose = FALSE),
    "non-finite")
  expect_error(
    fitSpiDE(spe, "condition", sigma = 20, covariates = "bad",
             random = "intercept", verbose = FALSE),
    "bad")
})
