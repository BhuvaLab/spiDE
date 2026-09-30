# The patient-level plots (0.99.34): each patient's slopes, and what the
# patient intercepts (or slopes) capture.

spe <- buildNiches(.toySPE(n_samples = 16, n_per = 150, seed = 7), sigma = 30)
fit <- fitSpiDE(spe, condition = "condition", sigma = 30, verbose = FALSE)
res <- testSpiDE(fit, procedure = "all")
builds <- function(p) {
  expect_true(inherits(p, "ggplot"))
  expect_no_warning(b <- ggplot2::ggplot_build(p))
  invisible(b)
}

test_that("plotPatientSlopes() draws one column per gene, weights summing to one", {
  p <- plotPatientSlopes(res, gene = c("G1", "G2"), index = "A", niche = "B")
  b <- builds(p)
  expect_true(all(c("gene", "patient", "slope", "se", "weight", "condition", "gene_lab") %in% names(p$data)))
  expect_equal(as.numeric(tapply(p$data$weight, p$data$gene, sum)), c(1, 1), tolerance = 1e-12)
  expect_equal(length(unique(b$layout$layout$gene_lab)), 2L)
  expect_setequal(levels(p$data$condition), c("Non-responder", "Responder"))
})

test_that("plotPatientSlopes() without a condition, and a patient without a slope", {
  f0 <- fitSpiDE(spe, sigma = 30, index = "A", verbose = FALSE)
  f0@index$A$beta["G1", "S3", ] <- NA
  p <- plotPatientSlopes(testSpiDE(f0), gene = "G1", index = "A", niche = "B")
  builds(p)
  expect_false("S3" %in% p$data$patient)
  expect_equal(levels(p$data$condition), "all patients")
})

test_that("plotPatientSlopes() refuses the sandwich engine and unknown triplets", {
  fs <- fitSpiDE(spe, condition = "condition", sigma = 30, index = "A", engine = "sandwich", verbose = FALSE)
  expect_error(plotPatientSlopes(testSpiDE(fs), "G1", "A", "B"), "slopes engine")
  expect_error(plotPatientSlopes(res, "G1", "A", "Z"), "niches tested in A")
  expect_error(plotPatientSlopes(res, "nope", "A", "B"), "not tested in A: nope")
})

test_that("plotPatientEffects() relates intercept PCs to patient and derived covariates", {
  p <- plotPatientEffects(res, index = "A")
  builds(p)
  expect_true(all(c("pc", "covariate", "r2", "group") %in% names(p$data)))
  expect_true(all(is.na(p$data$r2) | (p$data$r2 >= 0 & p$data$r2 <= 1)))
  expect_true(all(c("condition", "cells (log10)", "depth (mean log)", "B density") %in% p$data$covariate))
  expect_false(any(c("patient", "sample_id") %in% p$data$covariate))
})

test_that("plotPatientEffects() on the slopes, as a PCA, and with chosen covariates", {
  builds(plotPatientEffects(res, index = "A", niche = "B"))
  p <- plotPatientEffects(res, index = "A", type = "pca", colour.by = "condition")
  builds(p)
  expect_true(all(c("patient", "PC1", "PC2", "condition") %in% names(p$data)))
  p <- plotPatientEffects(fit, index = "A", covariates = c("condition", "B density"))
  expect_setequal(unique(as.character(p$data$covariate)), c("condition", "B density"))
  expect_error(plotPatientEffects(res, index = "A", covariates = "nope"), "nope")
})

test_that("plotPatientEffects() refuses a fit without intercepts, and slopes of the sandwich engine", {
  f <- fit
  f@index$A$intercept <- NULL
  expect_error(plotPatientEffects(f, index = "A"), "refit")
  fs <- fitSpiDE(spe, condition = "condition", sigma = 30, index = "A", engine = "sandwich", verbose = FALSE)
  builds(plotPatientEffects(fs, index = "A"))
  expect_error(plotPatientEffects(fs, index = "A", niche = "B"), "slopes engine")
})

test_that("plotPatientEffects() keeps its patients when many genes were not fitted", {
  f <- fit
  bad <- f@index$A$genes[1:6]
  f@index$A$intercept[bad, ] <- NA
  p <- plotPatientEffects(f, index = "A", type = "pca")
  expect_equal(nrow(p$data), length(f@index$A$patients))
})
