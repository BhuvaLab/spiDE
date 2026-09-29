# The toy fixture (R/toydata.R): the composition confound the between-patient
# niche-abundance test is checked on.

test_that(".toySPE(composition = 0) is unchanged and composition plants a between-sample effect", {
  a <- spiDE:::.toySPE()
  b <- spiDE:::.toySPE(composition = 0)
  expect_identical(SummarizedExperiment::assay(a, "counts"),
                   SummarizedExperiment::assay(b, "counts"))

  cs <- spiDE:::.toySPE(n_samples = 12, composition = 3)
  cd <- SummarizedExperiment::colData(cs)
  y <- SummarizedExperiment::assay(cs, "counts")["G2", ]
  isA <- cd$cell_type == "A"
  resp <- cd$condition == "Responder"

  # G2 in Responders' A cells differs BETWEEN samples ...
  m <- tapply(y[isA & resp], droplevels(factor(cd$sample_id[isA & resp])), mean)
  expect_gt(max(m) / min(m), 1.5)
  # ... in the same order as how close the sample's A cells sit to the B-rich
  # region (the per-sample shift the confound is planted on), so the sample's
  # MEAN B-niche density around its A cells is what G2 tracks
  mx <- tapply(cd$x[isA & resp], droplevels(factor(cd$sample_id[isA & resp])), mean)
  expect_gt(cor(as.numeric(m), as.numeric(mx[names(m)]), method = "spearman"), 0.5)
  # ... and the DEFAULT fixture draws none of this: no sample's A cells are
  # shifted away from the left edge of the field
  cd0 <- SummarizedExperiment::colData(a)
  isA0 <- cd0$cell_type == "A"
  expect_true(all(tapply(cd0$x[isA0], cd0$sample_id[isA0], min) < 0.3 * 500))
})

# BiocCheck flags the set.seed() inside .localSeed(); it must stay, because
# the helper restores the caller's RNG stream. withr::local_seed() (3.0.2) does
# not restore it in this use (measured 2026-09-30), so it is no substitute.
test_that("a seeded toy generator leaves the user's random stream untouched", {
  set.seed(3); x1 <- runif(1)
  set.seed(3); invisible(spiDE:::.toySPE(n_samples = 2, n_per = 20, n_genes = 2, seed = 1)); x2 <- runif(1)
  expect_identical(x1, x2)
})
