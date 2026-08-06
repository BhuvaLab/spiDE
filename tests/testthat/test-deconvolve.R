# deconvolveSpillover(): does inverting the mixing operator actually recover
# the pre-contamination counts?
#
# The round-trip is the real test. .toySpill() APPLIES exactly the operator
# deconvolveSpillover() inverts, so applying one then the other at the same
# kappa and radius must return (close to) the original -- and if the two ever
# drift apart, this fails loudly.

# .toySPE()'s default density is unrepresentative for this operator: 80 cells
# over a 500 um field gives a median nearest-neighbour distance of 26.6 um, so
# at a one-cell-diameter radius only ~34% of cells have ANY neighbour and there
# is almost nothing to mix. The real cohort sits at ~9.6 um. These tests
# therefore use a denser fixture. Nothing here fits a model, so the bandwidth
# tuning CLAUDE.md warns about does not apply.
.denseSPE <- function(...) .toySPE(n_per = 400, ...)

test_that("deconvolution inverts .toySpill on the same kappa and radius", {
  spe <- .denseSPE()
  orig <- as.matrix(SummarizedExperiment::assay(spe, "counts"))

  dirty <- spiDE:::.toySpill(spe, kappa = 0.3, radius = 25)
  back <- deconvolveSpillover(dirty, kappa = 0.3, radius = 25,
                              verbose = FALSE)
  rec <- as.matrix(SummarizedExperiment::assay(back, "counts_deconv"))

  # contamination genuinely moved the data
  dirty_m <- as.matrix(SummarizedExperiment::assay(dirty, "counts"))
  expect_gt(mean(dirty_m != orig), 0.15)

  # The inverse recovers most of it. Entry-level error only roughly halves
  # rather than vanishing, because rounding to integers costs up to 0.5 per
  # entry and this fixture's counts are 0/1/2/3 -- quantisation dominates at
  # that depth. The real cohort has a median of 1412 counts per cell, where
  # rounding is negligible. Correlation with the truth is the sharper check.
  err_dirty <- mean(abs(dirty_m - orig))
  err_rec <- mean(abs(rec - orig))
  expect_lt(err_rec, 0.7 * err_dirty)
  expect_gt(cor(as.vector(rec), as.vector(orig)), 0.99)
})

test_that("the deconvolved assay is non-negative integers, written separately", {
  spe <- spiDE:::.toySpill(.denseSPE(), kappa = 0.3, radius = 25)
  out <- deconvolveSpillover(spe, kappa = 0.3, radius = 25, verbose = FALSE)

  # the original assay must survive untouched -- the operation is lossy
  expect_true(all(c("counts", "counts_deconv") %in%
                    SummarizedExperiment::assayNames(out)))
  expect_identical(SummarizedExperiment::assay(out, "counts"),
                   SummarizedExperiment::assay(spe, "counts"))

  x <- as.matrix(SummarizedExperiment::assay(out, "counts_deconv"))
  expect_true(all(x >= 0))
  expect_true(all(x == round(x)))
})

test_that("the QC record reports kappa, radius and the floored fraction", {
  spe <- spiDE:::.toySpill(.denseSPE(), kappa = 0.3, radius = 25)
  out <- deconvolveSpillover(spe, kappa = 0.3, radius = 25, verbose = FALSE)
  qc <- S4Vectors::metadata(out)$spiDE_deconvolution

  expect_equal(qc$kappa, 0.3)
  expect_equal(qc$radius, 25)
  # `floored` counts ENTRIES and tracks sparsity -- on this fixture ~35% floor,
  # of which 98% were observed zeros, where clipping is the right answer. The
  # metric that bounds actual damage is the count MASS destroyed.
  expect_gt(qc$floored, 0)
  expect_lt(qc$floored_mass, 0.10)
  expect_lt(qc$floored_mass, qc$floored)   # mass damage far below entry count
})

test_that("kappa = 0 is a no-op that still writes the assay", {
  spe <- .toySPE()
  out <- deconvolveSpillover(spe, kappa = 0, radius = 20, verbose = FALSE)
  expect_identical(SummarizedExperiment::assay(out, "counts_deconv"),
                   SummarizedExperiment::assay(spe, "counts"))
  expect_equal(S4Vectors::metadata(out)$spiDE_deconvolution$floored, 0)
})

test_that("kappa outside [0, 0.5) is refused", {
  spe <- .toySPE()
  # at kappa >= 0.5 the mixing operator stops being diagonally dominant
  expect_error(deconvolveSpillover(spe, kappa = 0.5, radius = 20),
               "diagonally dominant")
  expect_error(deconvolveSpillover(spe, kappa = -0.1, radius = 20),
               "\\[0, 0.5\\)")
  expect_error(deconvolveSpillover(spe, kappa = 0.2, radius = 0),
               "positive")
})

test_that("cells are never mixed across samples", {
  # Exact isolation: deconvolving the whole object must give a sample bit-for-bit
  # the same answer as deconvolving that sample alone. Any cross-sample edge in
  # the neighbour operator -- samples share a coordinate frame, so cells in
  # different samples sit on top of each other -- would break this immediately.
  spe <- .denseSPE()
  full <- deconvolveSpillover(spe, kappa = 0.3, radius = 25, verbose = FALSE)

  s1 <- which(spe$sample_id == "S1")
  alone <- deconvolveSpillover(spe[, s1], kappa = 0.3, radius = 25,
                               verbose = FALSE)

  expect_equal(
    as.matrix(SummarizedExperiment::assay(full, "counts_deconv"))[, s1],
    as.matrix(SummarizedExperiment::assay(alone, "counts_deconv")),
    ignore_attr = TRUE)
})

test_that("isolated cells pass through untouched", {
  # A cell with no neighbour inside the radius donates nothing and receives
  # nothing, so a tiny radius must be an identity.
  spe <- spiDE:::.toySpill(.denseSPE(), kappa = 0.3, radius = 25)
  out <- deconvolveSpillover(spe, kappa = 0.3, radius = 1e-6, verbose = FALSE)
  expect_identical(
    as.matrix(SummarizedExperiment::assay(out, "counts_deconv")),
    as.matrix(SummarizedExperiment::assay(spe, "counts")))
})

test_that("heavy-tailed genes are flagged by their detection loss", {
  # The operator is a sharpener, so a gene with technical outliers gets its
  # tail pushed further out while its middle is FLOORED -- which shows up as
  # lost detection, not as an inflated maximum. On the cohort this ranked the
  # known-bad genes 1st, 2nd and 7th of 10,422; an earlier guard based on the
  # maximum missed all of them, because deconvolution raises every gene's
  # maximum by roughly the same factor.
  spe <- .denseSPE()
  cnt <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  # G2: low everywhere, enormous in a scattered handful -- the shape that
  # bimodalises. Its neighbours' counts then floor to zero.
  cnt["G2", ] <- 1L
  spikes <- seq(1, ncol(cnt), length.out = 40)
  cnt["G2", round(spikes)] <- 5000L
  SummarizedExperiment::assay(spe, "counts") <- cnt

  # This fixture is deliberately extreme, so the count-mass warning fires as
  # well; collect all warnings and assert the detection one is among them.
  ws <- capture_warnings(
    deconvolveSpillover(spe, kappa = 0.3, radius = 25, verbose = TRUE))
  expect_true(any(grepl("detection", ws)))

  out <- suppressWarnings(
    deconvolveSpillover(spe, kappa = 0.3, radius = 25, verbose = FALSE))
  qc <- S4Vectors::metadata(out)$spiDE_deconvolution
  expect_true("G2" %in% qc$amplified)
  # a well-behaved gene must not be flagged
  expect_false("G5" %in% qc$amplified)
  expect_equal(length(qc$det_drop), nrow(spe))
  expect_gt(qc$det_drop[["G2"]], 0.02)
})
