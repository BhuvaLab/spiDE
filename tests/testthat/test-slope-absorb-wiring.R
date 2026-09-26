# Wiring the per-sample absorption grouping from the design through to the
# polish and the inference.
#
# SpaNorm::nbNewtonSolver() can absorb a random block that is block-diagonal by sample,
# but nothing reaches that path while the call sites select the absorbed columns
# with `re_group == "SampleCellTypeInt"`. That literal is also exactly what the
# architecture notes warn against: it silently skips a block that a new mode
# adds. The design knows which sample each random column belongs to, so it
# should say so, and the call sites should ask.

# one small slope fit, shared by every test in this file that needs a fit
# (the fit is ~1.5 min of the file's run time; the polish and inference ~1 s)
.slopeWiringFit <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(21)
      spe <- buildNiches(spiDE:::.toySPE(n_genes = 5, n_per = 40), sigma = 30)
      f <- fitSpiDE(spe, condition = "condition", sigma = 30, random = "slope",
                    re.maxit = 1L, verbose = FALSE)@fits[[1]]
      Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))[rownames(f@alpha), ]
      pen <- f@penalty
      if (length(pen) == 1L) pen <- rep(pen, ncol(f@W))
      cache <<- list(f = f, Y = Y, pen = pen)
    }
    cache
  }
})

test_that(".buildRandomEffects reports the sample of every random column", {
  set.seed(4)
  n <- 40
  smp <- factor(rep(c("S1", "S2", "S3", "S4"), each = 10))
  ct <- factor(rep(c("A", "B"), times = 20))
  base <- matrix(stats::rnorm(n * 2), n, 2,
                 dimnames = list(NULL, c("CellTypeA:nicheB", "CellTypeB:nicheA")))

  re <- spiDE:::.buildRandomEffects(smp, base, "slope", cell_type_vec = ct)

  expect_length(re$re_sample, ncol(re$Z))
  expect_false(anyNA(re$re_sample))
  # every random column belongs to exactly one sample, and the sample it names
  # is the one whose cells it is non-zero on
  for (j in seq_len(ncol(re$Z))) {
    hit <- unique(as.character(smp[re$Z[, j] != 0]))
    expect_length(hit, 1L)
    expect_identical(hit, as.character(re$re_sample[j]))
  }
})

test_that("the design's random block is block-orthogonal under its own grouping", {
  # the property SpaNorm::nbNewtonSolver() relies on: no cell loads on two blocks
  spe <- buildNiches(spiDE:::.toySPE(n_genes = 6, n_per = 40), sigma = 30)
  des <- spiDE:::.buildNicheDesign(spe, "condition", 30, random = "slope")

  expect_length(des$re_sample, ncol(des$W))
  zi <- which(!is.na(des$re_sample))
  expect_gt(length(zi), 0L)
  # fixed columns carry no sample
  expect_true(all(is.na(des$re_sample[is.na(des$re_group)])))

  Z <- des$W[, zi, drop = FALSE]
  blk <- as.integer(factor(des$re_sample[zi]))
  hits <- integer(nrow(Z))
  for (b in unique(blk)) {
    hits <- hits + (rowSums(abs(Z[, blk == b, drop = FALSE])) > 0)
  }
  expect_equal(max(hits), 1L)
})

test_that("a slope fit carries the grouping so the polish can absorb it", {
  f1 <- .slopeWiringFit()$f
  expect_length(f1@re_sample, ncol(f1@W))
  expect_true(any(!is.na(f1@re_sample)))
  # and the columns it groups are exactly the random ones
  expect_identical(is.na(f1@re_sample), is.na(f1@re_group))
})

test_that(".absorbSpec picks the whole random block for a slope fit and the nested block alone otherwise", {
  spe <- buildNiches(spiDE:::.toySPE(n_genes = 6, n_per = 40), sigma = 30)

  fi <- fitSpiDE(spe, condition = "condition", sigma = 30, random = "intercept",
                 re.maxit = 1L, verbose = FALSE)@fits[[1]]
  si <- spiDE:::.absorbSpec(fi)
  # intercept-only: the nested indicators, 1x1 blocks, i.e. the logical path
  expect_true(is.logical(si))
  expect_equal(sum(si), sum(fi@re_group == "SampleCellTypeInt", na.rm = TRUE))
  # and the shared-factor batched solver gets that same logical
  expect_identical(spiDE:::.absorbBatchSpec(fi), si)

  fs <- .slopeWiringFit()$f
  ss <- spiDE:::.absorbSpec(fs)
  # slope: every random column, grouped by sample
  expect_false(is.logical(ss))
  expect_identical(is.na(ss), is.na(fs@re_group))
  expect_gt(length(unique(ss[!is.na(ss)])), 1L)
})

# ---------------------------------------------------------------------------
# The changed code, run end to end on a slope fit. The structural tests above
# say the grouping is there; these say the polish and the inference USE it and
# get the answer the dense path gets. Absorption is an exact rearrangement, so
# the standard is agreement to rounding, not "close".
# ---------------------------------------------------------------------------

# Put the random block FIRST. The design appends it after the fixed columns, so
# the fixed columns are a prefix of both the nested-only and the whole-block
# dense sets and the two `sel_x` indices coincide -- a test on the design as
# built cannot tell them apart. Moving the random block to the front makes the
# fixed columns' positions differ between the two, so indexing the CPU
# covariance with the wrong one picks the wrong entries.
.randomFirst <- function(f) {
  rnd <- which(!is.na(f@re_group))
  perm <- c(rnd, setdiff(seq_len(ncol(f@W)), rnd))
  f@W <- f@W[, perm, drop = FALSE]
  f@alpha <- f@alpha[, perm, drop = FALSE]
  f@covtype <- f@covtype[perm]
  f@coefmap <- f@coefmap[perm, , drop = FALSE]
  f@re_group <- f@re_group[perm]
  f@re_sample <- f@re_sample[perm]
  if (length(f@penalty) > 1L) f@penalty <- f@penalty[perm]
  f
}

test_that("the polish of a slope fit with the whole random block absorbed matches the dense polish", {
  d <- .slopeWiringFit()
  f <- d$f
  spec <- spiDE:::.absorbSpec(f)
  expect_false(is.logical(spec))                # the per-sample grouping path
  dense <- rep(FALSE, ncol(f@W))
  st <- spiDE:::.testedStartCols(f)
  for (eng in c("gene", "batch")) {
    ab <- SpaNorm::polishNB(d$Y, f@W, f@alpha, f@psi, lambda.a = d$pen,
                            absorb = spec,
                            absorb.batch = spiDE:::.absorbBatchSpec(f),
                            start.cols = st, engine = eng)
    dn <- SpaNorm::polishNB(d$Y, f@W, f@alpha, f@psi, lambda.a = d$pen,
                            absorb = dense, start.cols = st, engine = eng)
    expect_true(all(ab$polish$polished), label = eng)
    expect_equal(ab$alpha, dn$alpha, tolerance = 1e-8, label = eng)
    expect_equal(ab$psi, dn$psi, tolerance = 1e-8, label = eng)
  }
})

test_that("inference on a slope fit with the whole random block absorbed matches the dense covariance", {
  d <- .slopeWiringFit()
  for (f in list(d$f, .randomFirst(d$f))) {
    got <- spiDE:::.blockedInference(f, d$Y)
    # absorption forced OFF: with no column tagged as a nested indicator,
    # .blockedInference() builds no `absorb` and inverts the full dense gram
    fd <- f
    fd@re_group[!is.na(fd@re_group) & fd@re_group == "SampleCellTypeInt"] <- "NotAbsorbed"
    ref <- spiDE:::.blockedInference(fd, d$Y)
    expect_true(any(is.finite(got@t_stat)))
    expect_equal(got@se, ref@se, tolerance = 1e-8)
    expect_equal(got@t_stat, ref@t_stat, tolerance = 1e-8)
  }
})

test_that("the batched (shared-factor) polish of a slope fit runs and matches the per-gene solver", {
  # The shared-factor batched solver is the device path: SpaNorm::polishNB()
  # builds it only when a GPU is active. It absorbs 1x1 blocks only, so a slope
  # fit must hand it the nested indicators (.absorbBatchSpec()) and not the
  # per-sample grouping the per-gene CPU solver takes (.absorbSpec()). Reach
  # that path on the CPU by reporting a GPU and making the device transfer the
  # identity -- every batched kernel runs on a base matrix too (SpaNorm's
  # test-polishEngine.R).
  d <- .slopeWiringFit()
  f <- d$f
  spec <- spiDE:::.absorbSpec(f)
  spec_batch <- spiDE:::.absorbBatchSpec(f)
  # the batched solver's absorption is the nested block alone, inside the
  # grouping the per-gene solver absorbs whole
  expect_true(is.logical(spec_batch))
  expect_identical(spec_batch,
                   !is.na(f@re_group) & f@re_group == "SampleCellTypeInt")
  expect_true(all(!is.na(spec[spec_batch])))
  st <- spiDE:::.testedStartCols(f)
  cpu <- SpaNorm::polishNB(d$Y, f@W, f@alpha, f@psi, lambda.a = d$pen,
                           absorb = spec, absorb.batch = spec_batch,
                           start.cols = st, engine = "batch", backend = "cpu")
  testthat::local_mocked_bindings(
    checkGPU = function(...) TRUE,
    toGPUMatrix = function(x, ...) x,
    # SpaNorm::polishNB() calls the internal .requireFloat64() with no
    # arguments, whose default `dtype = getBackendDtype()` is evaluated at
    # that call (test-polish-backend.R relies on the same mechanism), so
    # mocking the exported getBackendDtype() reaches the real refusal check
    # instead of mocking the internal directly.
    getBackendDtype = function(...) "float64",
    .package = "SpaNorm"
  )
  dev <- SpaNorm::polishNB(d$Y, f@W, f@alpha, f@psi, lambda.a = d$pen,
                           absorb = spec, absorb.batch = spec_batch,
                           start.cols = st, engine = "batch", backend = "gpu")
  expect_true(all(dev$polish$polished))
  # a shared factorisation refreshes on a different schedule from the per-gene
  # one, so the two converge to the same optimum rather than along one path
  expect_equal(dev$alpha, cpu$alpha, tolerance = 1e-5)
  expect_equal(dev$psi, cpu$psi, tolerance = 1e-5)
})

test_that("a SpiDEFit whose re_sample and re_group lengths disagree is invalid", {
  # .absorbSpec() falls back to the nested-only absorption on a length
  # mismatch, so subsetting one slot without the other must error, not degrade
  f <- .slopeWiringFit()$f
  expect_true(validObject(f))
  bad <- f
  bad@re_sample <- bad@re_sample[-1]
  expect_error(validObject(bad), "re_sample")
  # absent re_sample (a fit saved before the slot existed) stays valid
  old <- f
  old@re_sample <- NULL
  expect_true(validObject(old))
})
