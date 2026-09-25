# Wiring the per-sample absorption grouping from the design through to the
# polish and the inference.
#
# .newtonSolver() can absorb a random block that is block-diagonal by sample,
# but nothing reaches that path while the call sites select the absorbed columns
# with `re_group == "SampleCellTypeInt"`. That literal is also exactly what the
# architecture notes warn against: it silently skips a block that a new mode
# adds. The design knows which sample each random column belongs to, so it
# should say so, and the call sites should ask.

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
  # the property .newtonSolver() relies on: no cell loads on two blocks
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
  spe <- buildNiches(spiDE:::.toySPE(n_genes = 6, n_per = 40), sigma = 30)
  fit <- fitSpiDE(spe, condition = "condition", sigma = 30, random = "slope",
                  re.maxit = 1L, verbose = FALSE)
  f1 <- fit@fits[[1]]
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

  fs <- fitSpiDE(spe, condition = "condition", sigma = 30, random = "slope",
                 re.maxit = 1L, verbose = FALSE)@fits[[1]]
  ss <- spiDE:::.absorbSpec(fs)
  # slope: every random column, grouped by sample
  expect_false(is.logical(ss))
  expect_identical(is.na(ss), is.na(fs@re_group))
  expect_gt(length(unique(ss[!is.na(ss)])), 1L)
})
