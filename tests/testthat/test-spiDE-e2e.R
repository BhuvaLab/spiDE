# End to end: the toy fixture plants G1 up in index type A, in Responders, in
# proportion to the local density of B. With enough patients both engines
# find it, in the pooled test and in the condition-specific test.

spe <- .toySPE(n_samples = 16, n_per = 150, n_genes = 12, seed = 7)

test_that("the slopes engine recovers the planted G1 / A / B effect", {
  res <- spiDE(spe, condition = "condition", sigma = 30, index = "A", procedure = "all",
               verbose = FALSE)
  cond <- results(res, test = "condition")
  pool <- results(res, test = "pooled")
  top <- cond[cond$gene == "G1", ]
  expect_equal(top$niche[which.min(top$p)], "B")
  expect_lt(min(top$q), 0.05)
  expect_gt(top$estimate[top$niche == "B"], 0)
  # planted in Responders only, so the pooled slope is about half the size: it
  # is positive, but the condition-specific test is the one that finds it
  expect_gt(pool$estimate[pool$gene == "G1" & pool$niche == "B"], 0)
})

test_that("the sandwich engine recovers it too", {
  res <- spiDE(spe, condition = "condition", sigma = 30, index = "A", engine = "sandwich",
               procedure = "all", verbose = FALSE)
  top <- results(res, test = "condition")
  top <- top[top$gene == "G1", ]
  expect_equal(top$niche[which.min(top$p)], "B")
  expect_lt(min(top$q), 0.05)
})

test_that("the sandwich engine refuses a condition it was not fitted with", {
  sp <- buildNiches(spe, sigma = 30, verbose = FALSE)
  fs <- fitSpiDE(sp, engine = "sandwich", index = "A", sigma = 30, verbose = FALSE)
  expect_error(testSpiDE(fs, condition = "condition"), "fitted with")
  expect_identical(unique(testSpiDE(fs)@table$test), "pooled")
})
