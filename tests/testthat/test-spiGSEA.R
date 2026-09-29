# spiGSEA() on the new engines (experimental): the competitive set test.

test_that("spiGSEA returns one row per usable (set, index, niche)", {
  spe <- .toySPE(n_samples = 12, n_per = 100, n_genes = 12, seed = 3)
  res <- spiDE(spe, condition = "condition", sigma = 30, index = "A", procedure = "all",
               verbose = FALSE)
  sets <- list(s1 = paste0("G", 1:4), s2 = paste0("G", 5:9), tiny = "G1")
  gs <- spiGSEA(res, sets, test = "pooled", min.size = 3)
  expect_true(all(c("set", "index", "niche", "test", "size", "t", "p", "q", "rho") %in% colnames(gs)))
  expect_false("tiny" %in% gs$set)
  expect_true(all(gs$p >= 0 & gs$p <= 1))
  expect_true(all(gs$rho >= 0 & gs$rho <= 1))
})
