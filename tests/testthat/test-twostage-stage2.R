test_that("SpiDEResults carries a diagnostics list slot", {
  r <- new("SpiDEResults", fits = list(), sigma = 30, condition = "condition",
           mode = "condition", index = "A", niche = "B",
           covariates = character(), coldata = S4Vectors::DataFrame(),
           gene.weights = NULL, p.cauchy.pos = NULL, p.cauchy.neg = NULL,
           results = data.frame(), fdr = 0.05, call = NULL,
           diagnostics = list(r2 = data.frame()))
  expect_identical(names(r@diagnostics), "r2")
  # default construction still works and defaults to an empty list
  r0 <- new("SpiDEResults", fits = list(), sigma = 30, condition = "c",
            mode = "condition", index = "A", niche = "B",
            covariates = character(), coldata = S4Vectors::DataFrame(),
            gene.weights = NULL, p.cauchy.pos = NULL, p.cauchy.neg = NULL,
            results = data.frame(), fdr = 0.05, call = NULL)
  expect_identical(r0@diagnostics, list())
})
