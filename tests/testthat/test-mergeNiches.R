test_that("mergeNiches sums grouped columns and carries leftovers", {
  spe <- .toySPE()
  spe <- buildNiches(spe, sigma = 20)
  before <- SingleCellExperiment::reducedDim(spe, "Niche20")

  spe <- mergeNiches(spe, groups = list(AC = c("A", "C")), sigma = 20)
  after <- SingleCellExperiment::reducedDim(spe, "Niche20")

  expect_setequal(colnames(after), c("AC", "B"))
  # merged column equals the sum of its members
  expect_equal(after[, "AC"], before[, "A"] + before[, "C"])
  # leftover B is unchanged
  expect_equal(after[, "B"], before[, "B"])
})

test_that("mergeNiches errors on a missing bandwidth", {
  spe <- .toySPE()
  spe <- buildNiches(spe, sigma = 20)
  expect_error(
    mergeNiches(spe, groups = list(AC = c("A", "C")), sigma = 99),
    "not found"
  )
})

test_that("mergeNiches with NULL sigma updates all niche reducedDims", {
  spe <- .toySPE()
  spe <- buildNiches(spe, sigma = c(10, 20))
  spe <- mergeNiches(spe, groups = list(AC = c("A", "C")))
  expect_setequal(
    colnames(SingleCellExperiment::reducedDim(spe, "Niche10")),
    c("AC", "B")
  )
  expect_setequal(
    colnames(SingleCellExperiment::reducedDim(spe, "Niche20")),
    c("AC", "B")
  )
})
