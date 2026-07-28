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

test_that("mergeNiches records the group membership in metadata", {
  spe <- .toySPE()
  spe <- buildNiches(spe, sigma = 20)
  spe <- mergeNiches(spe, groups = list(AC = c("A", "C")), sigma = 20)

  gm <- S4Vectors::metadata(spe)$spiDE_niche_groups[["Niche20"]]
  expect_setequal(names(gm), c("AC", "B"))
  expect_setequal(gm[["AC"]], c("A", "C"))
  expect_equal(gm[["B"]], "B")
})

test_that("re-merging composes membership back to fine cell types", {
  spe <- .toySPE()
  spe <- buildNiches(spe, sigma = 20)
  spe <- mergeNiches(spe, groups = list(AC = c("A", "C")), sigma = 20)
  spe <- mergeNiches(spe, groups = list(ACB = c("AC", "B")), sigma = 20)

  gm <- S4Vectors::metadata(spe)$spiDE_niche_groups[["Niche20"]]
  expect_setequal(names(gm), "ACB")
  expect_setequal(gm[["ACB"]], c("A", "B", "C"))
})

test_that("a mergeNiches run makes nicheDesign drop member interactions", {
  spe <- .toySPE()
  spe <- buildNiches(spe, sigma = 20)
  spe <- mergeNiches(spe, groups = list(AC = c("A", "C")), sigma = 20)

  des <- nicheDesign(spe, condition = "condition", sigma = 20)
  cm <- des$coefmap
  ac <- cm[!is.na(cm$niche) & cm$niche == "AC", ]
  expect_false(any(ac$index %in% c("A", "C")))
  expect_true(any(ac$index == "B" & ac$type == "ResponseNiche"))
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
