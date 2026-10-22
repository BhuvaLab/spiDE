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

test_that("a merge after a merge composes the members at every bandwidth, with one metadata record", {
  spe <- buildNiches(.toySPE(), sigma = c(10, 20), verbose = FALSE)
  fine <- SingleCellExperiment::reducedDim(spe, "Niche20")
  m1 <- mergeNiches(spe, groups = list(AC = c("A", "C")))
  m2 <- mergeNiches(m1, groups = list(ACB = c("AC", "B")))
  for (rd in c("Niche10", "Niche20")) {
    gm <- S4Vectors::metadata(m2)$spiDE_niche_groups[[rd]]
    expect_identical(names(gm), "ACB")
    expect_setequal(gm$ACB, c("A", "B", "C"))
  }
  expect_equal(unname(SingleCellExperiment::reducedDim(m2, "Niche20")[, "ACB"]), unname(rowSums(fine)))
  # a second merge of one bandwidth leaves the other bandwidth's record alone
  m3 <- mergeNiches(m1, groups = list(ACB = c("AC", "B")), sigma = 20)
  gm <- S4Vectors::metadata(m3)$spiDE_niche_groups
  expect_setequal(names(gm), c("Niche10", "Niche20"))
  expect_setequal(gm$Niche20$ACB, c("A", "B", "C"))
  expect_setequal(gm$Niche10$AC, c("A", "C"))
})
