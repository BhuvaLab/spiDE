test_that("nicheDesign tags covariates and drops self-interactions", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  des <- nicheDesign(spe, condition = "condition", sigma = 20)

  expect_true(all(levels(des$covtype) %in%
    c("CellType", "Niche", "Response", "ResponseNiche", "Other")))
  # no ResponseNiche/Niche column may have index == niche
  cm <- des$coefmap
  self <- !is.na(cm$index) & !is.na(cm$niche) & cm$index == cm$niche
  expect_false(any(self))
  # 3 cell types -> 3*3 - 3 self = 6 ResponseNiche columns
  expect_equal(sum(cm$type == "ResponseNiche"), 6)
  # design is full rank
  expect_equal(qr(des$W)$rank, ncol(des$W))
})

test_that("restricting niche cell types reduces interaction columns", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  des <- nicheDesign(spe, condition = "condition", sigma = 20, niche = "B")
  cm <- des$coefmap
  # only B as a niche -> index in {A,C} -> 2 ResponseNiche columns
  expect_setequal(cm$niche[cm$type == "ResponseNiche"], "B")
  expect_equal(sum(cm$type == "ResponseNiche"), 2)
})

test_that("restricting index cell types reduces interaction columns", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  des <- nicheDesign(spe, condition = "condition", sigma = 20, index = "A")
  cm <- des$coefmap
  expect_setequal(cm$index[cm$type == "ResponseNiche"], "A")
})

test_that("nicheDesign errors when niches are not built", {
  spe <- .toySPE()
  expect_error(nicheDesign(spe, condition = "condition", sigma = 20), "buildNiches")
})

test_that("nicheDesign errors on a non-two-level condition", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  spe$cell_type[1] <- "A" # unchanged; add a bogus 3-level condition
  spe$bad <- rep(c("x", "y", "z"), length.out = ncol(spe))
  expect_error(nicheDesign(spe, condition = "bad", sigma = 20), "two levels")
})
