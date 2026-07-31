test_that("nicheDesign tags covariates and drops self-interactions", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  des <- nicheDesign(spe, condition = "condition", sigma = 20)

  # ResponseCellType is the cell-means CellType:condition term added with the
  # celltype-response design; the bare "Response" main effect no longer exists.
  expect_true(all(levels(des$covtype) %in%
    c("CellType", "Niche", "Response", "ResponseCellType", "ResponseNiche",
      "Other")))
  # no ResponseNiche/Niche column may have index == niche
  cm <- des$coefmap
  self <- !is.na(cm$index) & !is.na(cm$niche) & cm$index == cm$niche
  expect_false(any(self))
  # 3 cell types -> 3*3 - 3 self = 6 ResponseNiche columns
  expect_equal(sum(cm$type == "ResponseNiche"), 6)
  # design is full rank
  expect_equal(qr(des$W)$rank, ncol(des$W))
})

test_that(".isSelfNiche: no map reduces to exact index == niche", {
  idx <- c("A", "A", "B", NA)
  nch <- c("A", "B", "B", NA)
  expect_equal(spiDE:::.isSelfNiche(idx, nch), c(TRUE, FALSE, TRUE, FALSE))
})

test_that(".isSelfNiche: member index against a merged niche is dropped", {
  gm <- list(AC = c("A", "C"), B = "B")
  # A and C are members of AC -> drop; B is not -> keep; self B:B -> drop
  idx <- c("A", "C", "B", "B")
  nch <- c("AC", "AC", "AC", "B")
  expect_equal(spiDE:::.isSelfNiche(idx, nch, gm),
               c(TRUE, TRUE, FALSE, TRUE))
})

test_that(".isSelfNiche: sanitises raw member/key names before matching", {
  gm <- list(`T cell` = c("T cell CD4", "T cell CD8"))
  idx <- c("T.cell.CD4", "T.cell.CD8", "B.cell")
  nch <- c("T.cell", "T.cell", "T.cell")
  expect_equal(spiDE:::.isSelfNiche(idx, nch, gm), c(TRUE, TRUE, FALSE))
})

test_that("nicheDesign drops member interactions of a merged niche", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  # simulate a mergeNiches result: AC column = A + C, B carried through
  nm <- SingleCellExperiment::reducedDim(spe, "Niche20")
  merged <- cbind(AC = nm[, "A"] + nm[, "C"], B = nm[, "B"])
  rownames(merged) <- rownames(nm)
  SingleCellExperiment::reducedDim(spe, "Niche20") <- merged
  S4Vectors::metadata(spe)$spiDE_niche_groups <-
    list(Niche20 = list(AC = c("A", "C"), B = "B"))

  des <- nicheDesign(spe, condition = "condition", sigma = 20)
  cm <- des$coefmap
  ac <- cm[!is.na(cm$niche) & cm$niche == "AC", ]

  # no A- or C-indexed interaction survives against the AC niche
  expect_false(any(ac$index %in% c("A", "C")))
  # the genuine cross-type B:AC effect survives, as both Niche and ResponseNiche
  expect_true(any(ac$index == "B" & ac$type == "Niche"))
  expect_true(any(ac$index == "B" & ac$type == "ResponseNiche"))
  # self B:B is still dropped
  # CellType:condition rows carry NA niche, so the niche side needs guarding
  # too or `any()` returns NA rather than FALSE.
  expect_false(any(!is.na(cm$niche) & !is.na(cm$index) &
                     cm$niche == "B" & cm$index == "B"))
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
