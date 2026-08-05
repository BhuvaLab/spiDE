# Condition-free (niche-only) designs: condition = NULL promotes the two-way
# CellType:niche interactions from nuisance terms to the tested set. See
# specs/2026-08-06-niche-only-design-design.md.

test_that(".tagCovtype with response_coef = NULL yields no Response tags", {
  cols <- c("CellTypeA", "CellTypeA:B", "B", "Age")
  tg <- spiDE:::.tagCovtype(cols, niche_cols = c("A", "B", "C"),
                            response_coef = NULL)
  expect_equal(tg$type, c("CellType", "Niche", "Other", "Other"))
  expect_equal(tg$index, c(NA, "A", NA, NA))
  expect_equal(tg$niche, c(NA, "B", NA, NA))
})

test_that("nicheDesign(condition = NULL) builds a niche-only design", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  des <- nicheDesign(spe, condition = NULL, sigma = 20)

  expect_identical(des$mode, "niche")
  cm <- des$coefmap
  # no response terms of any kind survive
  expect_equal(sum(grepl("^Response", cm$type)), 0L)
  # 3 cell types -> 3 CellType main effects, 3*3 - 3 self = 6 CellType:niche
  expect_equal(sum(cm$type == "CellType"), 3L)
  expect_equal(sum(cm$type == "Niche"), 6L)
  # the self-niche drop still applies
  expect_false(any(!is.na(cm$index) & !is.na(cm$niche) & cm$index == cm$niche))
  # and the design is still full rank
  expect_equal(qr(des$W)$rank, ncol(des$W))
})

test_that("nicheDesign keeps condition mode when a condition is supplied", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  des <- nicheDesign(spe, condition = "condition", sigma = 20)
  expect_identical(des$mode, "condition")
  expect_equal(sum(des$coefmap$type == "ResponseNiche"), 6L)
})

test_that("condition-free designs drop merged-niche member interactions", {
  # .isSelfNiche() is mode-independent, but the drop must still be reached from
  # the condition-free branch -- the spec lists merged niches explicitly.
  spe <- buildNiches(.toySPE(), sigma = 20)
  nm <- SingleCellExperiment::reducedDim(spe, "Niche20")
  merged <- cbind(AC = nm[, "A"] + nm[, "C"], B = nm[, "B"])
  SingleCellExperiment::reducedDim(spe, "Niche20") <- merged
  S4Vectors::metadata(spe)[["spiDE_niche_groups"]] <-
    list(Niche20 = list(AC = c("A", "C"), B = "B"))

  des <- nicheDesign(spe, condition = NULL, sigma = 20)
  cm <- des$coefmap
  nc <- cm[cm$type == "Niche", ]
  # A and C are members of AC -> both dropped against it; B:B dropped too.
  # Survivors: A:B, C:B, B:AC
  expect_equal(nrow(nc), 3L)
  expect_false(any(nc$index %in% c("A", "C") & nc$niche == "AC"))
  expect_false(any(nc$index == "B" & nc$niche == "B"))
})

test_that("condition-free designs honour index and niche restrictions", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  des <- nicheDesign(spe, condition = NULL, sigma = 20, index = "A")
  cm <- des$coefmap
  expect_equal(sum(cm$type == "Niche"), 2L)      # A:B and A:C
  expect_true(all(cm$index[cm$type == "Niche"] == "A"))
})

test_that("the niche-only design carries no bare niche main effects", {
  # niche_n = sum_c CellType_c:niche_n, so keeping the main effects would make
  # model.matrix() alias away one interaction per niche -- always the first
  # cell type's -- and that cell type could then never be tested.
  spe <- buildNiches(.toySPE(), sigma = 20)
  des <- nicheDesign(spe, condition = NULL, sigma = 20)
  expect_false(any(colnames(des$W) %in% c("A", "B", "C")))
  # every index cell type keeps a slope against every non-self niche
  cm <- des$coefmap[des$coefmap$type == "Niche", ]
  expect_setequal(cm$index, c("A", "B", "C"))
  expect_equal(nrow(cm), 6L)
})

test_that(".testedCols / .nicheTestCols resolve the tested set by mode", {
  ct <- factor(
    c("CellType", "Niche", "ResponseNiche", "ResponseCellType", "Other",
      "Random"),
    levels = c("CellType", "Niche", "Response", "ResponseNiche",
               "ResponseCellType", "Other", "Random"))

  expect_equal(spiDE:::.testedCols(ct, "condition"),
               c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE))
  expect_equal(spiDE:::.nicheTestCols(ct, "condition"),
               c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
  expect_equal(spiDE:::.testedCols(ct, "niche"),
               c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_equal(spiDE:::.nicheTestCols(ct, "niche"),
               c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE))
})

test_that("condition-mode .testedCols reproduces the historical predicate", {
  # The code being replaced is grepl("Response", covtype), which also admits
  # the (currently unpopulated) bare "Response" level. Pin that equivalence so
  # the refactor in the next task cannot silently narrow the tested set.
  ct <- c("CellType", "Niche", "Response", "ResponseNiche", "ResponseCellType",
          "Other", "Random")
  expect_equal(spiDE:::.testedCols(ct, "condition"), grepl("Response", ct))
})

test_that(".fitMode falls back to condition for objects lacking the slot", {
  expect_identical(spiDE:::.fitMode(structure(list(), class = "foo")),
                   "condition")
})

test_that("mode defaults to condition and survives updateObject", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  res <- fitSpiDE(spe, condition = "condition", sigma = 20, verbose = FALSE)
  expect_identical(res@mode, "condition")
  expect_identical(fits(res)[[1]]@mode, "condition")

  # emulate an object serialised before @mode existed
  f <- fits(res)[[1]]
  attr(f, "mode") <- NULL
  expect_false("mode" %in% names(attributes(f)))
  expect_identical(updateObject(f)@mode, "condition")
})
