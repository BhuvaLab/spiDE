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

# --- shared condition-free pipeline run ------------------------------------
# Two bandwidths, not three: this still exercises the cross-bandwidth Cauchy
# combination, but test-spiDE-e2e.R's three-bandwidth run already costs ~475 s
# of Bioconductor's 600 s budget on its own.
spe_niche <- .toyNiche()
res_niche <- spiDE(spe_niche, condition = NULL, sigma = c(10, 30),
                   covariates = "Age", verbose = FALSE)
tab_niche <- results(res_niche)

test_that(".toyNiche carries no condition column", {
  cd <- SummarizedExperiment::colData(spe_niche)
  expect_false("condition" %in% colnames(cd))
  expect_true(all(c("sample_id", "cell_type", "Age", "Area") %in% colnames(cd)))
})

test_that("spiDE runs condition-free and populates the results table", {
  expect_s4_class(res_niche, "SpiDEResults")
  expect_identical(res_niche@mode, "niche")
  expect_true(is.na(res_niche@condition))
  expect_identical(fits(res_niche)[[1]]@mode, "niche")

  expect_true(nrow(tab_niche) > 0)
  expect_true(all(c(
    "gene", "ct_index", "ct_niche", "bandwidth.max", "coef", "t",
    "DirectionGene", "DirectionIndex", "DirectionNiche",
    "fdr.gene", "fdr.index", "fdr.niche"
  ) %in% colnames(tab_niche)))
  expect_true(all(tab_niche$fdr.niche >= 0 & tab_niche$fdr.niche <= 1))
})

test_that("the planted condition-free niche effect is recovered and specific", {
  ab <- tab_niche[tab_niche$gene == "G1" & tab_niche$ct_index == "A" &
                  tab_niche$ct_niche == "B", ]
  expect_equal(nrow(ab), 1)
  expect_equal(ab$DirectionNiche, "Up")
  expect_gt(abs(ab$t), 5)

  # Specificity is NOT asserted by argmax here. Index A has exactly two
  # competing niches (B and C), so "B is strongest" reduces to |t_B| > |t_C|,
  # which is close to a coin flip at this effect size: measured over 10 seeds,
  # B wins 7/10 with median |t| 6.85 against C's 4.96. The signal is real but
  # the argmax is not a reliable single-seed test, and whether it passes
  # depends on which realisation the shipped fixture happens to be. The
  # direction and magnitude assertions above test recovery directly and do not
  # have that failure mode.
  g1a <- tab_niche[tab_niche$gene == "G1" & tab_niche$ct_index == "A", ]
  expect_true(all(c("B", "C") %in% g1a$ct_niche))
})

test_that("condition-free results carry no celltype or patient layer", {
  expect_equal(nrow(results(res_niche, type = "celltype")), 0L)
  expect_equal(nrow(results(res_niche, type = "patient")), 0L)
  expect_length(fits(res_niche)[[1]]@se_patient, 0L)
})

test_that("condition-free results resolve the index and niche cell type sets", {
  # These are read off the tested columns, which in niche mode are tagged
  # "Niche" -- a bare == "ResponseNiche" lookup is all-FALSE here and would
  # leave both slots empty.
  expect_setequal(res_niche@index, c("A", "B", "C"))
  expect_setequal(res_niche@niche, c("A", "B", "C"))
})

test_that("show() reports the niche-only mode and the cell type sets", {
  txt <- paste(capture.output(show(res_niche)), collapse = " ")
  expect_match(txt, "niche-only")
  expect_no_match(txt, "Condition: NA")
  expect_match(txt, "Index cell types: A, B, C")
  expect_match(txt, "Niche cell types: A, B, C")
})

test_that("checkSample accepts a condition-free call", {
  expect_true(spiDE:::checkSample(spe_niche, condition = NULL))
  # the sample-constant-covariate rejection still applies: Age is constant
  # within sample, so the per-sample random intercept would absorb it
  expect_error(
    spiDE:::checkSample(spe_niche, condition = NULL, covariates = "Age"),
    "constant within sample")
})

# --- gene-set layer on a condition-free fit --------------------------------

test_that("spiGSEA type='niche' works on a condition-free fit", {
  sets <- list(setA = rownames(spe_niche)[1:5],
               setB = rownames(spe_niche)[6:10])
  gs <- spiGSEA(res_niche, spe_niche, sets, min.size = 3, fdr = 1,
                verbose = FALSE)
  expect_s3_class(gs, "data.frame")
  expect_true(nrow(gs) > 0)
  expect_true(all(c("geneset", "ct_index", "ct_niche") %in% colnames(gs)))
})

test_that("spiGSEA type='celltype' is refused on a condition-free fit", {
  sets <- list(setA = rownames(spe_niche)[1:5])
  expect_error(
    spiGSEA(res_niche, spe_niche, sets, type = "celltype", min.size = 3,
            fdr = 1, verbose = FALSE),
    "ResponseCellType")
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
