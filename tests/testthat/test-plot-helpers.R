# The plotting layer's theme, colours, checkers and shared data builders (0.99.34).

spe <- buildNiches(.toySPE(n_samples = 16, n_per = 150, seed = 7), sigma = 30)
fit <- fitSpiDE(spe, condition = "condition", sigma = 30, index = "A", verbose = FALSE)
res <- testSpiDE(fit, procedure = "all")

test_that("theme_spiDE() is bhuvad_theme with ticks, and refuses a bad rl", {
  th <- theme_spiDE()
  expect_true(ggplot2::is_theme(th))
  expect_equal(th$panel.border$colour, "black")
  expect_true(inherits(th$panel.grid, "element_blank") || inherits(th$panel.grid, "ggplot2::element_blank"))
  expect_equal(th$legend.title$face, "italic")
  expect_error(theme_spiDE(-1), "rl")
})

test_that("spiDEColours() gives each role its colours; the magnitude ramp darkens", {
  expect_equal(unname(spiDEColours("direction")), c("#B8791A", "#5B4B9A"))
  expect_equal(unname(spiDEColours("condition")), c("#0A9396", "#AE2012"))
  expect_length(spiDEColours("magnitude", 9), 9)
  lum <- function(h) colSums(grDevices::col2rgb(h) * c(0.2126, 0.7152, 0.0722))
  expect_true(all(diff(lum(spiDEColours("magnitude", 9))) < 0))
  expect_length(spiDEColours("nuisance", 10), 10)
  expect_named(spiDEColours(), c("up", "down", "first", "second", "mid", "faint", "light", "grey"))
})

test_that("plotmath labels survive quotes and backslashes in names", {
  lab <- .itTriplet('a"b\\c', "T cell", "B|x")
  expect_silent(e <- .parseLabels(lab))
  expect_true(is.expression(e))
  expect_equal(.fmtq(c(0.0234, 2e-8, NA)), c("0.023", "2e-08", "NA"))
})

test_that("the slope weights are the pooled test's, over the whole gene family", {
  sw <- .slopeWeights(fit@index$A)
  vp <- sweep(fit@index$A$v_model, c(2, 3), fit@index$A$factor, "*")
  expect_equal(sw$v, vp)
  expect_equal(sw$tau, .dlTau2(fit@index$A$beta, vp))
})

test_that("cell-type means are the dense CP10k means, and niche folds pool merged members", {
  M <- .cellTypeMeans(spe, c("G1", "G2"), "counts", "cell_type")
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  lib <- colSums(Y)
  keep <- lib > 0
  ref <- sapply(c("A", "B", "C"), function(k) {
    i <- keep & spe$cell_type == k
    rowMeans(t(t(Y[c("G1", "G2"), i]) / lib[i]) * 1e4)
  })
  expect_equal(unname(M), unname(ref), tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(unname(attr(M, "ncells")[c("A", "B", "C")]),
               as.integer(table(spe$cell_type[keep])[c("A", "B", "C")]))
  sm <- mergeNiches(spe, groups = list(BC = c("B", "C")), sigma = 30)
  fm <- fitSpiDE(sm, sigma = 30, index = "A", verbose = FALSE)
  tb <- testSpiDE(fm)@table
  f <- .nicheFold(tb[tb$gene == "G1", ], sm, fm)
  n <- attr(M, "ncells")
  pooled <- (M["G1", "B"] * n[["B"]] + M["G1", "C"] * n[["C"]]) / (n[["B"]] + n[["C"]])
  expect_equal(f[1], (pooled + 0.01) / (M["G1", "A"] + 0.01), tolerance = 1e-10)
})

test_that("expression bands split each index type's genes into fifths", {
  b <- .expressionBands(fit)
  expect_equal(nlevels(b$band), 5L)
  expect_lte(diff(range(table(b$band[b$index == "A"]))), 1L)
})

test_that("index-cell data are the fit's cells, and a foreign spe is refused", {
  d <- .indexCellData(spe, fit, "A", "B", c("G1", "G2"))
  expect_equal(dim(d$counts), c(2L, length(d$patient)))
  expect_setequal(unique(d$patient), fit@index$A$patients)
  other <- spe[, spe$sample_id %in% c("S1", "S2")]
  other$sample_id <- paste0("X", other$sample_id)
  expect_error(.indexCellData(other, fit, "A", "B", "G1"), "the object the fit was made from")
})

test_that("the checkers name what is wrong", {
  expect_error(checkResults(fit, "f()"), "SpiDEResults")
  expect_silent(checkResults(fit, "f()", fit.ok = TRUE))
  expect_error(checkTriplet(fit, "Z"), "fitted index types: A")
  expect_error(checkTriplet(fit, "A", "Z"), "niches tested in A")
  expect_error(checkTriplet(fit, "A", "B", c("G1", "nope")), "not tested in A: nope")
  expect_error(.resolveTest(testSpiDE(fitSpiDE(spe, sigma = 30, index = "A", verbose = FALSE)), "condition"),
               "no condition-specific test")
})

test_that("the niche fold names genes that are not in spe", {
  tb <- res@table[res@table$gene == "G1", ]
  tb$gene[1] <- "nope"
  expect_error(.nicheFold(tb, spe, fit), "not in 'spe': nope")
})
