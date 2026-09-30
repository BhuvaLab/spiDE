# The niche plots (0.99.34): the niche covariate in tissue, and the index
# cells' expression against it.

spe <- buildNiches(.toySPE(n_samples = 16, n_per = 150, seed = 7), sigma = 30)
fit <- fitSpiDE(spe, condition = "condition", sigma = 30, index = "A", verbose = FALSE)
res <- testSpiDE(fit, procedure = "all")
builds <- function(p) {
  expect_true(inherits(p, "ggplot"))
  expect_no_warning(b <- ggplot2::ggplot_build(p))
  invisible(b)
}

test_that("plotNicheMap() colours one sample's cells by a niche's density", {
  p <- plotNicheMap(spe, niche = "B", sample = "S1")
  builds(p)
  expect_true(all(c("x", "y", "cell_type", "density", "section", "focal") %in% names(p$data)))
  expect_equal(nrow(p$data), sum(spe$sample_id == "S1"))
  expect_equal(p$data$density, log1p(SingleCellExperiment::reducedDim(spe, "Niche30")[spe$sample_id == "S1", "B"]),
               ignore_attr = TRUE)
  p <- plotNicheMap(spe, niche = "B", sample = "S1", index = "A", unit = "µm")
  builds(p)
  expect_equal(sum(p$data$focal), sum(spe$sample_id == "S1" & spe$cell_type == "A"))
  # the scale bar's label sits below the bar: the panel must leave room for it
  b <- ggplot2::ggplot_build(p)
  yr <- range(p$data$y)
  expect_lt(b$layout$panel_params[[1]]$y.range[1], yr[1] - 0.1 * diff(yr))
  expect_error(plotNicheMap(spe, niche = "Z", sample = "S1"), "not a column")
  expect_error(plotNicheMap(spe, niche = "B", sample = "nope"), "nope")
})

test_that("plotNicheResponse() adjusted: one panel per gene, patients' bins and group means", {
  p <- plotNicheResponse(spe, res, gene = c("G1", "G2"), index = "A", niche = "B")
  b <- builds(p)
  expect_true(all(c("gene", "patient", "bin", "expr", "density", "condition", "gene_lab") %in% names(p$data)))
  expect_equal(length(unique(b$layout$layout$gene_lab)), 2L)
  expect_true(all(is.finite(p$data$expr)))
  # thousands marked with a comma, not the scales default of a space
  labs <- unlist(lapply(b$layout$panel_params, function(pp) pp$y$get_labels()))
  expect_false(any(grepl("[0-9] [0-9]", labs)))
})

test_that("plotNicheResponse() raw, the sandwich engine, and a fit without intercepts", {
  builds(plotNicheResponse(spe, res, gene = "G1", index = "A", niche = "B", adjust = FALSE))
  fs <- fitSpiDE(spe, condition = "condition", sigma = 30, index = "A", engine = "sandwich", verbose = FALSE)
  builds(plotNicheResponse(spe, testSpiDE(fs), gene = "G1", index = "A", niche = "B"))
  f <- fit
  f@index$A$intercept <- NULL
  expect_error(plotNicheResponse(spe, f, "G1", "A", "B"), "refit")
  builds(plotNicheResponse(spe, f, "G1", "A", "B", adjust = FALSE))
})

test_that("plotNicheResponse() refuses an spe that is not the fit's", {
  other <- spe[, spe$sample_id %in% c("S1", "S2")]
  expect_error(plotNicheResponse(other, res, "G1", "A", "B"), "the object the fit was made from")
})

test_that("plotNicheResponse() drops a gene the fit could not fit, and stops if none is left", {
  f <- fit
  f@index$A$intercept["G2", ] <- NA
  f@index$A$beta["G2", , ] <- NA
  expect_message(p <- plotNicheResponse(spe, f, c("G1", "G2"), "A", "B"), "G2")
  builds(p)
  expect_setequal(unique(p$data$gene), "G1")
  expect_error(suppressMessages(plotNicheResponse(spe, f, "G2", "A", "B")), "no gene")
})

test_that("plotNicheResponse() references the patients' median intercept, robust to one extreme", {
  p1 <- plotNicheResponse(spe, fit, "G1", "A", "B")
  f <- fit
  f@index$A$intercept["G1", "S2"] <- -50
  p2 <- plotNicheResponse(spe, f, "G1", "A", "B")
  s1 <- function(p) p$data$expr[p$data$patient == "S1"]
  expect_lt(max(abs(log(s1(p2) / s1(p1)))), 0.5)
})

test_that("plotNicheMap() puts each section of a patient on its own origin", {
  s2 <- spe
  s2$section <- ifelse(SpatialExperiment::spatialCoords(s2)[, 1] > 250, "right", "left")
  p <- plotNicheMap(s2, niche = "B", sample = "S1", section = "section")
  builds(p)
  expect_setequal(unique(p$data$section), c("left", "right"))
  expect_equal(as.numeric(tapply(p$data$x, p$data$section, min)), c(0, 0))
  expect_equal(as.numeric(tapply(p$data$y, p$data$section, min)), c(0, 0))
})
