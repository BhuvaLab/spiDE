# Every plot draws on a PDF device (0.99.34): its text must be Latin-1 or
# plotmath, since the PDF device cannot convert symbols such as a Greek Delta
# or <= from UTF-8 (R CMD check runs the examples on one).

spe <- buildNiches(.toySPE(n_samples = 16, n_per = 150, seed = 7), sigma = 30)
res <- testSpiDE(fitSpiDE(spe, condition = "condition", sigma = 30, verbose = FALSE), procedure = "all")
draws <- function(p) {
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({ grDevices::dev.off(); unlink(f) })
  expect_no_warning(expect_no_error(print(p)))
}

test_that("every plot draws on a PDF device", {
  gs <- spiGSEA(res, list(one = paste0("G", 1:5), two = paste0("G", 6:12), three = paste0("G", 13:20)),
                min.size = 3)
  draws(plotPatientSlopes(res, gene = c("G1", "G2"), index = "A", niche = "B"))
  draws(plotPatientEffects(res, index = "A"))
  draws(plotNicheMap(spe, niche = "B", sample = "S1", index = "A", unit = "µm"))
  draws(plotNicheResponse(spe, res, gene = "G1", index = "A", niche = "B"))
  draws(plotCallMap(res, test = "pooled"))
  draws(plotCallMap(res, test = "pooled", style = "graph"))
  draws(plotTopTriplets(res, test = "pooled", n = 10, spe = spe))
  draws(plotSpillover(res, spe))
  draws(plotTripletHeatmap(res, test = "pooled", n = 8))
  draws(plotTripletHeatmap(gs, test = "pooled"))
  draws(plotPValues(res))
})
