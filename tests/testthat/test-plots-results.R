# The results plots (0.99.34): where the calls are, the strongest triplets,
# the spillover signature, triplet heatmaps and p-value bands.

spe <- buildNiches(.toySPE(n_samples = 16, n_per = 150, seed = 7), sigma = 30)
fit <- fitSpiDE(spe, condition = "condition", sigma = 30, verbose = FALSE)
res <- testSpiDE(fit, procedure = "all")
builds <- function(p) {
  expect_true(inherits(p, "ggplot"))
  expect_no_warning(b <- ggplot2::ggplot_build(p))
  invisible(b)
}

test_that("plotCallMap() counts the calls per pair and direction", {
  p <- plotCallMap(res, test = "pooled", fdr = 0.05)
  builds(p)
  tab <- res@table[res@table$test == "pooled" & !is.na(res@table$q) & res@table$q <= 0.05, ]
  expect_equal(sum(p$data$calls, na.rm = TRUE), nrow(tab))
  expect_true(all(c("index", "niche", "direction", "calls", "tested") %in% names(p$data)))
  expect_false(any(p$data$tested & p$data$index == p$data$niche))
  g <- plotCallMap(res, test = "condition", style = "graph")
  b <- builds(g)
  # node labels sit outside their nodes (the nodes are on the unit circle)
  txt <- b$data[[which(vapply(g$layers, function(l) inherits(l$geom, "GeomText"), logical(1)))]]
  expect_true(all(sqrt(txt$x^2 + txt$y^2) >= 1.25))
  expect_true(all(txt$vjust[txt$y > 0.5] == 0) && all(txt$vjust[txt$y < -0.5] == 1))
})

test_that("plotCallMap() draws when nothing is called", {
  builds(plotCallMap(res, test = "pooled", fdr = 1e-12))
  builds(plotCallMap(res, test = "pooled", fdr = 1e-12, style = "graph"))
})

test_that("plotTopTriplets() ranks, filters by gene, and flags niche markers with spe", {
  p <- plotTopTriplets(res, test = "pooled", n = 10)
  builds(p)
  expect_equal(nrow(p$data), 10L)
  expect_true(all(c("gene", "index", "niche", "estimate", "lo", "hi", "called", "direction", "marker",
                    "label") %in% names(p$data)))
  p <- plotTopTriplets(res, gene = c("G1", "G3"), spe = spe)
  builds(p)
  expect_setequal(unique(p$data$gene), c("G1", "G3"))
  expect_true(is.logical(p$data$marker))
  expect_error(plotTopTriplets(res, gene = "nope"), "no triplet")
})

test_that("plotSpillover() places every tested triplet by its niche fold", {
  p <- plotSpillover(res, spe)
  builds(p)
  n <- sum(res@table$test == "pooled" & is.finite(res@table$p))
  expect_equal(nrow(p$data), n)
  expect_true(all(c("fold", "z", "status") %in% names(p$data)))
})

test_that("merged niches: the plots accept the group, the fold pools its members", {
  sm <- mergeNiches(spe, groups = list(BC = c("B", "C")), sigma = 30)
  rmg <- testSpiDE(fitSpiDE(sm, condition = "condition", sigma = 30, index = "A", verbose = FALSE))
  builds(plotSpillover(rmg, sm))
  builds(plotTopTriplets(rmg, spe = sm, n = 5))
  builds(plotCallMap(rmg))
})

test_that("plotTripletHeatmap() draws genes, marks calls, greys untested cells", {
  p <- plotTripletHeatmap(res, test = "pooled", n = 8)
  builds(p)
  expect_true(all(c("feature", "index", "niche", "t", "called") %in% names(p$data)))
  expect_equal(length(unique(p$data$feature)), 8L)
  r2 <- res
  r2@table <- r2@table[!(r2@table$gene == "G2" & r2@table$index == "B"), ]
  p <- plotTripletHeatmap(r2, features = c("G1", "G2"), test = "pooled")
  b <- builds(p)
  expect_true(any(vapply(b$data, function(l) nrow(l) > 0 && "linetype" %in% names(l), logical(1))))
})

test_that("plotTripletHeatmap() draws spiGSEA() sets", {
  gs <- spiGSEA(res, list(one = paste0("G", 1:5), two = paste0("G", 6:12), three = paste0("G", 13:20)),
                min.size = 3)
  p <- plotTripletHeatmap(gs, test = "pooled")
  builds(p)
  expect_setequal(unique(as.character(p$data$feature)), c("one", "two", "three"))
  expect_error(plotTripletHeatmap(data.frame(a = 1)), "SpiDEResults or a spiGSEA")
})

test_that("plotPValues() draws each test in expression fifths", {
  p <- plotPValues(res)
  b <- builds(p)
  expect_true(all(c("p", "test", "band") %in% names(p$data)))
  expect_equal(nlevels(p$data$band), 5L)
  expect_setequal(as.character(unique(p$data$test)), c("pooled", "condition"))
})

test_that("the p-value and heatmap plots stop readably when nothing has a p-value", {
  r0 <- res
  r0@table$p <- NA_real_
  expect_error(plotPValues(r0), "no tested triplet")
  expect_error(plotTripletHeatmap(r0, test = "pooled"), "no tested triplet")
})

test_that("the call graph wraps long cell-type names", {
  r3 <- res
  r3@table$niche[r3@table$niche == "B"] <- "Airway epithelium B"
  r3@table$index[r3@table$index == "B"] <- "Airway epithelium B"
  g <- plotCallMap(r3, test = "pooled", style = "graph")
  b <- builds(g)
  txt <- b$data[[which(vapply(g$layers, function(l) inherits(l$geom, "GeomText"), logical(1)))]]
  expect_true(all(nchar(unlist(strsplit(txt$label, "\n"))) <= 12))
})

test_that("plots of a test without finite p-values stop readably", {
  r0 <- res
  r0@table$p[r0@table$test == "condition"] <- NA_real_
  r0@table$q[r0@table$test == "condition"] <- NA_real_
  expect_error(plotCallMap(r0), "no tested triplet")
  expect_error(plotSpillover(r0, spe, test = "condition"), "no tested triplet")
})
