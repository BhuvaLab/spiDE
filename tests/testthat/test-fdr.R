test_that(".geneIndexFDR gates genes then index cell types", {
  genes <- paste0("G", 1:5)
  idx <- c("A", "B")
  # G1 strongly significant at index A; others null
  p.pos <- rbind(
    G1 = c(Gene = 1e-6, A = 1e-6, B = 0.9),
    G2 = c(0.9, 0.9, 0.9),
    G3 = c(0.9, 0.9, 0.9),
    G4 = c(0.9, 0.9, 0.9),
    G5 = c(0.9, 0.9, 0.9)
  )
  colnames(p.pos) <- c("Gene", idx)
  p.neg <- matrix(0.9, 5, 3, dimnames = list(genes, c("Gene", idx)))

  gi <- spiDE:::.geneIndexFDR(p.pos, p.neg, fdr = 0.05)
  expect_equal(gi$gene, "G1")
  expect_equal(gi$ct_index, "A")
  expect_equal(gi$DirectionGene, "Up")
})

test_that(".geneIndexFDR marks a two-directional hit as Both", {
  idx <- c("A", "B")
  sig <- c(Gene = 1e-8, A = 1e-8, B = 0.9)
  null <- c(0.9, 0.9, 0.9)
  p.pos <- rbind(G1 = sig, G2 = null)
  p.neg <- rbind(G1 = sig, G2 = null)
  colnames(p.pos) <- colnames(p.neg) <- c("Gene", idx)

  gi <- spiDE:::.geneIndexFDR(p.pos, p.neg, fdr = 0.05)
  expect_equal(gi$DirectionGene[gi$gene == "G1"], "Both")
})

test_that(".geneIndexFDR returns NULL when nothing passes", {
  p <- matrix(0.9, 3, 3, dimnames = list(paste0("G", 1:3), c("Gene", "A", "B")))
  expect_null(spiDE:::.geneIndexFDR(p, p, fdr = 0.05))
})

test_that(".hierarchicalFDR returns the empty schema when nothing passes", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  res <- fitSpiDE(spe, condition = "condition", sigma = 20, verbose = FALSE)
  fitl <- fits(res)
  gene.w <- spiDE:::.geneWeights(fitl)
  # all-null p-values (gene, plus the three index cell types)
  genes <- rownames(fitl[[1]]@alpha)
  p <- matrix(0.9, length(genes), 4,
    dimnames = list(genes, c("Gene", "A", "B", "C")))
  out <- spiDE:::.hierarchicalFDR(fitl, p, p, gene.w, fdr = 0.05)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0)
  expect_true(all(c("gene", "ct_index", "ct_niche", "fdr.niche") %in% names(out)))
})
