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
  res <- fitSpiDE(spe, condition = "condition", sigma = 20, random = "none", verbose = FALSE)
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

test_that(".nicheRecords returns two-sided p-values", {
  # ACAT needs U(0,1) input; one-sided p on (0, 0.5) makes tan((0.5-p)*pi)
  # strictly positive, so the cross-bandwidth combination is a half-Cauchy.
  spe <- buildNiches(.toySPE(), sigma = 20)
  res <- testSpiDE(fitSpiDE(spe, condition = "condition", sigma = 20,
                            random = "none", verbose = FALSE), spe = spe, verbose = FALSE)
  fitl <- fits(res)
  genes <- rownames(fitl[[1]]@alpha)[1:5]
  recs <- spiDE:::.nicheRecords(fitl, genes)
  expect_true(max(recs$p) > 0.5)
  expect_equal(recs$p, 2 * pnorm(-abs(recs$t)), tolerance = 1e-10)
})

test_that(".nicheLevelFDR combines bandwidths on the two-sided scale", {
  spe <- buildNiches(.toySPE(), sigma = c(20, 40))
  res <- testSpiDE(fitSpiDE(spe, condition = "condition", sigma = c(20, 40),
                            random = "none", verbose = FALSE), spe = spe, verbose = FALSE)
  fitl <- fits(res)
  gene.w <- spiDE:::.geneWeights(fitl)
  genes <- rownames(fitl[[1]]@alpha)[1:5]
  recs <- spiDE:::.nicheRecords(fitl, genes)
  gi <- unique(recs[, c("gene", "ct_index")])

  out <- spiDE:::.nicheLevelFDR(fitl, gi, gene.w, fdr = 1)

  # BH at the largest rank of a group returns that group's largest p, so the
  # largest fdr.niche IS the largest combined p -- an exact anchor.
  sigmas <- vapply(fitl, function(f) f@sigma, numeric(1))
  trip <- paste(recs$gene, recs$ct_index, recs$ct_niche, sep = "\r")
  ref <- vapply(split(seq_len(nrow(recs)), trip), function(ix) {
    s <- recs[ix, , drop = FALSE]
    p2 <- 2 * pnorm(-abs(s$t))
    w <- gene.w[s$gene[1], match(s$bandwidth, sigmas)]
    spiDE:::.cauchyCombine(matrix(p2, nrow = 1), matrix(w, nrow = 1))
  }, numeric(1))
  key <- paste(out$gene, out$ct_index, sep = "\r")
  got <- vapply(split(out$fdr.niche, key), max, numeric(1))
  want <- vapply(split(ref, sub("\r[^\r]*$", "", names(ref))), max, numeric(1))
  expect_equal(got[names(want)], want, tolerance = 1e-8)
})
