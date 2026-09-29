# depth = "spatial_spline": a library-size spline within each patient,
# absorbed with the patient's intercept. The inference residualises on each
# patient's block (Frisch-Waugh-Lovell); these tests check it against the
# full design.

spe_d <- buildNiches(.toySPE(n_samples = 8, n_per = 650, n_genes = 4, seed = 11), sigma = 30,
                     verbose = FALSE)
# fitSpiDE() leaves out cells with no counts (no library size); so do these tests
lib_d <- Matrix::colSums(SummarizedExperiment::assay(spe_d, "counts"))
ct_d <- ifelse(lib_d > 0, as.character(spe_d$cell_type), NA_character_)

test_that("depth blocks: df by section size, padded compactly, R2 diagnostic", {
  ct <- ct_d; smp <- as.character(spe_d$sample_id)
  ik <- spiDE:::.indexCells(ct, smp, "A", 10L)
  NM <- SingleCellExperiment::reducedDim(spe_d, "Niche30")
  L <- log1p(NM[ik, c("B", "C")])
  ell <- log(Matrix::colSums(SummarizedExperiment::assay(spe_d, "counts")))[ik]
  pat <- factor(smp[ik])
  blk <- spiDE:::.depthBlocks(ell, SpatialExperiment::spatialCoords(spe_d), smp, ik, pat, L)
  n <- as.numeric(table(pat))
  # one section per patient here: intercept + l + l * B (df^2 columns)
  expect_equal(vapply(blk$Zs, ncol, integer(1)), 1L + ifelse(n >= 200, 10L, ifelse(n >= 100, 5L, 1L)))
  expect_equal(ncol(blk$Zc), max(vapply(blk$Zs, ncol, integer(1))))
  for (s in seq_len(nlevels(pat))) {
    i <- which(as.integer(pat) == s)
    w <- ncol(blk$Zs[[s]])
    expect_equal(unname(blk$Zc[i, seq_len(w)]), unname(blk$Zs[[s]]))
    if (w < ncol(blk$Zc)) expect_true(all(blk$Zc[i, -seq_len(w)] == 0))
    expect_true(all(blk$Zs[[s]][, 1] == 1))
  }
  expect_true(all(blk$r2 >= 0 & blk$r2 <= 1, na.rm = TRUE))
})

# With multi-column blocks nested in patients, clubSandwich's own CR2 on the
# absorbed design and on the full design differ (SE by up to 2.4% on this
# fixture, 2026-09-30), unlike the intercept-only case where they coincide; the
# engine's absorbed CR2 lies between the two. The Bell-McCaffrey df agree to
# two decimals. So the SE is held to 3% of the full-design answer, the df to 0.5%.
test_that("spatial_spline: the sandwich engine's CR2 and df agree with clubSandwich on the full design", {
  skip_if_not_installed("clubSandwich")
  ct <- ct_d; smp <- as.character(spe_d$sample_id)
  ik <- spiDE:::.indexCells(ct, smp, "A", 10L)
  NM <- SingleCellExperiment::reducedDim(spe_d, "Niche30")
  nc <- spiDE:::.nicheColumns(NM, "A", ik)
  L <- log1p(NM[ik, nc$cols, drop = FALSE])
  ell <- log(Matrix::colSums(SummarizedExperiment::assay(spe_d, "counts")))[ik]
  pat <- factor(smp[ik])
  trt <- as.numeric(spe_d$condition[ik] == "Responder")
  blk <- spiDE:::.depthBlocks(ell, SpatialExperiment::spatialCoords(spe_d), smp, ik, pat, L)
  des <- spiDE:::.indexDesign(L, NULL, pat, trt = trt, tested = nc$tested, blocks = blk)
  Y <- SummarizedExperiment::assay(spe_d, "counts")[c("G1", "G2"), ik, drop = FALSE]
  fit <- spiDE:::.fitIndexGLM(Y, des)
  got <- spiDE:::.sandwichCR2(fit, des, Y)
  # the full design: each patient's own block columns, zero elsewhere, then the dense columns
  pid <- as.integer(pat)
  Zf <- do.call(cbind, lapply(seq_along(des$Zs), function(s) {
    z <- matrix(0, length(pid), ncol(des$Zs[[s]])); z[pid == s, ] <- des$Zs[[s]]; z }))
  Xd <- des$W$X
  colnames(Xd) <- make.names(colnames(Xd))
  for (g in rownames(Y)) {
    eta <- spiDE:::.linPred(des$W, fit$alpha[g, ]); mu <- exp(eta)
    w <- mu / (1 + fit$psi[g] * mu); z <- eta + (as.numeric(Y[g, ]) - mu) / mu
    m <- stats::lm(z ~ 0 + Zf + Xd, weights = w)
    V <- clubSandwich::vcovCR(m, cluster = pat, type = "CR2", inverse_var = TRUE)
    cs <- clubSandwich::coef_test(m, vcov = V, test = "Satterthwaite")
    o <- cs[match(paste0("Xd", make.names(colnames(des$W)[des$tested])), cs$Coef), ]
    mine <- got[got$gene == g, ]
    expect_equal(mine$se, unname(o$SE), tolerance = 0.03)
    expect_equal(mine$df, unname(o$df_Satt), tolerance = 0.005)
  }
})

test_that("spatial_spline: a patient slope is one Fisher step with the block partialled out", {
  fit <- fitSpiDE(spe_d, index = "A", sigma = 30, depth = "spatial_spline", verbose = FALSE)
  x <- fit@index$A
  expect_true(all(x$depth_r2 >= 0 & x$depth_r2 <= 1, na.rm = TRUE))
  ct <- ct_d; smp <- as.character(spe_d$sample_id)
  ik <- spiDE:::.indexCells(ct, smp, "A", 10L)
  NM <- SingleCellExperiment::reducedDim(spe_d, "Niche30")
  nc <- spiDE:::.nicheColumns(NM, "A", ik)
  L <- log1p(NM[ik, c(nc$tested, setdiff(nc$cols, nc$tested)), drop = FALSE])
  ell <- log(Matrix::colSums(SummarizedExperiment::assay(spe_d, "counts")))[ik]
  pat <- factor(smp[ik])
  blk <- spiDE:::.depthBlocks(ell, SpatialExperiment::spatialCoords(spe_d), smp, ik, pat, L)
  des <- spiDE:::.indexDesign(L, NULL, pat, tested = nc$tested, blocks = blk)
  Y <- SummarizedExperiment::assay(spe_d, "counts")[x$genes, ik, drop = FALSE]
  f <- spiDE:::.fitIndexGLM(Y, des)
  g <- 1L; s <- 2L
  i <- which(des$patient == s)
  wr <- spiDE:::.workingWR(as.numeric(Y[g, ]), des$W, f$alpha[g, ], f$psi[g])
  Z <- des$Zs[[s]]
  Lt <- stats::lm.wfit(Z, L[i, , drop = FALSE], wr$w[i])$residuals
  step <- solve(crossprod(Lt * sqrt(wr$w[i])), colSums(Lt * wr$r[i]))
  b <- f$alpha[g, des$niche_cols] + step
  expect_equal(unname(x$beta[g, s, ]), unname(b[seq_along(nc$tested)]), tolerance = 1e-6)
  r <- testSpiDE(fit)
  expect_true(any(is.finite(r@table$t)))
})
