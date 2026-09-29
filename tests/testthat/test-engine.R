# The per-patient engines' building blocks, checked against direct formulas.

spe16 <- buildNiches(.toySPE(n_samples = 16, n_per = 120, n_genes = 12, seed = 2), sigma = 30,
                     verbose = FALSE)

test_that("the index design absorbs a patient block and tags the tested columns", {
  L <- matrix(rnorm(40), 20, 2, dimnames = list(NULL, c("B", "C")))
  pat <- factor(rep(c("p1", "p2"), each = 10))
  d0 <- spiDE:::.indexDesign(L, NULL, pat, tested = "B")
  expect_equal(d0$npat, 2L)
  expect_equal(colnames(d0$W)[d0$tested], "B")
  expect_true(all(d0$absorb[seq_len(2)]) && !any(d0$absorb[-seq_len(2)]))
  d1 <- spiDE:::.indexDesign(L, NULL, pat, trt = rep(0:1, each = 10), tested = c("B", "C"))
  expect_equal(colnames(d1$W)[d1$tested], c("condition:B", "condition:C"))
})

test_that("an index type is never tested against its own (merged) niche", {
  spe <- mergeNiches(buildNiches(.toySPE(), sigma = 20, verbose = FALSE),
                     groups = list(AC = c("A", "C")), sigma = 20)
  NM <- SingleCellExperiment::reducedDim(spe, "Niche20")
  ct <- as.character(spe$cell_type)
  map <- S4Vectors::metadata(spe)$spiDE_niche_groups[["Niche20"]]
  nA <- spiDE:::.nicheColumns(NM, "A", which(ct == "A"), group_map = map)
  nB <- spiDE:::.nicheColumns(NM, "B", which(ct == "B"), group_map = map)
  expect_false("AC" %in% nA$cols)
  expect_true("AC" %in% nB$tested)
})

fit16 <- fitSpiDE(spe16, index = "A", sigma = 30, verbose = FALSE)

test_that("fitSpiDE (slopes) returns per-patient slopes with both variances", {
  x <- fit16@index$A
  expect_s4_class(fit16, "SpiDEFit")
  expect_equal(dim(x$beta), c(length(x$genes), length(x$patients), length(x$niches)))
  expect_true(all(x$v_tile[is.finite(x$v_tile)] > 0))
  expect_true(all(x$factor >= 1))
})

test_that("a patient slope is one Fisher-scoring step from the pooled fit", {
  x <- fit16@index$A
  spe <- spe16
  ct <- as.character(spe$cell_type); smp <- as.character(spe$sample_id)
  ik <- spiDE:::.indexCells(ct, smp, "A", 10L)
  NM <- SingleCellExperiment::reducedDim(spe, "Niche30")
  nc <- spiDE:::.nicheColumns(NM, "A", ik)
  L <- log1p(NM[ik, c(nc$tested, setdiff(nc$cols, nc$tested)), drop = FALSE])
  Y <- SummarizedExperiment::assay(spe, "counts")[x$genes, ik, drop = FALSE]
  cov <- scale(matrix(log(Matrix::colSums(SummarizedExperiment::assay(spe, "counts")))[ik],
                      dimnames = list(NULL, "loglib")), scale = FALSE)
  des <- spiDE:::.indexDesign(L, cov, factor(smp[ik]), tested = nc$tested)
  fit <- spiDE:::.fitIndexGLM(Y, des)
  g <- 1L; s <- 3L
  i <- which(des$patient == s)
  wr <- spiDE:::.workingWR(as.numeric(Y[g, ]), des$W, fit$alpha[g, ], fit$psi[g])
  Lt <- sweep(L[i, , drop = FALSE], 2, colSums(L[i, , drop = FALSE] * wr$w[i]) / sum(wr$w[i]))
  step <- solve(crossprod(Lt * sqrt(wr$w[i])), colSums(Lt * wr$r[i]))
  b <- fit$alpha[g, des$niche_cols] + step
  expect_equal(unname(x$beta[g, s, ]), unname(b[seq_along(nc$tested)]), tolerance = 1e-6)
})

test_that("the pooled and condition tests equal a direct limma fit with the effective df", {
  x <- fit16@index$A
  j <- 1L
  b <- matrix(x$beta[, , j], nrow = length(x$genes))
  vt <- matrix(x$v_tile[, , j], nrow = length(x$genes))
  tau <- spiDE:::.dlTau2(x$beta, x$v_tile)[, j]
  trt <- as.numeric(fit16@patients$condition[match(x$patients, fit16@patients$patient)] == "Responder")
  got <- spiDE:::.slopeColumnTest(b, vt, tau, x$mean_expr, trt = trt)
  w <- 1 / (vt + tau); w[!is.finite(w) | !is.finite(b)] <- NA
  w <- w / rowMeans(w, na.rm = TRUE)
  f <- limma::lmFit(b, cbind(1, trt), weights = w)
  f$Amean <- log(x$mean_expr + 1e-3)
  f <- limma::eBayes(f, robust = TRUE, trend = TRUE)
  neff <- rowSums(w, na.rm = TRUE)^2 / rowSums(w^2, na.rm = TRUE)
  df <- pmin(f$df.total, pmax(neff - 2, 1))
  expect_equal(got$estimate, unname(f$coefficients[, 2]), tolerance = 1e-10)
  expect_equal(got$df, unname(df), tolerance = 1e-10)
  ok <- is.finite(got$t)
  expect_equal(got$p[ok], unname(2 * pt(-abs(f$t[ok, 2]), df[ok])), tolerance = 1e-10)
})

test_that("the pooled test and its filter never look at the condition", {
  r1 <- testSpiDE(fit16, condition = "condition")
  pt <- fit16@patients
  set.seed(4)
  lab <- stats::setNames(sample(as.character(pt$condition)), pt$patient)
  fit_perm <- fit16
  fit_perm@patients$perm <- factor(lab[pt$patient])
  r2 <- testSpiDE(fit_perm, condition = "perm")
  p1 <- r1@table[r1@table$test == "pooled", ]
  p2 <- r2@table[r2@table$test == "pooled", ]
  expect_identical(p1[, c("gene", "niche", "p", "q")], p2[, c("gene", "niche", "p", "q")])
})

test_that("results are identical serially and in parallel", {
  skip_on_os("windows")
  f2 <- fitSpiDE(spe16, index = "A", sigma = 30, verbose = FALSE,
                 BPPARAM = BiocParallel::MulticoreParam(2))
  expect_equal(f2@index$A$beta, fit16@index$A$beta, tolerance = 1e-8)
  expect_equal(f2@index$A$v_tile, fit16@index$A$v_tile, tolerance = 1e-8)
})

test_that("testSpiDE runs the pooled test alone without a condition", {
  r0 <- testSpiDE(fit16)
  expect_identical(unique(r0@table$test), "pooled")
  expect_identical(unique(results(r0)$test), "pooled")
  r1 <- testSpiDE(fit16, condition = "condition", procedure = "filtered")
  cond <- r1@table[r1@table$test == "condition", ]
  pool <- r1@table[r1@table$test == "pooled", ]
  pass <- paste(pool$gene, pool$niche)[pool$q < 0.05]
  expect_true(all(paste(cond$gene, cond$niche)[cond$in_family] %in% pass))
})

test_that("a mixed-model object from spiDE <= 0.99.22 is refused with the archive pointer", {
  legacy <- fit16
  attr(legacy, "covtype") <- c("Niche", "ResponseNiche")
  expect_error(testSpiDE(legacy), "spiDEmixed::readSpiDE")
  expect_output(show(legacy), "legacy spiDE mixed-model fit")
})

test_that("a gene the solver could not fit drops out; a bound dispersion is refitted at the bound", {
  set.seed(9)
  pat <- factor(rep(sprintf("p%d", 1:6), each = 60))
  L <- matrix(rnorm(360), 360, 1, dimnames = list(NULL, "B"))
  des <- spiDE:::.indexDesign(L, NULL, pat, tested = "B")
  mu <- exp(1 + 0.3 * L[, 1])
  Y <- rbind(pois = rpois(360, mu), nb = rnbinom(360, size = 2, mu = mu))
  rownames(Y) <- c("pois", "nb")
  fit <- spiDE:::.fitIndexGLM(Y, des)
  expect_equal(unname(fit$status["genes"]), 2)
  expect_equal(unname(fit$status["psi_at_bound"]), 1)   # the Poisson gene
  expect_equal(unname(fit$psi["pois"]), 1e-3)
  expect_equal(unname(fit$alpha["pois", des$niche_cols]), 0.3, tolerance = 0.1)
  expect_gt(unname(fit$psi["nb"]), 0.1)
  # a gene marked unpolished keeps no estimate
  fake <- fit
  fake$polish$polished <- c(TRUE, FALSE)
  fake$polish$psi_bound <- c(FALSE, FALSE)
  out <- spiDE:::.fitStatus(fake, Y, des, NULL)
  expect_true(all(is.na(out$alpha["nb", ])))
  expect_true(is.na(out$psi["nb"]))
  expect_equal(unname(out$status["not_fitted"]), 1)
})
