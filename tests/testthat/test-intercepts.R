# The fit keeps what the plots read of each patient (0.99.34): its intercept per
# gene in the shared, condition-free fit, its mean log1p niche densities and its
# mean log library size.

spe <- buildNiches(.toySPE(n_samples = 16, n_per = 150, seed = 7), sigma = 30)

test_that("the slopes engine keeps each patient's intercept, niche means and depth", {
  fit <- fitSpiDE(spe, sigma = 30, index = "A", verbose = FALSE)
  xi <- fit@index$A
  expect_equal(dim(xi$intercept), c(length(xi$genes), length(xi$patients)))
  expect_equal(dimnames(xi$intercept), list(xi$genes, xi$patients))
  expect_equal(dim(xi$niche_mean), c(length(xi$patients), length(xi$niches)))
  expect_equal(rownames(xi$niche_mean), xi$patients)
  expect_named(xi$loglib_mean, xi$patients)
  expect_true(all(is.finite(xi$intercept)))
  expect_true(all(is.finite(xi$loglib_mean)))
})

test_that("the stored intercepts are the shared fit's patient intercepts", {
  fit <- fitSpiDE(spe, sigma = 30, index = "A", verbose = FALSE)
  xi <- fit@index$A
  cd <- SummarizedExperiment::colData(spe)
  Y <- SummarizedExperiment::assay(spe, "counts")
  cov <- .cellCovariates(cd, character(), "loglib", Y)
  usable <- is.finite(cov[, "loglib"])
  ct <- ifelse(usable, as.character(cd$cell_type), NA_character_)
  smp <- as.character(cd$sample_id)
  ik <- .indexCells(ct, smp, "A", 10L)
  NM <- as.matrix(SingleCellExperiment::reducedDim(spe, "Niche30"))
  nc <- .nicheColumns(NM, "A", ik)
  L <- log1p(NM[ik, c(nc$tested, setdiff(nc$cols, nc$tested)), drop = FALSE])
  des <- .indexDesign(L, .indexCovariates(cov[ik, , drop = FALSE], "loglib"), factor(smp[ik]),
                      tested = nc$tested)
  f0 <- .fitIndexGLM(Y[xi$genes, ik, drop = FALSE], des)
  expect_equal(unname(xi$intercept), unname(f0$alpha[, des$intercept_cols]), tolerance = 1e-10)
  expect_equal(unname(xi$niche_mean[, "B"]),
               as.numeric(tapply(L[, "B"], factor(smp[ik]), mean)), tolerance = 1e-12)
})

test_that("the sandwich engine keeps them too, with the patients' cell counts", {
  fit <- fitSpiDE(spe, condition = "condition", sigma = 30, index = "A", engine = "sandwich",
                  verbose = FALSE)
  xi <- fit@index$A
  expect_equal(dim(xi$intercept), c(length(xi$genes), length(xi$patients)))
  expect_named(xi$ncells, xi$patients)
  expect_equal(sum(xi$ncells), sum(spe$cell_type == "A"))
})

test_that("patientIntercepts() returns them in long form with the patient columns", {
  fit <- fitSpiDE(spe, condition = "condition", sigma = 30, index = "A", verbose = FALSE)
  d <- patientIntercepts(fit, gene = c("G1", "G2"))
  expect_equal(nrow(d), 2L * length(fit@index$A$patients))
  expect_true(all(c("gene", "index", "patient", "intercept", "ncells", "condition") %in% colnames(d)))
  expect_equal(d$intercept[d$gene == "G1" & d$patient == "S1"], fit@index$A$intercept["G1", "S1"])
  res <- testSpiDE(fit)
  expect_equal(patientIntercepts(res, gene = "G1"), patientIntercepts(fit, gene = "G1"))
})

test_that("a fit without intercepts (spiDE <= 0.99.33) is refused with a pointer to refit", {
  fit <- fitSpiDE(spe, sigma = 30, index = "A", verbose = FALSE)
  fit@index$A$intercept <- NULL
  expect_error(patientIntercepts(fit), "refit")
})
