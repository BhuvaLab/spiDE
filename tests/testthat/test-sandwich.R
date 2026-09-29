# The sandwich engine's CR2 sandwich and Bell-McCaffrey df, against
# clubSandwich on the working (IRLS) linear model of the same converged NB fit:
# weighted least squares of the working response z = eta + (y - mu) / mu on the
# full design (patient indicators included), weights w = mu / (1 + psi mu),
# clusters = patients. The engine's low-rank CR2 on the within-patient
# (absorbed) design must equal clubSandwich's CR2 on that full design -- the
# Pustejovsky & Tipton (2018) absorption result (clubSandwich itself gives the
# same answer on the absorbed and the full design).
# (clubSandwich on a glm() object is not the oracle: it linearises the GLM
# differently, and differs from both by up to ~20% SE on a high-count gene.)
#
# The IRLS weights are inverse variances, so clubSandwich must be told so
# (inverse_var = TRUE): left to infer it for a weighted lm, it takes an identity
# working variance, which moves the Bell-McCaffrey df by up to ~10% on a
# low-count gene (13.1 against 11.8 on G3 here) while barely moving the SE.
# With it, the df agree to two decimals and the SE within 0.3% (2026-09-30).

test_that("CR2 and Bell-McCaffrey df agree with clubSandwich on the working model", {
  skip_if_not_installed("clubSandwich")
  spe <- buildNiches(.toySPE(n_samples = 16, n_per = 150, n_genes = 6, seed = 5), sigma = 30,
                     verbose = FALSE)
  ct <- as.character(spe$cell_type); smp <- as.character(spe$sample_id)
  ik <- spiDE:::.indexCells(ct, smp, "A", 5L)
  NM <- SingleCellExperiment::reducedDim(spe, "Niche30")
  nc <- spiDE:::.nicheColumns(NM, "A", ik)
  L <- log1p(NM[ik, nc$cols, drop = FALSE])
  Y <- SummarizedExperiment::assay(spe, "counts")[c("G1", "G2", "G3"), ik, drop = FALSE]
  pat <- factor(smp[ik])
  trt <- as.numeric(spe$condition[ik] == "Responder")
  des <- spiDE:::.indexDesign(L, NULL, pat, trt = trt, tested = nc$tested)
  fit <- spiDE:::.fitIndexGLM(Y, des)
  got <- spiDE:::.sandwichCR2(fit, des, Y)
  Xd <- des$W[, -seq_len(des$npat), drop = FALSE]
  colnames(Xd) <- make.names(colnames(Xd))
  for (g in rownames(Y)) {
    eta <- as.numeric(des$W %*% fit$alpha[g, ])
    mu <- exp(eta)
    w <- mu / (1 + fit$psi[g] * mu)
    z <- eta + (as.numeric(Y[g, ]) - mu) / mu
    m <- stats::lm(z ~ 0 + pat + Xd, weights = w)
    V <- clubSandwich::vcovCR(m, cluster = pat, type = "CR2", inverse_var = TRUE)
    cs <- clubSandwich::coef_test(m, vcov = V, test = "Satterthwaite")
    o <- cs[match(paste0("Xd", make.names(colnames(des$W)[des$tested])), cs$Coef), ]
    mine <- got[got$gene == g, ]
    expect_equal(mine$estimate, unname(o$beta), tolerance = 1e-3)
    expect_equal(mine$se, unname(o$SE), tolerance = 0.005)
    expect_equal(mine$df, unname(o$df_Satt), tolerance = 0.005)
  }
})
