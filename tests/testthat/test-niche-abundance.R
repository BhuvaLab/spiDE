# testNicheAbundance(): the between-patient niche-abundance association (up to
# spiDE 0.99.22, compositionTest()). This is what the per-patient intercepts of
# fitSpiDE() ABSORB; here it is tested on its own terms, with S units and a
# limma moderated t.

# 200 cells per sample and a bandwidth of 50: at lower density the niche
# covariate's per-sample MEAN is dominated by cell-placement noise (between-
# sample sd 0.07 against a within-sample 0.25 at 60 cells), and no estimator
# can see a between-sample effect the covariate does not carry. Here the ratio
# is 0.7 and the planted confound reaches the patient-level test at t ~ 4.6.
spe_c <- buildNiches(.toySPE(n_samples = 12, n_per = 200, n_genes = 12,
                             composition = 3, seed = 3), sigma = 50)

test_that("testNicheAbundance returns the documented tidy schema", {
  ct <- testNicheAbundance(spe_c, condition = "condition", sigma = 50, verbose = FALSE)
  expect_s3_class(ct, "data.frame")
  expect_true(all(c("gene", "index", "niche", "term", "estimate", "t", "p",
                    "n_patients", "leverage", "downweighted", "q", "q.global") %in% names(ct)))
  expect_setequal(unique(ct$term), c("niche", "condition:niche"))
  # an index type is never tested against its own niche
  expect_false(any(ct$index == ct$niche))
  expect_true(all(ct$p >= 0 & ct$p <= 1))
  expect_true(all(ct$q >= ct$p - 1e-12))
})

test_that("the planted between-sample confound is a patient-level finding", {
  # .toySPE(composition = 3) shifts G2's baseline in Responders' A cells in
  # proportion to the sample's B-cell prevalence, constant within (sample, A).
  # That is a between-sample association of expression in A with B density
  # around A -- exactly what this test exists to report, on the
  # condition:niche term since the shift is condition-specific.
  ct <- testNicheAbundance(spe_c, condition = "condition", sigma = 50, verbose = FALSE)
  g2 <- ct[ct$gene == "G2" & ct$index == "A" & ct$niche == "B" &
             ct$term == "condition:niche", ]
  expect_equal(nrow(g2), 1L)
  expect_gt(g2$t, 2)
  # and it is the strongest interaction in the A x B pair
  ab <- ct[ct$index == "A" & ct$niche == "B" & ct$term == "condition:niche", ]
  expect_equal(ab$gene[which.max(ab$t)], "G2")
})

test_that("without a condition only the pooled niche term is reported", {
  ct <- testNicheAbundance(spe_c, condition = NULL, sigma = 50, verbose = FALSE)
  expect_equal(unique(ct$term), "niche")
})

test_that("min.cells drops thin samples and a sample-level covariate enters the design", {
  ct_all <- testNicheAbundance(spe_c, condition = "condition", sigma = 50,
                            min.cells = 1L, verbose = FALSE)
  # ~55-85 A cells per sample, so min.cells = 70 drops some samples while
  # leaving enough to fit
  ct_strict <- testNicheAbundance(spe_c, condition = "condition", sigma = 50,
                               min.cells = 70L, verbose = FALSE)
  expect_true(min(ct_strict$n_patients) <= min(ct_all$n_patients))
  expect_true(all(ct_strict$n_patients >= 3L))
  # Age is constant within sample; it must be accepted here (fitSpiDE() rejects
  # it because the per-patient intercepts absorb it)
  ct_cov <- testNicheAbundance(spe_c, condition = "condition", sigma = 50,
                            covariates = "Age", verbose = FALSE)
  expect_true(nrow(ct_cov) > 0)
})

test_that("index and niche restrictions are honoured", {
  ct <- testNicheAbundance(spe_c, condition = "condition", sigma = 50,
                        index = "A", niche = c("B", "C"), verbose = FALSE)
  expect_equal(unique(ct$index), "A")
  expect_setequal(unique(ct$niche), c("B", "C"))
})

test_that("with a condition, niche is the mean of the two conditions' associations, condition:niche their difference", {
  # no leverage cap, no covariates: the coefficients are those of per-condition simple regressions
  set.seed(4)
  Y <- matrix(rnorm(5 * 12), 5, 12, dimnames = list(paste0("g", 1:5), NULL))
  df <- data.frame(niche = rnorm(12), condition = factor(rep(c("a", "b"), each = 6), levels = c("a", "b")))
  r <- spiDE:::.abundancePair(Y, df, max.leverage = Inf)
  slope <- function(i) apply(Y[, i], 1, function(y) stats::coef(stats::lm(y ~ df$niche[i]))[[2]])
  sa <- slope(1:6); sb <- slope(7:12)
  expect_equal(r$terms$niche$estimate, unname((sa + sb) / 2), tolerance = 1e-10)
  expect_equal(r$terms$`condition:niche`$estimate, unname(sb - sa), tolerance = 1e-10)
  # one condition with a single sample: only the pooled association, across all samples
  df1 <- transform(df, condition = factor(c(rep("a", 11), "b"), levels = c("a", "b")))
  r1 <- spiDE:::.abundancePair(Y, df1, max.leverage = Inf)
  expect_equal(names(r1$terms), "niche")
  expect_equal(r1$terms$niche$estimate,
               unname(apply(Y, 1, function(y) stats::coef(stats::lm(y ~ df$niche))[[2]])), tolerance = 1e-10)
})

test_that("q.global is one Benjamini-Hochberg family per term", {
  ct <- testNicheAbundance(spe_c, condition = "condition", sigma = 50, verbose = FALSE)
  for (tm in unique(ct$term)) {
    i <- ct$term == tm
    expect_equal(ct$q.global[i], stats::p.adjust(ct$p[i], "BH"))
  }
})

test_that("the leverage cap down-weights an influential sample, and only from the design", {
  set.seed(5)
  X <- cbind(1, c(rnorm(19), 12))                       # one sample far out on the niche axis
  cap <- spiDE:::.capLeverage(X, max.leverage = 3)
  expect_gt(cap$leverage, 3)
  expect_equal(cap$downweighted, 1L)
  expect_lt(cap$w[20], 1)
  # every sample ends at or below the cap, or at the weight floor (1/20)
  hw <- stats::hat(X * sqrt(cap$w), intercept = FALSE)
  expect_true(all(hw <= 3 * ncol(X) / nrow(X) * 1.001 | cap$w <= 0.05 + 1e-12))
  # a milder outlier is brought exactly to the cap
  X2 <- cbind(1, c(rnorm(19), 6))
  c2 <- spiDE:::.capLeverage(X2, max.leverage = 3)
  expect_gt(min(c2$w), 0.05)
  expect_lte(max(stats::hat(X2 * sqrt(c2$w), intercept = FALSE)), 3 * ncol(X2) / nrow(X2) * 1.001)
  expect_identical(spiDE:::.capLeverage(X, Inf)$w, rep(1, 20))
  # the weighted fit equals weighted least squares with those weights
  Y <- matrix(rnorm(3 * 20), 3, 20, dimnames = list(paste0("g", 1:3), NULL))
  r <- spiDE:::.abundancePair(Y, data.frame(niche = X[, 2]), max.leverage = 3)
  expect_equal(r$downweighted, 1L)
  expect_equal(r$terms$niche$estimate,
               unname(apply(Y, 1, function(y) stats::coef(stats::lm(y ~ X[, 2], weights = cap$w))[[2]])), tolerance = 1e-10)
  # no influential sample: nothing is down-weighted
  expect_equal(spiDE:::.capLeverage(cbind(1, seq(-1, 1, length.out = 20)), 3)$downweighted, 0L)
  expect_error(testNicheAbundance(spe_c, sigma = 50, max.leverage = 1, verbose = FALSE), "max.leverage")
})

test_that("a pair the cap would leave without support drops out instead of stopping the run", {
  # 14 samples, one group of two with an extreme sample: the cap pushes its weight to the floor
  set.seed(6)
  Y <- matrix(rnorm(4 * 14), 4, 14, dimnames = list(paste0("g", 1:4), NULL))
  df <- data.frame(niche = c(rnorm(12), 0, 25), condition = factor(c(rep("a", 12), "b", "b"), levels = c("a", "b")))
  cap <- spiDE:::.capLeverage(model.matrix(~ niche * condition, transform(df, condition = as.numeric(condition == "b") - 0.5)), 2)
  expect_gte(min(cap$w), 0.05)
  r <- expect_no_error(spiDE:::.abundancePair(Y, df, max.leverage = 2))
  if (!is.null(r)) expect_true(all(vapply(r$terms, function(x) all(is.finite(x$p)), TRUE)))
})

test_that("gene sets: a set's score is its genes' mean standardised expression, less the others' when competitive", {
  set.seed(7)
  Y <- matrix(rnorm(30 * 10, 5), 30, 10, dimnames = list(paste0("g", 1:30), paste0("s", 1:10)))
  sets <- list(a = paste0("g", 1:6), b = paste0("g", 4:12), small = c("g1", "g2"), absent = c("x", "y", "z", "w", "v"))
  Z <- t(scale(t(Y)))
  sc <- spiDE:::.abundanceSetScores(Y, sets, "self-contained", 5L, 500L)
  expect_equal(rownames(sc$S), c("a", "b"))                       # too small and absent sets dropped
  expect_equal(unname(sc$size), c(6L, 9L))
  expect_equal(unname(sc$S["a", ]), unname(colMeans(Z[1:6, ])), tolerance = 1e-12)
  cp <- spiDE:::.abundanceSetScores(Y, sets, "competitive", 5L, 500L)
  expect_equal(unname(cp$S["b", ]), unname(colMeans(Z[4:12, ]) - colMeans(Z[-(4:12), ])), tolerance = 1e-12)
})

test_that("testNicheAbundance(genesets = ) tests sets like genes, with the same terms and families", {
  sets <- list(with_G2 = c("G2", "G3", "G4", "G5", "G6"), other = c("G7", "G8", "G9", "G10", "G11"))
  ab <- testNicheAbundance(spe_c, condition = "condition", sigma = 50, genesets = sets, verbose = FALSE)
  expect_true(all(c("set", "size", "index", "niche", "term", "estimate", "t", "p", "n_patients",
                    "leverage", "downweighted", "q", "q.global", "direction") %in% names(ab)))
  expect_false("gene" %in% names(ab))
  expect_setequal(unique(ab$term), c("niche", "condition:niche"))
  for (tm in unique(ab$term)) expect_equal(ab$q.global[ab$term == tm], stats::p.adjust(ab$p[ab$term == tm], "BH"))
  # the planted G2 confound (A cells, B prevalence, Responders) reaches the set holding G2
  x <- ab[ab$index == "A" & ab$niche == "B" & ab$term == "condition:niche", ]
  expect_gt(x$t[x$set == "with_G2"], x$t[x$set == "other"])
  # min.size above every set: nothing to test
  expect_error(testNicheAbundance(spe_c, sigma = 50, genesets = sets, min.size = 50L, verbose = FALSE), "min.size")
  expect_error(testNicheAbundance(spe_c, sigma = 50, genesets = c("G1", "G2"), verbose = FALSE), "genesets")
})
