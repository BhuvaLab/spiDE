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
  expect_lte(max(stats::hat(X * sqrt(cap$w), intercept = FALSE)), 3 * ncol(X) / nrow(X) * 1.001)
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
