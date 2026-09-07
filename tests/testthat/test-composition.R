# compositionTest(): the between-sample composition association, at the
# patient level. This is what the nested intercept ABSORBS in fitSpiDE(); here
# it is tested on its own terms, with S units and a limma moderated t.

# 200 cells per sample and a bandwidth of 50: at lower density the niche
# covariate's per-sample MEAN is dominated by cell-placement noise (between-
# sample sd 0.07 against a within-sample 0.25 at 60 cells), and no estimator
# can see a between-sample effect the covariate does not carry. Here the ratio
# is 0.7 and the planted confound reaches the patient-level test at t ~ 4.6.
spe_c <- buildNiches(.toySPE(n_samples = 12, n_per = 200, n_genes = 12,
                             composition = 3, seed = 3), sigma = 50)

test_that("compositionTest returns the documented tidy schema", {
  ct <- compositionTest(spe_c, condition = "condition", sigma = 50, verbose = FALSE)
  expect_s3_class(ct, "data.frame")
  expect_true(all(c("gene", "ct_index", "ct_niche", "term", "coef", "t", "p",
                    "fdr", "fdr.global", "n_samples") %in% names(ct)))
  expect_setequal(unique(ct$term), c("niche", "condition:niche"))
  # an index type is never tested against its own niche
  expect_false(any(ct$ct_index == ct$ct_niche))
  expect_true(all(ct$p >= 0 & ct$p <= 1))
  expect_true(all(ct$fdr >= ct$p - 1e-12))
})

test_that("the planted between-sample confound is a patient-level finding", {
  # .toySPE(composition = 3) shifts G2's baseline in Responders' A cells in
  # proportion to the sample's B-cell prevalence, constant within (sample, A).
  # That is a between-sample association of expression in A with B density
  # around A -- exactly what this test exists to report, on the
  # condition:niche term since the shift is condition-specific.
  ct <- compositionTest(spe_c, condition = "condition", sigma = 50, verbose = FALSE)
  g2 <- ct[ct$gene == "G2" & ct$ct_index == "A" & ct$ct_niche == "B" &
             ct$term == "condition:niche", ]
  expect_equal(nrow(g2), 1L)
  expect_gt(g2$t, 2)
  # and it is the strongest interaction in the A x B pair
  ab <- ct[ct$ct_index == "A" & ct$ct_niche == "B" & ct$term == "condition:niche", ]
  expect_equal(ab$gene[which.max(ab$t)], "G2")
})

test_that("without a condition only the pooled niche term is reported", {
  ct <- compositionTest(spe_c, condition = NULL, sigma = 50, verbose = FALSE)
  expect_equal(unique(ct$term), "niche")
})

test_that("min.cells drops thin samples and a sample-level covariate enters the design", {
  ct_all <- compositionTest(spe_c, condition = "condition", sigma = 50,
                            min.cells = 1L, verbose = FALSE)
  # ~55-85 A cells per sample, so min.cells = 70 drops some samples while
  # leaving enough to fit
  ct_strict <- compositionTest(spe_c, condition = "condition", sigma = 50,
                               min.cells = 70L, verbose = FALSE)
  expect_true(min(ct_strict$n_samples) <= min(ct_all$n_samples))
  expect_true(all(ct_strict$n_samples >= 3L))
  # Age is constant within sample; it must be accepted here (fitSpiDE rejects
  # it under random != "none" because the per-sample intercept absorbs it)
  ct_cov <- compositionTest(spe_c, condition = "condition", sigma = 50,
                            covariates = "Age", verbose = FALSE)
  expect_true(nrow(ct_cov) > 0)
})

test_that("index and niche restrictions are honoured", {
  ct <- compositionTest(spe_c, condition = "condition", sigma = 50,
                        index = "A", niche = c("B", "C"), verbose = FALSE)
  expect_equal(unique(ct$ct_index), "A")
  expect_setequal(unique(ct$ct_niche), c("B", "C"))
})
