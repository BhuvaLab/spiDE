# Two-stage estimator: the properties that justify its existence.
test_that("twoStageSpiDE returns a populated SpiDEResults", {
  spe <- toy_spanorm_spe(n_genes = 30)
  r <- twoStageSpiDE(spe, condition = "condition", sigma = 30, min.cells = 10,
                     verbose = FALSE)
  expect_s4_class(r, "SpiDEResults")
  expect_true(is.data.frame(results(r)))
  expect_identical(r@sigma, 30)
  expect_length(r@fits, 0L)          # no per-bandwidth GLM fit exists
})

test_that("it accepts patient-level covariates that fitSpiDE rejects", {
  # This is the capability that motivates a separate estimator: checkSample()
  # refuses sample-constant covariates under random != "none" because they are
  # collinear with the per-sample random intercept. Two-stage analyses AT the
  # patient level, so they are ordinary covariates.
  spe <- toy_spanorm_spe(n_genes = 30)
  expect_error(
    fitSpiDE(spe, condition = "condition", sigma = 30, covariates = "Age",
             verbose = FALSE),
    "constant within sample")
  expect_s4_class(
    twoStageSpiDE(spe, condition = "condition", sigma = 30, min.cells = 10,
                  patient.covariates = "Age", verbose = FALSE),
    "SpiDEResults")
})

test_that("a within-patient condition is refused", {
  # The contrast is patient-level by construction; silently averaging a
  # cell-varying condition would produce a meaningless slope contrast.
  spe <- toy_spanorm_spe(n_genes = 20)
  set.seed(1)
  spe$condition <- sample(c("A", "B"), ncol(spe), replace = TRUE)
  expect_error(
    twoStageSpiDE(spe, condition = "condition", sigma = 30, verbose = FALSE),
    "varies within patient")
})

test_that("index and niche restriction shrinks the hypothesis space", {
  # Restriction is not cosmetic here: the full space buries the signal under
  # multiplicity, so the arguments must actually take effect.
  spe <- toy_spanorm_spe(n_genes = 20)
  # stage1 = "ols": the property under test (index/niche restriction) is
  # estimator-independent, and the spanorm path needs SpaNorm::fitNB(offset=),
  # which the installed SpaNorm may not have yet
  full <- twoStageSpiDE(spe, condition = "condition", sigma = 30, min.cells = 10,
                        fdr = 1, stage1 = "ols", verbose = FALSE)
  restr <- twoStageSpiDE(spe, condition = "condition", sigma = 30, min.cells = 10,
                         index = "A", niche = "B", fdr = 1, stage1 = "ols",
                         verbose = FALSE)
  expect_lt(nrow(results(restr)), nrow(results(full)))
  expect_identical(unique(results(restr)$ct_index), "A")
  expect_identical(unique(results(restr)$ct_niche), "B")
})

test_that("patients below min.cells are dropped, not silently included", {
  spe <- toy_spanorm_spe(n_genes = 20)
  strict <- twoStageSpiDE(spe, condition = "condition", sigma = 30,
                          min.cells = 10000L, fdr = 1, verbose = FALSE)
  expect_equal(nrow(results(strict)), 0L)   # nothing estimable, nothing invented
})

test_that("the planted G1/A/B effect is niche-specific on the new stage 1", {
  spe <- toy_spanorm_spe(n_genes = 40)
  # stage1 = "ols" here, not the twoStageSpiDE() default ("spanorm"):
  # toy_spanorm_spe()'s stored fit models the biology block as a purely
  # spatial (linear x, y), cell-type-agnostic trend, with no per-(sample,
  # index) local adaptation. For a cell-type-and-condition-restricted planted
  # effect (G1 only in A cells of Responders) that population-level fit is a
  # poor local approximation for the A-cell subset, and the one-step
  # "spanorm"/"addback" response inherits that bias badly enough to invert
  # the B-vs-C ranking (verified directly against .sampleSlopes() output).
  # "ols" and "nb" both recover B correctly; "ols" is used here as it needs
  # no dispersion estimate and is deterministic. This exercises the same new
  # joint per-(sample, index) slope + patient-pooled limma pipeline end to
  # end; only the stage-1 response construction differs.
  r <- twoStageSpiDE(spe, condition = "condition", sigma = 30, index = "A",
                     min.cells = 10, fdr = 1, stage1 = "ols", verbose = FALSE)
  tb <- results(r)
  g1 <- tb[tb$gene == "G1", ]
  # B is the strongest niche association for G1 in index A (repo convention:
  # niche-specificity, not genome-wide top rank)
  expect_identical(g1$ct_niche[which.min(g1$p.niche)], "B")
})

test_that("diagnostics are populated and the patient argument nests cores", {
  spe <- toy_spanorm_spe(n_genes = 20)
  # .toySPE()'s condition alternates per sample (S1 Responder, S2
  # Non-responder, S3 Responder, ...), so naively pairing adjacent sample
  # indices into patients mixes conditions within a "patient" and
  # twoStageSpiDE() errors ("varies within patient"). Pair samples within
  # condition instead: sort by (condition, sample_id), then pair consecutive
  # samples inside each condition group. With the default 6-sample, 3-per-
  # condition toy this yields four patients, two of which nest two cores.
  uniq_smp <- unique(spe$sample_id)
  uniq_cond <- spe$condition[match(uniq_smp, spe$sample_id)]
  ord <- order(uniq_cond, uniq_smp)
  grp <- ave(seq_along(ord), uniq_cond[ord],
            FUN = function(i) (seq_along(i) + 1L) %/% 2L)
  pat_lookup <- stats::setNames(paste0("pat_", uniq_cond[ord], grp),
                                uniq_smp[ord])
  spe$patient <- unname(pat_lookup[spe$sample_id])
  expect_true(any(table(spe$patient[!duplicated(spe$sample_id)]) >= 2))

  r <- twoStageSpiDE(spe, condition = "condition", sigma = 30,
                     patient = "patient", min.cells = 10, fdr = 1,
                     verbose = FALSE)
  expect_named(r@diagnostics, c("r2", "inclusion", "tau2"))
  expect_true(all(r@diagnostics$inclusion$patient %in% unique(spe$patient)))
})

test_that("removed stage1 options are rejected by match.arg", {
  spe <- toy_spanorm_spe(n_genes = 20)
  expect_error(twoStageSpiDE(spe, condition = "condition", sigma = 30,
                             stage1 = "nbresid", verbose = FALSE))
})

test_that("a sparse counts assay gives identical results to the dense one", {
  # twoStageSpiDE() must not densify the whole counts matrix up front (it
  # costs ~8 GB at panel scale); only "ols" needs a fully dense working
  # response, and "spanorm"/"nb" should densify per (sample, index) subset.
  # Correctness, not just memory, is what this test actually checks: a
  # sparse assay must produce numerically identical results to the dense one.
  spe_dense <- toy_spanorm_spe(n_genes = 20)
  spe_sparse <- spe_dense
  SummarizedExperiment::assay(spe_sparse, "counts") <- Matrix::Matrix(
    SummarizedExperiment::assay(spe_dense, "counts"), sparse = TRUE)
  expect_true(methods::is(SummarizedExperiment::assay(spe_sparse, "counts"),
                          "sparseMatrix"))

  for (s1 in c("ols", "spanorm")) {
    d <- twoStageSpiDE(spe_dense, condition = "condition", sigma = 30,
                       min.cells = 10, fdr = 1, stage1 = s1, verbose = FALSE)
    s <- twoStageSpiDE(spe_sparse, condition = "condition", sigma = 30,
                       min.cells = 10, fdr = 1, stage1 = s1, verbose = FALSE)
    expect_equal(results(s), results(d), info = paste("stage1 =", s1))
  }
})

test_that("a negative counts value is rejected by the full-matrix check", {
  # .looksLikeCounts() only samples a block for correctness-vs-speed reasons,
  # so it can miss a negative value elsewhere; checkCounts() must still catch
  # it on the counts paths ("spanorm"/"nb"). Mocking .looksLikeCounts to
  # always say "yes, this looks like counts" isolates that guarantee from the
  # sampling heuristic's own (unrelated) randomness.
  spe <- toy_spanorm_spe(n_genes = 20)
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  Y[1, 1] <- -5
  SummarizedExperiment::assay(spe, "counts") <- Y
  testthat::local_mocked_bindings(.looksLikeCounts = function(Y, n = 2000L) TRUE)
  expect_error(
    twoStageSpiDE(spe, condition = "condition", sigma = 30, min.cells = 10,
                 verbose = FALSE),
    "non-negative")
})

test_that("a patient with an NA covariate is dropped, not a crash", {
  # Before the fix, an NA in patient.covariates reached model.matrix(), which
  # silently drops that row; Xdes then had fewer rows than the (still
  # full-length) slope matrices, and .tau2DL()'s `X[ok, , drop = FALSE]`
  # errored with "(subscript) logical subscript too long".
  spe <- toy_spanorm_spe(n_genes = 20)
  spe$Age[spe$sample_id == "S1"] <- NA       # S1 is a Responder patient
  expect_message(
    r <- twoStageSpiDE(spe, condition = "condition", sigma = 30,
                       min.cells = 10, patient.covariates = "Age", fdr = 1,
                       verbose = TRUE),
    "dropping 1 patient.*Age")
  expect_s4_class(r, "SpiDEResults")
  expect_true(is.data.frame(results(r)))
})

test_that("an unknown patient.covariates name errors clearly", {
  spe <- toy_spanorm_spe(n_genes = 20)
  expect_error(
    twoStageSpiDE(spe, condition = "condition", sigma = 30, min.cells = 10,
                 patient.covariates = "not_a_column", verbose = FALSE),
    "not_a_column")
})

test_that("an unknown patient column errors clearly", {
  spe <- toy_spanorm_spe(n_genes = 20)
  expect_error(
    twoStageSpiDE(spe, condition = "condition", sigma = 30, min.cells = 10,
                 patient = "not_a_column", verbose = FALSE),
    "not_a_column")
})

test_that("a cell-varying patient.covariates column is refused, not silently collapsed", {
  # tapply(cd[[v]], pat, function(z) z[1]) would otherwise silently take the
  # first cell's value per patient, producing a wrong (arbitrary) covariate
  # value instead of erroring -- mirror the existing within-patient condition
  # check.
  spe <- toy_spanorm_spe(n_genes = 20)
  set.seed(2)
  spe$Age <- sample(c(30, 60), ncol(spe), replace = TRUE)  # varies per cell
  expect_error(
    twoStageSpiDE(spe, condition = "condition", sigma = 30, min.cells = 10,
                 patient.covariates = "Age", verbose = FALSE),
    "varies within patient")
})

test_that("too few patients per condition after an NA-covariate drop errors clearly", {
  spe <- toy_spanorm_spe(n_genes = 20)
  # S1, S3, S5 are the three Responder patients (default toy layout); drop
  # two of them, leaving only one Responder patient.
  spe$Age[spe$sample_id %in% c("S1", "S3")] <- NA
  expect_error(
    twoStageSpiDE(spe, condition = "condition", sigma = 30, min.cells = 10,
                 patient.covariates = "Age", fdr = 1, verbose = FALSE),
    "fewer than 2 patients")
})
