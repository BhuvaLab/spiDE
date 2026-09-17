# A between-patient contrast cannot have more degrees of freedom than the
# patients that carry it.
#
# The v11 arm violates this badly and in a specific direction. Measured on
# data/Niche30Fitv11.rds (55 patients, S - 2 = 53):
#
#   Tumor        27,764 cells    df    306
#   Fibroblast   13,380 cells    df    714
#   ...
#   Smooth.m.     1,073 cells    df 34,544
#   NK             785 cells     df 25,225
#
#   Spearman cor(df, cells) = -0.978 over 13 compartments
#
# So the df is NOT "cell-level" -- a cell-level df would rank Tumor top. It
# runs INVERSELY to the evidence: the rarer the compartment, the larger the
# reference df, which is anti-conservative exactly where the data is thinnest.
# The mechanism is that a rare compartment's (patient x cell type) groups hold
# few cells, their random effects shrink hard, the variance-component gradients
# collapse, varse2 shrinks, and df = 2 vjj^2 / varse2 explodes.
#
# The fix is a BOUND rather than a recalibration: whatever the Satterthwaite
# approximation returns, a contrast that is constant within a patient is
# identified by patients, so its df is bounded by the patients contributing to
# it. The existing clamp is at ncells, which is why 34,544 survives.

test_that("a between-patient contrast's df is bounded by the patients", {
  spe <- buildNiches(spiDE:::.toyClustered(n_samples = 16, n_per = 30,
                                           n_genes = 10, sd_patient = 0.20),
                     sigma = 30)
  f <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                re.celltype = TRUE, df.method = "satterthwaite",
                verbose = FALSE, backend = "cpu")
  ff <- fits(f)[[1]]
  ct <- as.character(ff@covtype)
  S <- sum(ff@re_group == "SampleInt", na.rm = TRUE)
  d <- ff@df[names(ff@df) %in% ff@coefmap$covariate[ct == "ResponseCellType"]]
  expect_gt(length(d), 0L)
  # the bound, with a little slack for the Satterthwaite approximation itself
  expect_lte(max(d), (S - 2) * 1.05)

  # and the within-patient layer must NOT be capped with it: the niche varies
  # inside a patient, so it legitimately carries more information
  dn <- ff@df[names(ff@df) %in% ff@coefmap$covariate[ct == "ResponseNiche"]]
  if (length(dn)) expect_gt(stats::median(dn), stats::median(d))
})

test_that("a rarer compartment does not earn MORE degrees of freedom", {
  # the inversion, reproduced small: thin one cell type down and check its df
  # does not overtake the abundant one's. This is the shape of the cohort
  # defect (cor(df, cells) = -0.978), and it is the assertion a fix has to
  # satisfy -- capping alone would pass the bound above while leaving the
  # ordering upside down.
  set.seed(4)
  spe <- buildNiches(spiDE:::.toyClustered(n_samples = 16, n_per = 60,
                                           n_genes = 10, sd_patient = 0.20),
                     sigma = 30)
  keep <- !(spe$cell_type == levels(factor(spe$cell_type))[1]) |
    (stats::runif(ncol(spe)) < 0.15)     # thin cell type 1 to ~15%
  spe <- spe[, keep]
  f <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                re.celltype = TRUE, df.method = "satterthwaite",
                verbose = FALSE, backend = "cpu")
  ff <- fits(f)[[1]]
  ct <- as.character(ff@covtype)
  rc <- which(ct == "ResponseCellType")
  d <- ff@df[names(ff@df) %in% ff@coefmap$covariate[rc]]
  lab <- ff@coefmap$index[rc][match(names(d), ff@coefmap$covariate[rc])]
  n <- vapply(lab, function(k) sum(spe$cell_type == k), numeric(1))
  ok <- !is.na(n) & n > 0
  expect_gte(sum(ok), 2L)
  # The inversion must be gone. Flattening counts as gone: a conservative bound
  # ties the capped columns at the patient-level df, which loses the legitimate
  # ordering (a compartment present in 55 patients should out-rank one present
  # in 30) but is not WRONG the way an inversion is. Refining the cap to the
  # per-compartment patient count is the better answer and is a separate step.
  rho <- suppressWarnings(stats::cor(d[ok], n[ok], method = "spearman"))
  expect_true(stats::sd(d[ok]) == 0 || is.na(rho) || rho > -0.5,
              info = sprintf("Spearman cor(df, cells) = %.3f, sd(df) = %.3g",
                             rho, stats::sd(d[ok])))
})
