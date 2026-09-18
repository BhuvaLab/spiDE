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

  # The three-way niche terms are bounded too. They are a between-patient
  # comparison OF a within-patient slope: each patient contributes one slope,
  # and more cells per patient sharpen it without creating more of them. On the
  # v11 arm they ran to a median df of 27,345 against a residual df of 76,347,
  # with cor(df, cells of the index compartment) = -0.722 -- the same inversion
  # as the two-way layer, which a genuinely within-patient quantity would not
  # show.
  dn <- ff@df[names(ff@df) %in% ff@coefmap$covariate[ct == "ResponseNiche"]]
  if (length(dn)) expect_lte(max(dn), (S - 2) * 1.05)
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

test_that("the bound survives the polish, which recomputes the df", {
  # The production path ALWAYS polishes, and .polishSpiDEFit() refreshes the
  # reference df at the reported penalty after the tau2 loop -- calling
  # .satterthwaiteDF() directly and assigning to @df. A bound applied only in
  # .fitNBmixed() is therefore computed and then overwritten, i.e. inert
  # exactly where it matters.
  #
  # This is the second time on this branch: .reprofilePsi() would have
  # overwritten the dispersion bisection at the last step the same way. Any
  # quantity the fit computes and the polish recomputes needs the correction at
  # BOTH sites, or at one site both call.
  spe <- buildNiches(spiDE:::.toyClustered(n_samples = 16, n_per = 30,
                                           n_genes = 8, sd_patient = 0.20),
                     sigma = 30)
  f <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                re.celltype = TRUE, df.method = "satterthwaite",
                verbose = FALSE, backend = "cpu")
  S <- sum(fits(f)[[1]]@re_group == "SampleInt", na.rm = TRUE)

  fp <- polishSpiDE(f, spe, verbose = FALSE)
  ff <- fits(fp)[[1]]
  ct <- as.character(ff@covtype)
  d <- ff@df[names(ff@df) %in% ff@coefmap$covariate[ct == "ResponseCellType"]]
  expect_gt(length(d), 0L)
  expect_lte(max(d), (S - 2) * 1.05)
})

test_that("the cap is per compartment, not a single flat value", {
  # A flat cap at S - 2 ties every violating compartment together, discarding a
  # real difference: a compartment present in all 16 patients carries more
  # information than one present in 9. The contributing-patient count per
  # compartment is recoverable from the nested (patient x cell type) groups --
  # there is one such column per NON-EMPTY pair, so counting the groups of a
  # compartment counts its patients.
  set.seed(11)
  spe <- buildNiches(spiDE:::.toyClustered(n_samples = 16, n_per = 60,
                                           n_genes = 8, sd_patient = 0.20),
                     sigma = 30)
  # remove one cell type from half the patients entirely, so it is present in
  # far fewer patients than the others
  types <- levels(factor(spe$cell_type))
  drop_in <- unique(spe$sample_id)[1:8]
  gone <- spe$cell_type == types[1] & spe$sample_id %in% drop_in
  spe <- spe[, !gone]

  f <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                re.celltype = TRUE, df.method = "satterthwaite",
                verbose = FALSE, backend = "cpu")
  ff <- fits(f)[[1]]
  ct <- as.character(ff@covtype)
  rc <- which(ct == "ResponseCellType")
  d <- ff@df[names(ff@df) %in% ff@coefmap$covariate[rc]]
  lab <- ff@coefmap$index[rc][match(names(d), ff@coefmap$covariate[rc])]

  # patients per compartment, straight from the nested group names
  nz <- which(!is.na(ff@re_group) & ff@re_group == "SampleCellTypeInt")
  gnm <- colnames(ff@W)[nz]
  per <- vapply(lab, function(k) sum(endsWith(gnm, paste0(".", k))), numeric(1))
  expect_true(all(per > 0))
  expect_gt(max(per) - min(per), 2)          # the fixture must be unbalanced

  # the thinned compartment must be capped BELOW the others
  thin <- which.min(per)
  expect_lt(d[thin], max(d))
  # and no compartment may exceed its own patients
  expect_true(all(d <= pmax(per - 2, 1) * 1.05 + 1e-8),
              info = paste(sprintf("%s: df %.1f vs patients %d", lab, d, per),
                           collapse = " | "))
})

test_that(".patientsPerTested matches labels that contain spaces", {
  # The nested columns carry the RAW cell-type label and coefmap$index the
  # sanitised one, so "B cell" vs "B.cell" silently missed on the real cohort:
  # three of thirteen compartments resolved to NA and quietly fell back to the
  # flat S - 2. The bound stayed valid, which is why nothing failed -- the
  # refinement just stopped refining.
  W <- matrix(0, 6, 4,
              dimnames = list(NULL, c("CellTypeB.cell:ResponseResponder",
                                      "SampleCellTypeP1.B cell",
                                      "SampleCellTypeP2.B cell",
                                      "SampleCellTypeP3.T cell")))
  W[1:2, 2] <- 1; W[3:4, 3] <- 1; W[5:6, 4] <- 1
  re_group <- c(NA, rep("SampleCellTypeInt", 3))
  coefmap <- data.frame(covariate = colnames(W),
                        index = c("B.cell", NA, NA, NA),
                        stringsAsFactors = FALSE)
  got <- spiDE:::.patientsPerTested(W, re_group, coefmap, tested = 1L)
  expect_equal(got, 2)          # B cell is in P1 and P2, not P3
})

# ---------------------------------------------------------------------------
# The default reference df (2026-09-18, phase 6).
#
# Satterthwaite was the default because it promised a per-column df that
# distinguished between-patient contrasts from within-patient ones. On this
# design it does not: every tested column is a Response term, the patient
# bound caps all of them at their own compartment's n - 2, and on the cohort
# design the bound binds on 168 of 168 columns. A per-column vector that is
# the bound everywhere carries nothing the scalar does not, while costing a
# variance-parameter covariance per gene and inviting callers to read
# structure into it. "between" is the default; satterthwaite stays available
# and is now bounded too, so the arm cannot mislead.
# ---------------------------------------------------------------------------

test_that("df.method defaults to between at every public entry point", {
  # source-level: all three signatures must agree, or a caller that goes
  # through spiDE() gets a different default from one that calls fitSpiDE()
  # the methods dispatch on spe = "ANY", and S4 rematching hides the real
  # formals inside .local(), so deparse the whole method rather than reading
  # formals() off the generic
  for (fn in c("fitSpiDE", "spiDE")) {
    # deparse wraps long signatures, so collapse runs of whitespace before
    # matching or the default reads as c("between",  "satterthwaite")
    src <- gsub("[[:space:]]+", " ", paste(deparse(getMethod(fn, "ANY")), collapse = " "))
    expect_true(grepl('df.method = c("between", "satterthwaite")', src, fixed = TRUE),
                info = paste(fn, "signature:", 
                             paste(regmatches(src, gregexpr("df\\.method = [^,)]*", src))[[1]],
                                   collapse = " ;; ")))
  }
  expect_identical(eval(formals(spiDE:::.fitNBmixed)$df.method), "between")
})

test_that("a mixed fit returns the scalar df by default and the vector on request", {
  spe <- buildNiches(spiDE:::.toyClustered(n_samples = 8, n_per = 30, n_genes = 6),
                     sigma = 30, verbose = FALSE)
  f_def <- fits(fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                         verbose = FALSE))[[1]]
  f_sat <- fits(fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                         df.method = "satterthwaite", verbose = FALSE))[[1]]
  n_patients <- length(unique(as.character(spe$sample_id)))
  expect_length(f_def@df, 1L)
  expect_equal(unname(f_def@df), n_patients - 2)
  # and the opt-in arm still produces one per tested column, all bounded
  expect_gt(length(f_sat@df), 1L)
  expect_true(all(f_sat@df <= n_patients - 2))
  expect_true(all(is.finite(f_sat@df)))
})
