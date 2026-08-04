# Precomputed fits for the unit tests.
#
# Bioconductor requires `R CMD check --no-build-vignettes` to finish in under
# 10 minutes. The mixed-effects tests were 2,200 s between them -- each one
# alone exceeded the whole budget -- because every test re-fitted a toy dataset
# through the PQL loop, which re-fits every gene up to `re.maxit` times.
#
# Shrinking the fixtures barely helped: `.toyClustered`'s strata are ~27 cells,
# below the `re.min.cells = 100` floor, so `re.prop` subsampling is disabled and
# every iteration fits all cells regardless of `n_per`. Cutting cells by 62%
# bought only 28%.
#
# So the fits are done ONCE here and shipped. The split is deliberate:
#
#   * CONTRACT assertions -- slot presence, @df length and names, covtype
#     levels, finiteness -- run against these fixtures in tests/. They verify
#     the object shape the rest of the package depends on, which is what a
#     fast check should catch.
#
#   * NUMERICAL assertions -- df ~ S-2, tau2 recovery, anti-conservative
#     deflation -- stay as LIVE fits in longtests/. Freezing those would mean a
#     regression in the fitting code went undetected, which is the one thing
#     they exist to catch.
#
#   Rscript data-raw/make_test_fixtures.R
devtools::load_all(".", quiet = TRUE)

dst <- "inst/extdata/testfits"
dir.create(dst, showWarnings = FALSE, recursive = TRUE)

# keep these in step with longtests/, which refits the same configurations live
cfgs <- list(
  s16 = list(n_samples = 16, n_per = 30, n_genes = 20),
  s12 = list(n_samples = 12, n_per = 30, n_genes = 20)
)

for (nm in names(cfgs)) {
  k <- cfgs[[nm]]
  spe <- buildNiches(
    spiDE:::.toyClustered(n_samples = k$n_samples, n_per = k$n_per,
                          n_genes = k$n_genes, sd_patient = 0.7),
    sigma = 30)
  for (dfm in c("between", "satterthwaite")) {
    f <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                  df.method = dfm, verbose = FALSE)
    # testSpiDE() too: fitSpiDE() alone leaves @t_stat / @se empty, so a
    # fit-only fixture cannot support the assertion that names(@df) align to
    # colnames(@t_stat) -- the alignment IS the contract being checked.
    f <- testSpiDE(f, spe = spe, fdr = 1)
    p <- file.path(dst, sprintf("fit_%s_intercept_%s.rds", nm, dfm))
    saveRDS(fits(f)[[1]], p)
    message(sprintf("wrote %s (%.1f MB)", p, file.size(p) / 1e6))
  }
  if (nm == "s16") {
    f <- fitSpiDE(spe, "condition", sigma = 30, random = "slope",
                  df.method = "satterthwaite", verbose = FALSE)
    f <- testSpiDE(f, spe = spe, fdr = 1)
    p <- file.path(dst, sprintf("fit_%s_slope_satterthwaite.rds", nm))
    saveRDS(fits(f)[[1]], p)
    message(sprintf("wrote %s (%.1f MB)", p, file.size(p) / 1e6))
  }
}

# a fixed-effects fit, for the "random='none' has no RE slots" contract test
spe <- buildNiches(spiDE:::.toySPE(), sigma = 20)
f <- fitSpiDE(spe, "condition", sigma = 20, verbose = FALSE)
f <- testSpiDE(f, spe = spe, fdr = 1)
p <- file.path(dst, "fit_toyspe_none.rds")
saveRDS(fits(f)[[1]], p)
message(sprintf("wrote %s (%.1f MB)", p, file.size(p) / 1e6))

# random = "slope" on the same toy data, for the RE-slot contract test. Slope
# fits are the most expensive thing the short suite did.
f <- fitSpiDE(spe, "condition", sigma = 20, random = "slope",
              df.method = "between", verbose = FALSE)
f <- testSpiDE(f, spe = spe, fdr = 1)
p <- file.path(dst, "fit_toyspe_slope.rds")
saveRDS(fits(f)[[1]], p)
message(sprintf("wrote %s (%.1f MB)", p, file.size(p) / 1e6))

message(sprintf("\ntotal fixture size: %.1f MB",
                sum(file.size(list.files(dst, full.names = TRUE))) / 1e6))
