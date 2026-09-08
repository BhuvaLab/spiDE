# The polish stage, the variance components, and retiring the two-stage estimator

Date: 2026-09-08. Status: **in progress** (approved).

## 1. The pipeline becomes fit -> polish -> test -> gsea

`fitSpiDE()` fits the shared model and nothing more: the per-gene convergence
step and its arguments (`converge`, `converge.maxit`, `converge.tol`,
`polish.psi`) leave it. `polishSpiDE()` is the stage the user runs, or skips,
according to their data:

1. per-gene damped Newton convergence at the fit's penalty (as now);
2. the dispersion rule, `psi = c("profile", "moderated")` — profile is the
   default again, because the moderated rule keeps whatever the shared fit
   left (fifteen times the converged value on the clustered fixture,
   `fdr-ordering/FINDINGS.md` 2026-09-08) and the weights of step 3 need a
   dispersion consistent with the converged mean;
3. **the variance components re-estimated from the converged fit** (mixed
   fits only, `tau2 = TRUE`): one Schall step with the polished coefficients
   in the numerator and the gene-averaged working weights at the polished
   mean and dispersion, then a re-polish at the new penalty, iterated to
   `tau2.tol` within `tau2.maxit`; `@tau2`, `@penalty` and the Satterthwaite
   `@df` are refreshed;
4. inference cleared, so `testSpiDE()` runs again.

`spiDE()` gains `polish = TRUE` in place of the convergence arguments and
calls the stage with its defaults. `testSpiDE()` and `spiGSEA()` are
unchanged. The test for step 3 is the long test's window (0.3–0.75 on the
clustered fixture that reads 10 today; the 2026-08-04 fixture stays near
0.43), added to the unit tests on a smaller fixture as well.

## 2. The two-stage estimator leaves the package

`twoStageSpiDE()` and its stages move to a standalone package in the research
tree, `research/twostage/` (`spiDEtwostage`), which imports spiDE for the
niches and checkers and returns its own result object (a list with the tidy
`results` table, `diagnostics` and `sigma`), so spiDE's `SpiDEResults` loses
the `diagnostics` slot and the validity relaxation that allowed empty fits.
The package need not be polished: the code and its tests are ported, the
vignette moves with it, the research harness's two-stage arms point at it.
Everything else in spiDE that mentions the estimator goes: the vignette, the
model and calibration vignettes' sections, `_pkgdown.yml`, NEWS (a removal
note), CLAUDE.md, the palette keeps its `twostage` colour for the research
site's archived rows.

## 3. Order

1. Polish stage and variance components, test-first, suite green.
2. Two-stage port and removal, suite green.
3. Docs, 0.99.19, R CMD check, merge to main, delete the branch.
4. Reruns resubmitted under the new pipeline (the moderated-psi benchmark
   array of 2026-09-08 was cancelled as no longer the shipped pipeline).
