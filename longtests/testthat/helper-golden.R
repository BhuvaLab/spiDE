# Shared generator for the golden polish+inference regression fixtures.
#
# longtests/ runs as its own testthat directory (see longtests/testthat.R and
# longtests/testthat/helper-gpu.R's own comment), so this helper is not seen
# by tests/testthat/ and devtools::test() never picks the golden test up. It
# lives here rather than in tests/testthat/ because the full generator
# (four production-pipeline fits, one of them mixed-effects with a nested
# variance component and one a random-slope fit) takes ~4-5 min and tolerance
# 0 cannot be expected to hold across CI's other platforms/BLAS -- see
# longtests/testthat/test-polish-golden.R's header for the version-mismatch
# handling.
#
# Captures the polish + inference outputs of the production pipeline BEFORE
# the polish machinery moves to SpaNorm (design/plans/2026-09-25-polish-to-spanorm.md).
# longtests/testthat/_golden/make_golden_polish.R calls this once to (re)write
# longtests/testthat/_golden/golden_polish.rds; test-polish-golden.R calls it
# again on every run and compares at tolerance 0.
#
# Regenerate the .rds ONLY on a deliberate numerical change to the fit/polish/
# inference stages, never to make the SpaNorm move "pass". Every fit below is
# seeded immediately before the call: fitNB() subsamples cells for dispersion
# above a size threshold with no internal seed, so an unseeded fit is
# non-reproducible across R sessions (see CLAUDE.md, "Related:" paragraph
# under "Supplying psi changes the objective").
golden_polish_fits <- function() {
  snap <- function(res) {
    fits <- res@fits
    list(
      fits = lapply(fits, function(f) {
        list(
          alpha = f@alpha, psi = f@psi, loglik = f@loglik, tau2 = f@tau2,
          penalty = f@penalty, df = f@df, t_stat = f@t_stat, se = f@se,
          polish = f@polish
        )
      }),
      results = results(res)
    )
  }

  data("toySpiDE", package = "spiDE", envir = environment())

  # Condition mode, two bandwidths: the cross-bandwidth Cauchy combination
  # (.combineBandwidths()) must be exercised, not just a single-bandwidth fit.
  set.seed(20260925)
  toy <- spiDE(toySpiDE, condition = "condition", sigma = c(30, 50),
               BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)

  # Mixed-effects with a nested (sample x cell type) variance component,
  # df.method = "satterthwaite": @df is the named per-tested-column vector
  # (not the "between" scalar), and this pins the post-polish Satterthwaite
  # refresh in .polishSpiDEFit().
  set.seed(20260925)
  cl <- .toyClustered(sd_nested = 0.3)
  clustered <- spiDE(cl, condition = "condition", sigma = 30,
                      df.method = "satterthwaite",
                      BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)

  set.seed(20260925)
  niche <- spiDE(toySpiDE, condition = NULL, sigma = 30,
                 BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)

  # random = "slope": per-sample random slopes on the CellType:niche bases,
  # absorbed by sample (not by the nested logical -- see .absorbSpec()) in
  # both the polish and the blocked inference. Same construction as
  # tests/testthat/test-slope-absorb-wiring.R's .slopeWiringFit() (the
  # smallest fixture already known to exercise this path), run through the
  # full production pipeline (fitSpiDE -> polishSpiDE -> testSpiDE) instead
  # of stopping at the fit.
  set.seed(20260925)
  speSlope <- buildNiches(.toySPE(n_genes = 5, n_per = 40), sigma = 30)
  fitSlope <- fitSpiDE(speSlope, condition = "condition", sigma = 30,
                        random = "slope", re.maxit = 1L,
                        BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
  fitSlope <- polishSpiDE(fitSlope, speSlope,
                           BPPARAM = BiocParallel::SerialParam(),
                           verbose = FALSE)
  slope <- testSpiDE(fitSlope, spe = speSlope,
                      BPPARAM = BiocParallel::SerialParam())

  list(
    toy = snap(toy),
    clustered = snap(clustered),
    nichemode = snap(niche),
    slope = snap(slope),
    meta = list(
      SpaNorm = as.character(utils::packageVersion("SpaNorm")),
      spiDE = as.character(utils::packageVersion("spiDE")),
      seed = 20260925
    )
  )
}
