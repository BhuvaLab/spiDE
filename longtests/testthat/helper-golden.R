# Shared generator for the golden polish+inference regression fixtures.
#
# longtests/ runs as its own testthat directory (see longtests/testthat.R and
# longtests/testthat/helper-gpu.R's own comment), so this helper is not seen
# by tests/testthat/ and devtools::test() never picks the golden test up. It
# lives here rather than in tests/testthat/ because the full generator (four
# production-pipeline fits, one of them mixed-effects with a nested variance
# component and one a random-slope fit) takes ~4-5 min; the gate itself is
# opt-in (SPIDE_RUN_GOLDEN=true) -- see longtests/testthat/test-polish-golden.R's
# header for why (tolerance 0 across BLAS kernels) and for the version-mismatch
# handling. Every call below pins `backend = "cpu"` (spiDE's `spiDE()`/
# `fitSpiDE()`/`testSpiDE()` default to `backend = "auto"`; `polishSpiDE()`
# already defaults to `"cpu"`), because the RDS was captured on this node
# where `"auto"` happened to resolve to the CPU -- running the gate on a GPU
# host must reproduce that, not whatever `"auto"` resolves to there.
#
# Captures the polish + inference outputs of the production pipeline BEFORE
# the polish machinery moves to SpaNorm (design/plans/2026-09-25-polish-to-spanorm.md).
# longtests/testthat/_golden/make_golden_polish.R calls this once to (re)write
# longtests/testthat/_golden/golden_polish.rds; test-polish-golden.R calls it
# again on every run and compares at tolerance 0.
#
# What this pins, precisely (do not read more into it): per `SpiDEFit`,
# `alpha`/`psi`/`loglik`/`tau2`/`penalty`/`df`/`t_stat`/`se`/`@polish`; per
# config, `results(res)` at the *default* `fdr = 0.05` -- which holds only 0
# or 1 row across these four configs, so that clause is a near-vacuous check
# of the FDR cascade, not a real one. NOT snapshotted: `p.combined.pos/neg`,
# `gene.weights`, `p.cauchy.pos/neg` -- the two-bandwidth Cauchy combination
# runs (it feeds `results()`), but nothing here compares its own outputs.
# NOT exercised by any of the four configs: a per-gene restart, a `psi_bound`
# clamp, or a singular-gene path (0 restarted and 0 psi_bound genes in every
# config); `random = "none"`, `psi = "moderated"`, or `engine = "gene"`. Those
# paths are pinned only by code-identity review, not by this golden.
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

  # Condition mode, two bandwidths, so the cross-bandwidth Cauchy combination
  # (.combineBandwidths()) runs rather than a single-bandwidth fit. Its own
  # outputs (p.combined.pos/neg, gene.weights) are not in the snapshot -- see
  # the header -- so this config exercises that code path without pinning it.
  set.seed(20260925)
  toy <- spiDE(toySpiDE, condition = "condition", sigma = c(30, 50),
               backend = "cpu",
               BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)

  # Mixed-effects with a nested (sample x cell type) variance component,
  # df.method = "satterthwaite": @df is the named per-tested-column vector
  # (not the "between" scalar), and this pins the post-polish Satterthwaite
  # refresh in .polishSpiDEFit().
  set.seed(20260925)
  cl <- .toyClustered(sd_nested = 0.3)
  clustered <- spiDE(cl, condition = "condition", sigma = 30,
                      df.method = "satterthwaite", backend = "cpu",
                      BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)

  set.seed(20260925)
  niche <- spiDE(toySpiDE, condition = NULL, sigma = 30,
                 backend = "cpu",
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
                        random = "slope", re.maxit = 1L, backend = "cpu",
                        BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
  fitSlope <- polishSpiDE(fitSlope, speSlope,
                           BPPARAM = BiocParallel::SerialParam(),
                           verbose = FALSE)
  slope <- testSpiDE(fitSlope, spe = speSlope, backend = "cpu",
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
