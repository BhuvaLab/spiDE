# Shared generator for the golden polish+inference regression fixtures.
#
# Captures the polish + inference outputs of the production pipeline BEFORE
# the polish machinery moves to SpaNorm (design/plans/2026-09-25-polish-to-spanorm.md).
# tests/testthat/_golden/make_golden_polish.R calls this once to (re)write
# tests/testthat/_golden/golden_polish.rds; test-polish-golden.R calls it again
# on every test run and compares at tolerance 0.
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

  # A single bandwidth (30): sigma = c(30, 50) together pushed the three fits
  # this generator runs to ~8 min wall clock (measured 2026-09-25), over the
  # ~3 min budget, so only the cheaper single-bandwidth toy fit is kept.
  set.seed(20260925)
  toy <- spiDE(toySpiDE, condition = "condition", sigma = 30,
               BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)

  set.seed(20260925)
  cl <- .toyClustered(sd_nested = 0.3)
  clustered <- spiDE(cl, condition = "condition", sigma = 30,
                      BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)

  set.seed(20260925)
  niche <- spiDE(toySpiDE, condition = NULL, sigma = 30,
                 BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)

  list(
    toy = snap(toy),
    clustered = snap(clustered),
    nichemode = snap(niche),
    meta = list(
      SpaNorm = as.character(utils::packageVersion("SpaNorm")),
      spiDE = as.character(utils::packageVersion("spiDE")),
      seed = 20260925
    )
  )
}
