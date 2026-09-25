# The batched and per-gene polish engines moved to SpaNorm with polishNB(), and
# their parity tests (a batch of one against the per-gene engine, batch
# boundaries, every per-gene path in one batch, a failing gene, the shared
# factorisation, the tensor path) moved with them, to SpaNorm's
# tests/testthat/test-polishEngine.R. What stays here is spiDE's side: that
# the top-level entry point reaches both engines.

test_that("spiDE() can reach the reference polish engine", {
  # The spec keeps the per-gene engine in the tree "as the reference
  # implementation and the test oracle, reachable through engine = 'gene'"
  # (SpaNorm::polishNB(engine = "gene") since the move). polishSpiDE()
  # exposes it; spiDE() did not forward it, so the top-level entry point could
  # not reach the reference implementation at all -- which is why the
  # deflation triage had to call the three stages by hand.
  #
  # The assertion is on the engine actually used: the polish says which it
  # took, so the message is the observation. Equality of the two would not do
  # -- they agree to 5e-13 by construction, so a silently ignored argument
  # would pass.
  spe <- buildNiches(.toySPE(n_genes = 6), sigma = 20)
  # backend = "cpu": this is about which POLISH engine spiDE() reaches, and on
  # a GPU node the default "auto" would send the fit to the device, where
  # SpaNorm's NB fitter fails on this fixture (H100, job 28552404)
  args <- list(spe, condition = "condition", sigma = 20, random = "intercept",
               re.maxit = 2L, fdr = 1, verbose = TRUE, backend = "cpu")
  expect_message(do.call(spiDE, c(args, list(engine = "gene"))), "per gene")
  expect_message(do.call(spiDE, c(args, list(engine = "batch"))),
                 "in batches of")
})
