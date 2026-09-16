# The batched per-gene Newton must be the per-gene one, restructured. These
# tests are written against .polishGene() as the oracle: it stays in the tree
# for exactly that purpose, reachable through engine = "gene".
#
# The risk being tested is not the arithmetic -- it is per-gene state leaking
# between slices of a batch. Genes converge at different iterations, take
# different numbers of line-search halvings, restart independently and fail
# independently, and every one of those is a per-gene branch that batching
# turns into an index set.

# a small design in the exact shape .polishFit() hands the solver: a dense
# block plus 0/1 indicators partitioning the cells
toy_batch <- function(n = 240L, G = 6L, seed = 11L) {
  set.seed(seed)
  x <- scale(rnorm(n))[, 1]
  ct <- rep(c(1, 2), length.out = n)
  X <- cbind(`(Intercept)` = 1, CellTypeB = as.numeric(ct == 2), niche = x)
  grp <- rep(seq_len(G), length.out = n)
  Z <- matrix(0, n, G, dimnames = list(NULL, paste0("SampleCellType", seq_len(G))))
  Z[cbind(seq_len(n), grp)] <- 1
  W <- cbind(X, Z)
  list(W = W, pen = c(0, 0, 0, rep(1 / 0.05, G)),
       nested = c(rep(FALSE, 3), rep(TRUE, G)),
       ct_cols = c(TRUE, TRUE, FALSE, rep(FALSE, G)))
}

# run the oracle over a set of genes, one at a time
gene_by_gene <- function(Yb, d, A0, psi0, solver, ...) {
  out <- lapply(seq_len(nrow(Yb)), function(i) {
    spiDE:::.polishGene(as.numeric(Yb[i, ]), d$W, A0[i, ], psi0[[i]], d$pen,
                        solver, ct_cols = d$ct_cols, ...)
  })
  list(alpha = do.call(rbind, lapply(out, `[[`, "alpha")),
       psi = vapply(out, `[[`, numeric(1), "psi"),
       loglik = vapply(out, `[[`, numeric(1), "loglik"),
       iterations = vapply(out, `[[`, integer(1), "iterations"),
       restarted = vapply(out, `[[`, logical(1), "restarted"),
       capped = vapply(out, `[[`, logical(1), "capped"),
       singular = vapply(out, `[[`, logical(1), "singular"),
       psi_bound = vapply(out, `[[`, logical(1), "psi_bound"),
       polished = vapply(out, `[[`, logical(1), "polished"))
}

expect_same_fit <- function(b, g, tol = 1e-12, fields = c("alpha", "psi", "loglik")) {
  for (f in fields) expect_equal(b[[f]], g[[f]], tolerance = tol, ignore_attr = TRUE)
  for (f in c("restarted", "capped", "singular", "psi_bound", "polished")) {
    expect_identical(unname(b[[f]]), unname(g[[f]]), info = f)
  }
}

test_that("a batch of one reproduces .polishGene() exactly", {
  d <- toy_batch()
  solver <- spiDE:::.newtonSolver(d$W, d$pen, d$nested)
  set.seed(1)
  y <- matrix(rnbinom(ncol(d$W) * 0 + nrow(d$W), mu = 6, size = 3), nrow = 1)
  A0 <- matrix(0, 1, ncol(d$W)); A0[1, 1] <- log(mean(y))
  b <- spiDE:::.polishBatch(y, d$W, A0, 0.4, d$pen, solver, ct_cols = d$ct_cols)
  g <- gene_by_gene(y, d, A0, 0.4, solver)
  expect_same_fit(b, g)
})

test_that("a batch reproduces the same genes run one at a time", {
  d <- toy_batch()
  solver <- spiDE:::.newtonSolver(d$W, d$pen, d$nested)
  set.seed(2)
  B <- 13L
  Yb <- matrix(rnbinom(B * nrow(d$W), mu = rep(c(2, 8, 30), length.out = B), size = 2),
               nrow = B)
  A0 <- matrix(0, B, ncol(d$W)); A0[, 1] <- log(pmax(rowMeans(Yb), 0.1))
  psi0 <- rep(0.4, B)
  b <- spiDE:::.polishBatch(Yb, d$W, A0, psi0, d$pen, solver, ct_cols = d$ct_cols)
  g <- gene_by_gene(Yb, d, A0, psi0, solver)
  expect_same_fit(b, g)
})

test_that("the batch boundary does not move a gene's answer", {
  # The batched sibling of test-polish.R's gene-blocking invariance, at a
  # tolerance rather than at machine precision, and the difference is real:
  # mu = A %*% t(W) and the penalty term A^2 %*% pen go through BLAS, which
  # blocks a 13-row GEMM differently from a 1-row one, so the summation order
  # depends on the batch. The per-gene engine has no such dependence and its
  # exact-invariance test still holds for it.
  #
  # Measured here at ~2e-12 relative on psi, which propagates from a profile
  # optimum found on a mu that differs in its last digits. 1e-9 is the bar: far
  # inside anything that could move a call, far outside the noise.
  d <- toy_batch()
  solver <- spiDE:::.newtonSolver(d$W, d$pen, d$nested)
  set.seed(3)
  B <- 13L
  Yb <- matrix(rnbinom(B * nrow(d$W), mu = rep(c(3, 12), length.out = B), size = 2), nrow = B)
  A0 <- matrix(0, B, ncol(d$W)); A0[, 1] <- log(pmax(rowMeans(Yb), 0.1))
  full <- spiDE:::.polishBatch(Yb, d$W, A0, rep(0.4, B), d$pen, solver, ct_cols = d$ct_cols)
  for (bs in c(1L, 2L, 5L)) {
    idx <- split(seq_len(B), ceiling(seq_len(B) / bs))
    parts <- lapply(idx, function(ii)
      spiDE:::.polishBatch(Yb[ii, , drop = FALSE], d$W, A0[ii, , drop = FALSE],
                           rep(0.4, length(ii)), d$pen, solver, ct_cols = d$ct_cols))
    got <- list(alpha = do.call(rbind, lapply(parts, `[[`, "alpha")),
                psi = unlist(lapply(parts, `[[`, "psi")),
                loglik = unlist(lapply(parts, `[[`, "loglik")))
    expect_equal(got$alpha, full$alpha, tolerance = 1e-9, ignore_attr = TRUE,
                 info = paste("batch size", bs))
    expect_equal(got$psi, full$psi, tolerance = 1e-9, ignore_attr = TRUE)
    expect_equal(got$loglik, full$loglik, tolerance = 1e-9, ignore_attr = TRUE)
  }
})

test_that("one batch holds genes taking every different path", {
  # THE test: a converged gene, a degenerate start, an all-zero gene, a gene
  # that needs several halvings, and a gene whose psi sits on a bound, in one
  # batch, against the same five run singly. This is what catches per-gene
  # state leaking across slices.
  d <- toy_batch()
  solver <- spiDE:::.newtonSolver(d$W, d$pen, d$nested)
  set.seed(5)
  n <- nrow(d$W)
  Yb <- rbind(
    rnbinom(n, mu = 8, size = 3),        # ordinary
    rnbinom(n, mu = 40, size = 5),       # bright
    rep(0L, n),                          # all zero: psi runs to a bound
    rpois(n, lambda = 4),                # under-dispersed: the other bound
    rnbinom(n, mu = 2, size = 1)         # noisy, hard line search
  )
  B <- nrow(Yb)
  A0 <- matrix(0, B, ncol(d$W)); A0[, 1] <- log(pmax(rowMeans(Yb), 0.1))
  A0[2, ] <- -30                          # a degenerate start: forces a restart
  psi0 <- c(0.4, 0.2, 0.5, 0.3, 0.8)
  b <- spiDE:::.polishBatch(Yb, d$W, A0, psi0, d$pen, solver, ct_cols = d$ct_cols)
  g <- gene_by_gene(Yb, d, A0, psi0, solver)
  expect_same_fit(b, g)
  expect_identical(b$iterations, g$iterations)
  expect_true(any(g$restarted))           # the fixture must actually exercise it
  expect_true(any(g$psi_bound))
})

test_that("a gene that cannot be polished does not poison its batch", {
  # .waldCauchyBlock() fails a whole sub-batch when one gene's Cholesky fails,
  # and its only recourse is to tell the user to shrink cov.batch. The polish
  # must not acquire that failure mode. The warm path's documented fallback --
  # a non-finite start returns fitNB's own fit, not the sane start -- is a real
  # code path, so no stub solver is needed to reach it.
  d <- toy_batch()
  solver <- spiDE:::.newtonSolver(d$W, d$pen, d$nested)
  set.seed(6)
  B <- 4L
  Yb <- matrix(rnbinom(B * nrow(d$W), mu = 7, size = 3), nrow = B)
  A0 <- matrix(0, B, ncol(d$W)); A0[, 1] <- log(rowMeans(Yb))
  A0[2, 3] <- NaN
  psi0 <- rep(0.4, B)
  b <- spiDE:::.polishBatch(Yb, d$W, A0, psi0, d$pen, solver,
                            ct_cols = d$ct_cols, warm = TRUE)
  g <- gene_by_gene(Yb, d, A0, psi0, solver, warm = TRUE)
  expect_false(b$polished[2])
  expect_identical(b$alpha[2, ], A0[2, ])     # fitNB's fit kept, not the sane start
  expect_true(all(b$polished[-2]))            # the neighbours are unaffected
  expect_same_fit(b, g)
  # and identical to those three run without the failing gene present at all
  ok <- spiDE:::.polishBatch(Yb[-2, , drop = FALSE], d$W, A0[-2, , drop = FALSE],
                             psi0[-2], d$pen, solver, ct_cols = d$ct_cols, warm = TRUE)
  expect_equal(b$alpha[-2, ], ok$alpha, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that(".polishBatch reports what a shared factorisation would cost", {
  # Phase 2e's design question, made measurable. newton() keeps a list of
  # per-gene factorisations under a per-gene staleness counter. One shared
  # TENSOR factorisation cannot do that: it must refresh the whole active stack
  # whenever any gene in it is stale. How much of Phase 0b's memoisation that
  # discards is an empirical question about the staleness trajectory, not
  # something to reason about -- so the engine counts both.
  #
  #   factorisations       what the per-gene policy actually built
  #   factorisations_sync  what refreshing the whole active set would have built
  #
  # The second is a counterfactual and costs one integer per iteration.
  d <- toy_batch()
  set.seed(4)
  B <- 8L
  mu <- exp(d$W %*% c(1.2, 0.4, 0.3, rep(0, 6)))
  Yb <- matrix(rnbinom(B * nrow(d$W), mu = rep(as.numeric(mu), each = B),
                       size = 1 / 0.4), nrow = B)
  A0 <- matrix(0, B, ncol(d$W)); A0[, 1] <- log(pmax(rowMeans(Yb), 0.1))
  solver <- spiDE:::.newtonSolver(d$W, d$pen, d$nested)

  out <- spiDE:::.polishBatch(Yb, d$W, A0, rep(0.4, B), d$pen, solver,
                              ct_cols = d$ct_cols)
  nf <- attr(out, "factorisations")
  expect_type(nf, "integer")
  expect_named(nf, c("pergene", "sync"))

  # every gene is factorised at least once before its first step
  expect_gte(nf[["pergene"]], B)
  # and a shared factorisation can never build fewer than the per-gene policy
  expect_gte(nf[["sync"]], nf[["pergene"]])
})

test_that("spiDE() can reach the reference polish engine", {
  # The spec keeps .polishGene() in the tree "as the reference implementation
  # and the test oracle, reachable through engine = 'gene'". polishSpiDE()
  # exposes it; spiDE() did not forward it, so the top-level entry point could
  # not reach the reference implementation at all -- which is why the
  # deflation triage had to call the three stages by hand.
  #
  # The assertion is on the engine actually used: the polish says which it
  # took, so the message is the observation. Equality of the two would not do
  # -- they agree to 5e-13 by construction, so a silently ignored argument
  # would pass.
  spe <- buildNiches(.toySPE(n_genes = 6), sigma = 20)
  args <- list(spe, condition = "condition", sigma = 20, random = "intercept",
               re.maxit = 2L, fdr = 1, verbose = TRUE)
  expect_message(do.call(spiDE, c(args, list(engine = "gene"))), "per gene")
  expect_message(do.call(spiDE, c(args, list(engine = "batch"))),
                 "in batches of")
})
