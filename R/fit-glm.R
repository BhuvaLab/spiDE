# One negative binomial GLM per gene for one index cell type, fitted with
# SpaNorm::polishNB(): damped Newton on each gene's own penalised likelihood,
# the patient intercepts absorbed by Schur complement, the dispersion by
# profile maximum likelihood at the converged mean. Genes are independent, so
# polishNB() may block and parallelise over them (BPPARAM); nothing here is
# shared across genes.

# A sane start per gene: each patient's log mean count (net of any offset) in
# the intercepts, zeros elsewhere.
.glmStart <- function(Yk, des, offset = NULL) {
  npc <- as.numeric(table(factor(des$patient, levels = seq_len(des$npat))))
  # per-patient sums as a product with a sparse indicator, so a sparse or
  # DelayedArray count matrix is never densified whole (genes x cells of a
  # WTA tumour compartment is several GB dense); only genes x patients is
  ind <- Matrix::sparseMatrix(i = seq_along(des$patient), j = des$patient, x = 1,
                              dims = c(length(des$patient), des$npat))
  lm <- as.matrix(Yk %*% ind) / rep(npc, each = nrow(Yk))
  a <- matrix(0, nrow(Yk), ncol(des$W))
  a[, des$intercept_cols] <- log(lm + 1e-3)
  if (!is.null(offset)) {
    mo <- if (is.matrix(offset)) {
      t(rowsum(t(offset), des$patient, reorder = TRUE)) / rep(npc, each = nrow(Yk))
    } else {
      matrix(tapply(offset, des$patient, mean), nrow(Yk), des$npat, byrow = TRUE)
    }
    a[, des$intercept_cols] <- a[, des$intercept_cols] - mo
  }
  a
}

#' @importFrom SpaNorm polishNB
.fitIndexGLM <- function(Yk, des, offset = NULL, BPPARAM = BiocParallel::SerialParam()) {
  a0 <- .glmStart(Yk, des, offset)
  rownames(a0) <- rownames(Yk)   # polishNB keys its per-gene diagnostics on these
  fit <- SpaNorm::polishNB(Yk, des$W, a0, rep(1, nrow(Yk)), lambda.a = des$pen,
                           offset = offset, absorb = des$absorb, start.cols = des$start,
                           psi.method = "profile", BPPARAM = BPPARAM)
  rownames(fit$alpha) <- rownames(Yk)
  names(fit$psi) <- rownames(Yk)
  .fitStatus(fit, Yk, des, offset, BPPARAM)
}

# polishNB() returns a gene it could not fit (singular, non-finite) at its
# start (patient log means, slopes 0) with polished = FALSE, and a gene whose
# profile dispersion optimum lies on the search bound with the input
# dispersion (1 here) and a mean converged at it. The first would be tested as
# if it had converged, so it drops out (NA). The second takes the bound its
# likelihood prefers (near-Poisson at the lower bound, a very overdispersed
# gene at the upper), since the constrained maximum is on that bound, and its
# mean is refitted at that dispersion. Counts are kept in fit$status.
.fitStatus <- function(fit, Yk, des, offset, BPPARAM = BiocParallel::SerialParam(),
                       psi.range = c(1e-3, 1e3)) {
  pol <- fit$polish
  bad <- if (is.null(pol$polished)) rep(FALSE, nrow(Yk)) else !pol$polished
  bad <- bad | !is.finite(fit$psi) | !apply(is.finite(fit$alpha), 1, all)
  fit$alpha[bad, ] <- NA_real_
  fit$psi[bad] <- NA_real_
  bnd <- which(!bad & (if (is.null(pol$psi_bound)) FALSE else pol$psi_bound))
  for (g in bnd) {
    o <- if (is.null(offset)) 0 else if (is.matrix(offset)) offset[g, ] else offset
    mu <- exp(.linPred(des$W, fit$alpha[g, ]) + o)
    y <- as.numeric(Yk[g, ])
    ll <- vapply(psi.range, function(p) sum(stats::dnbinom(y, size = 1 / p, mu = mu, log = TRUE)), numeric(1))
    fit$psi[g] <- psi.range[which.max(ll)]
  }
  if (length(bnd)) {
    ob <- if (is.matrix(offset)) offset[bnd, , drop = FALSE] else offset
    ref <- SpaNorm::polishNB(Yk[bnd, , drop = FALSE], des$W, fit$alpha[bnd, , drop = FALSE], fit$psi[bnd],
                             lambda.a = des$pen, offset = ob, absorb = des$absorb, start.cols = des$start,
                             psi.method = "fixed", BPPARAM = BPPARAM)
    ok <- if (is.null(ref$polish$polished)) rep(TRUE, length(bnd)) else ref$polish$polished
    ok <- ok & apply(is.finite(ref$alpha), 1, all)
    fit$alpha[bnd[ok], ] <- ref$alpha[ok, , drop = FALSE]
    fit$alpha[bnd[!ok], ] <- NA_real_
    fit$psi[bnd[!ok]] <- NA_real_
    bad[bnd[!ok]] <- TRUE
    bnd <- bnd[ok]
  }
  fit$status <- c(genes = nrow(Yk), not_fitted = sum(bad), psi_at_bound = length(bnd))
  fit
}

# Per-gene working quantities at the converged fit: weights w and working
# residuals r of the NB score, w = mu / (1 + psi mu), r = (y - mu) / (1 + psi mu).
.workingWR <- function(y, W, alpha, psi, offset = 0) {
  mu <- exp(.linPred(W, alpha) + offset)
  list(w = mu / (1 + psi * mu), r = (y - mu) / (1 + psi * mu), mu = mu)
}

# Split genes into chunks and run FUN over them with BiocParallel, one BLAS
# thread per worker: forked workers inherit the parent's BLAS thread count, and
# n workers x m threads on n cores oversubscribes the node (measured 9x slower
# on the mixed model's polish stage; see research/fdr-ordering/FINDINGS.md,
# 2026-09-10).
.bpGenes <- function(G, FUN, BPPARAM = BiocParallel::SerialParam()) {
  nw <- max(1L, BiocParallel::bpnworkers(BPPARAM))
  nch <- min(G, nw * 4L)
  # cut() refuses a single interval, so one chunk (e.g. a one-gene index type) is built directly
  chunks <- if (nch <= 1L) list(seq_len(G)) else split(seq_len(G), cut(seq_len(G), nch, labels = FALSE))
  single <- function(idx) {
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      old <- RhpcBLASctl::blas_get_num_procs()
      RhpcBLASctl::blas_set_num_threads(1L)
      on.exit(RhpcBLASctl::blas_set_num_threads(old), add = TRUE)
    }
    FUN(idx)
  }
  res <- if (nw > 1L) BiocParallel::bplapply(chunks, single, BPPARAM = BPPARAM) else lapply(chunks, single)
  res
}
