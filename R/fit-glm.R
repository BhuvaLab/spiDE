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
  a[, seq_len(des$npat)] <- log(lm + 1e-3)
  if (!is.null(offset)) {
    mo <- if (is.matrix(offset)) {
      t(rowsum(t(offset), des$patient, reorder = TRUE)) / rep(npc, each = nrow(Yk))
    } else {
      matrix(tapply(offset, des$patient, mean), nrow(Yk), des$npat, byrow = TRUE)
    }
    a[, seq_len(des$npat)] <- a[, seq_len(des$npat)] - mo
  }
  a
}

#' @importFrom SpaNorm polishNB
.fitIndexGLM <- function(Yk, des, offset = NULL, BPPARAM = BiocParallel::SerialParam()) {
  a0 <- .glmStart(Yk, des, offset)
  rownames(a0) <- rownames(Yk)   # polishNB keys its per-gene diagnostics on these
  fit <- SpaNorm::polishNB(Yk, des$W, a0, rep(1, nrow(Yk)), lambda.a = des$pen,
                           offset = offset, absorb = des$absorb, start.cols = des$absorb,
                           psi.method = "profile", BPPARAM = BPPARAM)
  rownames(fit$alpha) <- rownames(Yk)
  names(fit$psi) <- rownames(Yk)
  fit
}

# Per-gene working quantities at the converged fit: weights w and working
# residuals r of the NB score, w = mu / (1 + psi mu), r = (y - mu) / (1 + psi mu).
.workingWR <- function(y, W, alpha, psi, offset = 0) {
  mu <- exp(as.numeric(W %*% alpha) + offset)
  list(w = mu / (1 + psi * mu), r = (y - mu) / (1 + psi * mu), mu = mu)
}

# Split genes into chunks and run FUN over them with BiocParallel, one BLAS
# thread per worker: forked workers inherit the parent's BLAS thread count, and
# n workers x m threads on n cores oversubscribes the node (measured 9x slower
# on the mixed model's polish stage; see research/fdr-ordering/FINDINGS.md,
# 2026-09-10).
.bpGenes <- function(G, FUN, BPPARAM = BiocParallel::SerialParam()) {
  nw <- max(1L, BiocParallel::bpnworkers(BPPARAM))
  chunks <- split(seq_len(G), cut(seq_len(G), min(G, nw * 4L), labels = FALSE))
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
