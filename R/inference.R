# Wald inference + within-gene p-value combination (Brown's method or Cauchy/
# ACAT, per the `combine` argument). For each gene the working weights, Wald
# covariance,
# standard errors and t-statistics of the Response / ResponseNiche coefficients
# are computed, then per-covariate p-values are combined across correlated
# ResponseNiche covariates (Cauchy/ACAT by default, or Brown's method) into a
# gene-level and a per-index-cell-type p-value, separately for up and down.
# Reproduces the per-gene loop in batch_nichede_v9.R. Each gene is independent,
# so the loop is chunked into gene-blocks and parallelised with BiocParallel;
# the counts are realised one block at a time (DelayedArray-friendly).

#' Partition gene indices into blocks
#' @noRd
.chunkGenes <- function(n, block.size = NULL) {
  if (is.null(block.size) || block.size >= n) {
    return(list(seq_len(n)))
  }
  block.size <- max(1L, as.integer(block.size))
  split(seq_len(n), ceiling(seq_len(n) / block.size))
}

#' Apply a t (or normal) tail to a statistic matrix, with per-column df
#'
#' \code{df} may be NULL (normal reference), a scalar (one df for all columns),
#' a length-\code{ncol(tmat)} vector (a df per tested column, broadcast down the
#' rows), or a matrix matching \code{tmat}. \code{tmat} may be a genes x k matrix
#' or a length-k vector (a single gene). Centralises the reference-distribution
#' logic so both the inference and FDR stages consume \code{@df} identically.
#' @noRd
.ptByCol <- function(tmat, df, lower.tail = TRUE) {
  if (is.null(df)) {
    return(stats::pnorm(tmat, lower.tail = lower.tail))
  }
  if (length(df) == 1L) {
    return(stats::pt(tmat, df = df, lower.tail = lower.tail))
  }
  dfx <- if (is.matrix(tmat) && !is.matrix(df)) {
    matrix(df, nrow(tmat), ncol(tmat), byrow = TRUE)
  } else {
    df               # vector tmat with length-k df, or df already a matrix
  }
  stats::pt(tmat, df = dfx, lower.tail = lower.tail)
}

#' Per-gene NB log-likelihood, computed gene-block-wise
#'
#' Sums the negative-binomial log-likelihood per gene without ever densifying
#' the whole counts matrix (genes are chunked, so only a block of \code{Y} is
#' realised at a time). Used for the cross-bandwidth Cauchy gene weights.
#'
#' @param Y counts (genes x cells; dense, sparse or DelayedArray).
#' @param alpha the fitted coefficients (genes x p).
#' @param W the design (cells x p).
#' @param psi per-gene NB dispersion.
#' @param block.size genes per block.
#' @return a numeric vector of per-gene log-likelihoods.
#' @importFrom stats dnbinom
#' @noRd
.blockLoglik <- function(Y, alpha, W, psi, block.size = 2000L) {
  ng <- nrow(alpha)
  blocks <- .chunkGenes(ng, block.size)
  ll <- numeric(ng)
  for (gi in blocks) {
    Yb <- as.matrix(Y[gi, , drop = FALSE])
    mub <- SpaNorm::calculateMu(rep(0, length(gi)),
                                alpha[gi, , drop = FALSE], W)
    ll[gi] <- rowSums(stats::dnbinom(Yb, mu = mub, size = 1 / psi[gi],
                                     log = TRUE))
  }
  ll
}

#' Wald inference + within-gene p-value combination for a single gene
#'
#' Computes per-column Wald t-statistics/SEs, then combines the correlated
#' ResponseNiche p-values into a gene-level and per-index-cell-type p-value
#' (separately per direction). The combiner is Brown's method (correlation-aware,
#' via \code{poolr}) or the Cauchy combination test (ACAT, correlation-agnostic).
#'
#' @param alpha_g numeric coefficients for the Response/ResponseNiche columns.
#' @param Wsub the Response/ResponseNiche sub-design (cells x k).
#' @param wt_g numeric working weights for this gene (length ncells).
#' @param psi_g the variance scale for the SEs: the NB dispersion for a
#'   fixed-effects fit, or the working (Pearson) dispersion for a mixed fit.
#' @param cov_niche logical marking the ResponseNiche columns within the subset.
#' @param index_ct index cell type of each ResponseNiche column.
#' @param uniq_index the unique index cell types (column order of the output).
#' @param combine one of "cauchy" (Cauchy/ACAT combination, default) or "brown"
#'   (Brown's method) for the within-gene combination.
#' @param W_full the full design (fixed + random columns); NULL for the
#'   fixed-effects fit, in which case the covariance is formed from \code{Wsub}
#'   alone and a normal reference is used.
#' @param penalty per-column ridge penalty aligned to \code{W_full} (the fitted
#'   \code{lambda.a}); only used when \code{W_full} is supplied.
#' @param sel integer positions of the tested (Response/ResponseNiche) columns
#'   within \code{W_full}; only used when \code{W_full} is supplied.
#' @param df residual degrees of freedom for the t reference; only used when
#'   \code{W_full} is supplied (else a normal reference).
#' @return a list with t_stat, se (length k) and p.pos, p.neg (length
#'   1 + n_index: gene-level then per-index-cell-type).
#' @importFrom stats pnorm pt cov2cor setNames
#' @noRd
.waldBrownGene <- function(alpha_g, Wsub, wt_g, psi_g, cov_niche, index_ct,
                           uniq_index, combine = c("cauchy", "brown"),
                           W_full = NULL, penalty = NULL,
                           sel = NULL, df = NULL) {
  combine <- match.arg(combine)
  if (is.null(W_full)) {
    # fixed-effects fit: covariance from the tested sub-design, normal reference
    varcov <- SpaNorm::invert_mat(crossprod(Wsub * wt_g, Wsub))
  } else {
    # mixed-effects fit: full penalised information (X'WX + Lambda), then the
    # fixed-effect (Response/ResponseNiche) block -> larger, between-sample SEs
    info <- crossprod(W_full * wt_g, W_full) + diag(penalty)
    varcov <- SpaNorm::invert_mat(info)[sel, sel, drop = FALSE]
  }
  # NULL df -> normal reference (fixed fit); scalar/vector df -> t (mixed fit)
  ptail <- function(t, lower.tail) {
    .ptByCol(t, if (is.null(W_full)) NULL else df, lower.tail)
  }
  se <- sqrt(psi_g * diag(varcov))
  t_stat <- alpha_g / se
  names(se) <- names(t_stat) <- colnames(Wsub)

  if (combine == "brown") {
    # correlation structure of the ResponseNiche coefficients -> Brown weights
    vc <- varcov[cov_niche, cov_niche, drop = FALSE]
    a <- diag(1 / diag(vc), nrow = nrow(vc))
    vc <- a %*% vc %*% t(a)
    vc <- stats::cov2cor(vc)
    vc <- poolr::mvnconv(vc, target = "m2lp", side = 1, cov2cor = FALSE)
    combineP <- function(xn) poolr::fisher(xn, side = 1, R = vc,
                                           adjust = "generalized")$p
    combinePsub <- function(xn, is_ct) {
      poolr::fisher(xn[is_ct], side = 1, R = vc[is_ct, is_ct, drop = FALSE],
                    adjust = "generalized")$p
    }
  } else {
    # Cauchy (ACAT) combination: correlation-agnostic, no poolr / mvnconv. The
    # tan transform blows up at p = 0/1, so clamp the one-sided p-values first.
    eps <- 1e-15
    clamp <- function(xn) pmin(pmax(xn, eps), 1 - eps)
    combineP <- function(xn) .cauchyCombine(matrix(clamp(xn), nrow = 1))
    combinePsub <- function(xn, is_ct) {
      .cauchyCombine(matrix(clamp(xn[is_ct]), nrow = 1))
    }
  }

  dirs <- lapply(c(p.pos = FALSE, p.neg = TRUE), function(lower.tail) {
    x <- ptail(t_stat, lower.tail)
    xn <- x[cov_niche]
    p_gene <- combineP(xn)
    p_ct <- vapply(uniq_index, function(ct) {
      combinePsub(xn, index_ct == ct)
    }, numeric(1))
    c(Gene = p_gene, p_ct)
  })

  list(t_stat = t_stat, se = se, p.pos = dirs$p.pos, p.neg = dirs$p.neg)
}

#' Select rows of a matrix or a torch tensor
#'
#' Row-subsetting helper for the gene sub-batching in
#' \code{.waldCauchyBlock()}: base R \code{[} does not carry over to torch
#' tensors (which need \code{torch_index_select()}), and the working weights
#' may be either depending on the backend.
#'
#' @param x a matrix or torch tensor.
#' @param ii integer row positions to keep.
#' @return \code{x} restricted to rows \code{ii}, same type as the input.
#' @noRd
.rowsOf <- function(x, ii) {
  if (SpaNorm::is_torch_tensor(x)) {
    idx <- torch::torch_tensor(as.integer(ii), dtype = torch::torch_long(),
                               device = x$device)
    return(torch::torch_index_select(x, 1, idx))
  }
  x[ii, , drop = FALSE]
}

#' Batched per-gene weighted Gram matrix, for a sub-batch of genes at once
#'
#' Computes the weighted information matrix \code{crossprod(W * wt_g, W)} =
#' \code{W' diag(wt_g) W} for every gene \code{g} in a sub-batch, returning a
#' \code{(batch, p, p)} array/tensor and optionally adding a shared ridge
#' penalty (mixed-effects) to every slice.
#'
#' On the accelerator this is a single batched matmul over the
#' \code{(batch, ncells, p)} weighted design; on CPU it is a per-gene
#' \code{crossprod()} loop (base R has no batched matmul, and the
#' construction is the cheap part relative to the inversion anyway).
#'
#' Peak memory here is \strong{linear in \code{p}}
#' (\code{batch x ncells x p}) and directly controllable via the sub-batch
#' size. An earlier implementation instead precomputed a Khatri-Rao /
#' face-splitting cross term of \code{W} (\code{ncells x p^2}, built once per
#' bandwidth and shared across all gene blocks), turning every gene's Gram
#' matrix into one row of a single large matmul. That is elegant, but it
#' scales \emph{quadratically} in \code{p} and -- being built once, outside
#' the gene loop -- could not be bounded by any choice of block size. On a
#' realistic mixed-effects design it is fatal: a 602-column random-intercept
#' design over 21,843 cells needs 63 GB, and a 4,906-column random-slope
#' design needs 4.2 TB, which made \code{combine = "cauchy"} unusable with
#' random effects on \emph{both} backends. Do not reintroduce that form
#' without bounding it by the gene sub-batch.
#'
#' @param W the design the Gram matrix is over (cells x p; a matrix or a
#'   torch tensor) -- \code{Wsub} for a fixed-effects fit, the full \code{W}
#'   for a mixed fit.
#' @param wt_block the sub-batch's working weights, \code{batch x ncells}.
#' @param penalty_diag if supplied, a length-\code{p} ridge penalty added to
#'   every slice's diagonal (mixed-effects only).
#' @param backend the resolved backend (unused on the base-R path; kept for
#'   signature symmetry with the other batched helpers).
#' @return a \code{(batch, p, p)} array (or torch tensor).
#' @noRd
.gramBatch <- function(W, wt_block, penalty_diag = NULL, backend = "cpu") {
  if (SpaNorm::is_torch_tensor(W)) {
    # (batch, ncells, p) weighted design. Weighting both factors by sqrt(wt)
    # (rather than one by wt) keeps the product exactly symmetric, which
    # linalg_cholesky() in invert_mat_batched() relies on; the working
    # weights 1/(1/mu + psi) are strictly positive, so the sqrt is safe.
    Wg <- torch::torch_sqrt(wt_block)$unsqueeze(3) * W$unsqueeze(1)
    info <- torch::torch_matmul(Wg$transpose(2, 3), Wg)
    if (!is.null(penalty_diag)) {
      pen <- torch::torch_tensor(as.numeric(penalty_diag),
                                 dtype = info$dtype, device = info$device)
      info <- info + torch::torch_diag(pen)$unsqueeze(1)
    }
    return(info)
  }

  p <- ncol(W)
  b <- nrow(wt_block)
  pen_mat <- if (is.null(penalty_diag)) NULL else diag(penalty_diag, nrow = p)
  info <- array(0, dim = c(b, p, p))
  for (g in seq_len(b)) {
    ig <- crossprod(W * wt_block[g, ], W)
    if (!is.null(pen_mat)) ig <- ig + pen_mat
    info[g, , ] <- ig
  }
  info
}

#' Subset a (block, p, p) array/tensor to (block, k, k) via a column index
#'
#' On a base R array, \code{varcov[, sel, sel]} already takes the
#' \code{sel x sel} cross-block (R's array indexing is orthogonal per
#' dimension). On a torch tensor this is NOT true -- tensor \code{[}-indexing
#' with two non-scalar index vectors pairs them elementwise (NumPy/PyTorch
#' "advanced indexing" semantics), so the two dimensions must be subset one
#' at a time via \code{torch_index_select()} instead.
#'
#' @param varcov a \code{(block, p, p)} array or torch tensor.
#' @param sel integer positions to keep along both of the last two dimensions.
#' @return a \code{(block, k, k)} array or torch tensor, \code{k = length(sel)}.
#' @noRd
.subsetBatch <- function(varcov, sel) {
  if (SpaNorm::is_torch_tensor(varcov)) {
    sel_t <- torch::torch_tensor(as.integer(sel), dtype = torch::torch_long(),
                                 device = varcov$device)
    varcov <- torch::torch_index_select(varcov, 2, sel_t)
    return(torch::torch_index_select(varcov, 3, sel_t))
  }
  varcov[, sel, sel, drop = FALSE]
}

#' Extract the per-slice diagonal of a (block, p, p) array/tensor
#'
#' @param varcov a \code{(block, p, p)} array or torch tensor.
#' @return a \code{block x p} matrix (or torch tensor) of diagonals.
#' @noRd
.batchDiag <- function(varcov) {
  if (SpaNorm::is_torch_tensor(varcov)) {
    return(torch::torch_diagonal(varcov, dim1 = 2, dim2 = 3))
  }
  b <- dim(varcov)[1]
  p <- dim(varcov)[2]
  idx <- cbind(rep(seq_len(b), p), rep(seq_len(p), each = b),
              rep(seq_len(p), each = b))
  matrix(varcov[idx], nrow = b, ncol = p)
}

#' Batched Wald covariance + Cauchy combination for a whole gene block
#'
#' Batched replacement for the per-gene \code{.waldBrownGene()} loop when
#' \code{combine == "cauchy"}: per-gene weighted Gram matrices are built for a
#' sub-batch of genes at a time (\code{.gramBatch()}), inverted in one batched
#' Cholesky call (\code{SpaNorm::invert_mat_batched()}), and the within-gene
#' Cauchy combination is vectorized across the block by looping only over
#' \code{uniq_index} (a handful of cell types) rather than over genes --
#' \code{.cauchyCombine()} already accepts a genes x k matrix, it is just
#' never fed more than one row at a time by the per-gene path. Only valid for
#' \code{combine == "cauchy"}: Brown's method needs a gene-specific
#' correlation matrix and cannot batch this way, so it keeps using
#' \code{.waldBrownGene()} per gene.
#'
#' @param alpha_block genes(block) x k coefficients for the Response/
#'   ResponseNiche columns.
#' @param Wsub the Response/ResponseNiche sub-design (cells x k); only used
#'   for its column names.
#' @param wtb the block's working weights, genes(block) x cells.
#' @param scale_block the variance scale per gene (length genes(block)): NB
#'   dispersion for a fixed-effects fit, working (Pearson) dispersion for a
#'   mixed fit.
#' @param cov_niche,index_ct,uniq_index as in \code{.waldBrownGene()}.
#' @param Wgram the design the Gram matrix is taken over: \code{Wsub} for a
#'   fixed-effects fit, the full \code{W} for a mixed fit (already on the
#'   device when the GPU backend is active).
#' @param W_full,penalty,sel,df as in \code{.waldBrownGene()}; \code{W_full}
#'   NULL for the fixed-effects fit.
#' @param backend the resolved backend, forwarded to the GPU-aware helpers.
#' @param cov.batch genes per covariance sub-batch (NULL = all of them at
#'   once); see \code{.covBatchSize()}.
#' @return a list with t_stat, se (genes(block) x k) and p.pos, p.neg
#'   (genes(block) x (1 + n_index)) -- the same shape as stacking
#'   \code{.waldBrownGene()}'s per-gene output over the block.
#' @noRd
.waldCauchyBlock <- function(alpha_block, Wsub, wtb, scale_block, cov_niche,
                             index_ct, uniq_index, Wgram, W_full = NULL,
                             penalty = NULL, sel = NULL, df = NULL,
                             backend = "cpu", cov.batch = NULL) {
  b <- nrow(alpha_block)

  # The only per-gene quantity the Wald statistics need out of the (p x p)
  # covariance is the diagonal of its Response/ResponseNiche sub-block -- a
  # length-k vector, with k = ncol(Wsub) small regardless of how wide the
  # full design is. So the (batch, p, p) Gram/inverse stack is built, used
  # and discarded one sub-batch at a time, and only the (b, k) diagonals are
  # accumulated. This is what bounds peak memory for wide mixed-effects
  # designs independently of the caller's gene block.size.
  diagB <- matrix(0, nrow = b, ncol = ncol(Wsub))
  on_cpu <- !SpaNorm::is_torch_tensor(wtb)
  if (on_cpu && !is.null(cov.batch) && as.integer(cov.batch) == 1L) {
    # Very wide designs (large p) can only afford one gene per sub-batch, and
    # at that size the batched machinery is pure overhead: the (1, p, p)
    # array, invert_mat_batched()'s vapply/aperm round-trip and the
    # .subsetBatch() copy each duplicate a p x p matrix (192 MB at p = 4906)
    # to do the work of a single invert_mat() call. Go direct instead --
    # measurably faster, and identical by construction, which the
    # cov.batch-invariance test pins down.
    pen_mat <- if (is.null(penalty)) NULL else diag(penalty, nrow = ncol(Wgram))
    for (g in seq_len(b)) {
      info <- crossprod(Wgram * wtb[g, ], Wgram)
      if (!is.null(pen_mat)) info <- info + pen_mat
      vc_diag <- diag(SpaNorm::invert_mat(info))
      diagB[g, ] <- if (is.null(W_full)) vc_diag else vc_diag[sel]
    }
  } else {
    for (ii in .chunkGenes(b, cov.batch)) {
      info <- .gramBatch(Wgram, .rowsOf(wtb, ii), penalty_diag = penalty,
                         backend = backend)
      vc <- SpaNorm::invert_mat_batched(info)
      if (!is.null(W_full)) {
        vc <- .subsetBatch(vc, sel)
      }
      diagB[ii, ] <- SpaNorm::toRMatrix(.batchDiag(vc))
    }
  }
  scale_block <- as.numeric(SpaNorm::toRMatrix(scale_block))
  se <- sqrt(scale_block * diagB)
  colnames(se) <- colnames(Wsub)
  t_stat <- SpaNorm::toRMatrix(alpha_block) / se
  colnames(t_stat) <- colnames(Wsub)

  ptail <- function(t, lower.tail) {
    .ptByCol(t, if (is.null(W_full)) NULL else df, lower.tail)
  }

  eps <- 1e-15
  dirs <- lapply(c(p.pos = FALSE, p.neg = TRUE), function(lower.tail) {
    p_all <- ptail(t_stat, lower.tail)
    p_niche <- pmin(pmax(p_all[, cov_niche, drop = FALSE], eps), 1 - eps)
    p_gene <- .cauchyCombine(p_niche)
    p_ct <- vapply(uniq_index, function(ct) {
      .cauchyCombine(p_niche[, index_ct == ct, drop = FALSE])
    }, numeric(b))
    # vapply drops to a plain (unnamed-dim) vector when b == 1 (FUN.VALUE has
    # length 1), instead of the (b, n_index) matrix it returns for b > 1 --
    # force the matrix shape unconditionally (a no-op reshape when b > 1,
    # since vapply's own column-major layout already matches) and restore the
    # column names vapply's auto-naming would otherwise have supplied.
    p_ct <- matrix(p_ct, nrow = b, ncol = length(uniq_index),
                   dimnames = list(NULL, uniq_index))
    cbind(Gene = p_gene, p_ct)
  })

  list(t_stat = t_stat, se = se, p.pos = dirs$p.pos, p.neg = dirs$p.neg)
}

#' Number of genes per inference block under a GPU memory budget
#'
#' Picks \code{block.size} for the GPU inference path so the per-block
#' genes x cells intermediates plus the batched covariance intermediates stay
#' within the accelerator's memory budget. Returns \code{NULL} (a single
#' block, the pre-existing behaviour) whenever the CPU path is in use -- so
#' the CPU path is never affected by GPU blocking, exactly mirroring
#' \code{SpaNorm::geneBlockCount()}'s own invariant.
#'
#' This sizes the \strong{NB math} only -- \code{Yb}, \code{mub}, \code{wtb},
#' the NB log-pmf internals and (mixed) the dispersion intermediates, all
#' \code{block x ncells} -- via \code{SPIDE_TENSOR_MULT_NB}. The covariance
#' side is \emph{not} sized here: \code{.waldCauchyBlock()} sub-batches the
#' \code{(batch, p, p)} Gram/inverse stack internally (see
#' \code{.covBatchSize()}), so it is bounded whatever block size is chosen,
#' and the two concerns stay independent. The multiplier is a rough,
#' deliberately generous estimate that needs empirical validation on real
#' hardware, exactly as noted for SpaNorm's own
#' \code{GPU_BLOCK_TENSOR_MULTIPLIER}.
#'
#' @param ng number of genes.
#' @param ncells number of cells.
#' @param p_eff the covariance dimension (\code{ncol(Wsub)} fixed effects,
#'   \code{ncol(W_full)} mixed); retained for signature stability and future
#'   use -- the covariance term is now bounded by \code{.covBatchSize()}.
#' @param backend the resolved backend.
#' @param gpu.mem.budget \code{NULL} (auto-detect via
#'   \code{SpaNorm::getGPUMemoryBudget()}) or a budget in bytes.
#' @return a single integer block size, or \code{NULL} for a single block.
#' @noRd
SPIDE_TENSOR_MULT_NB <- 12
SPIDE_TENSOR_MULT_COV <- 6
# fraction of the memory budget each of the two independently-bounded stages
# (NB math per gene block, covariance per sub-batch) may claim
SPIDE_BUDGET_FRACTION <- 0.5
# default cap on the covariance stack when there is no GPU budget to consult
# (the CPU path); override with options(spiDE.cov.mem.budget = <bytes>)
SPIDE_COV_MEM_BUDGET_CPU <- 2e9

.inferenceBlockSize <- function(ng, ncells, p_eff, backend,
                                gpu.mem.budget = NULL) {
  if (!(backend %in% c("gpu", "auto") && SpaNorm::checkGPU())) {
    return(NULL)
  }
  budget <- SpaNorm::getGPUMemoryBudget(gpu.mem.budget)
  if (!is.finite(budget)) {
    return(NULL)
  }
  bytes <- SpaNorm::gpuDtypeBytes()
  per_gene <- bytes * as.numeric(ncells) * SPIDE_TENSOR_MULT_NB
  b <- max(1L, as.integer(floor(budget * SPIDE_BUDGET_FRACTION / per_gene)))
  if (b >= ng) NULL else b
}

#' Number of genes per covariance sub-batch
#'
#' Bounds the \code{(batch, p, p)} Gram/inverse stack (and, on the GPU, the
#' \code{(batch, ncells, p)} weighted design feeding it) inside
#' \code{.waldCauchyBlock()}. This is the knob that makes wide mixed-effects
#' designs tractable: peak covariance memory is linear in \code{p} and scales
#' with this batch, so a design too wide to process all at once is handled by
#' shrinking the sub-batch rather than failing.
#'
#' Applies on \strong{both} backends. The CPU path needs it just as much as
#' the GPU path -- a single block of 13,348 genes at \code{p = 4906} would
#' otherwise try to allocate a \code{(13348, 4906, 4906)} array -- and unlike
#' \code{.inferenceBlockSize()} it therefore does not return NULL for CPU.
#'
#' @param ncells number of cells.
#' @param p the covariance dimension (\code{ncol} of the Gram design).
#' @param backend the resolved backend.
#' @param gpu.mem.budget \code{NULL} (auto-detect) or a budget in bytes; only
#'   consulted on the GPU path.
#' @return a single integer, genes per covariance sub-batch (at least 1).
#' @noRd
.covBatchSize <- function(ncells, p, backend, gpu.mem.budget = NULL) {
  gpu_active <- backend %in% c("gpu", "auto") && SpaNorm::checkGPU()
  budget <- if (gpu_active) {
    SpaNorm::getGPUMemoryBudget(gpu.mem.budget)
  } else {
    getOption("spiDE.cov.mem.budget", SPIDE_COV_MEM_BUDGET_CPU)
  }
  if (!is.finite(budget)) {
    budget <- SPIDE_COV_MEM_BUDGET_CPU
  }
  bytes <- if (gpu_active) SpaNorm::gpuDtypeBytes() else 8
  ncells <- as.numeric(ncells)
  p <- as.numeric(p)

  # GPU: the (batch, ncells, p) weighted design plus its transpose view, then
  # the (batch, p, p) Gram/Cholesky/inverse stack. CPU: the (batch, p, p)
  # stack only -- construction there is one p x p crossprod at a time.
  per_gene <- if (gpu_active) {
    bytes * (ncells * p * 2 + p^2 * SPIDE_TENSOR_MULT_COV)
  } else {
    bytes * p^2 * SPIDE_TENSOR_MULT_COV
  }
  max(1L, as.integer(floor(budget * SPIDE_BUDGET_FRACTION / per_gene)))
}

#' Compute Wald inference + within-gene combination for a SpiDEFit, block-wise
#'
#' Fills the \code{t_stat}, \code{se}, \code{p.combined.pos}, \code{p.combined.neg}
#' and (recomputed) \code{loglik} slots. Genes are processed in blocks so the
#' genes x cells fitted-mean / weight matrices are never fully materialised, and
#' blocks are dispatched with BiocParallel.
#'
#' The block-wide NB math (fitted means, working weights, log-likelihood) and
#' the within-gene Wald covariance run on the accelerator when \code{backend}
#' resolves to one: per-gene weighted Gram matrices are built in a batched
#' matmul and inverted in one batched Cholesky call
#' (\code{.waldCauchyBlock()}), rather than one \code{crossprod()}/inverse per
#' gene. This batched path is used for \code{combine == "cauchy"} only;
#' \code{combine == "brown"} needs a gene-specific correlation matrix and so
#' keeps the per-gene \code{.waldBrownGene()} loop (its NB math is still
#' GPU-accelerated).
#'
#' Two memory concerns are bounded independently: \code{block.size} (via
#' \code{.inferenceBlockSize()}) sizes the \code{block x ncells} NB math,
#' while \code{.waldCauchyBlock()} sub-batches the \code{(batch, p, p)}
#' covariance stack internally (via \code{.covBatchSize()}). Keeping these
#' separate is what lets a very wide mixed-effects design run at any
#' \code{block.size}.
#'
#' @param fit a SpiDEFit (with the NB fit populated).
#' @param Y the counts matrix (dense, sparse, or DelayedArray).
#' @param block.size a numeric, genes per block (NULL = a single block on the
#'   CPU path, or an auto-selected memory-bounded size on the GPU path).
#' @param combine one of "cauchy" (default) or "brown".
#' @param backend one of "cpu" (default), "auto" or "gpu".
#' @param gpu.mem.budget NULL (auto-detect) or a GPU memory budget in bytes.
#' @param BPPARAM a BiocParallelParam. Forced to \code{SerialParam()} (with a
#'   warning) when the GPU backend is active and more than one worker is
#'   requested, to avoid multiple processes contending for one device.
#' @return the input \code{fit} with inference slots populated.
#' @importFrom stats dnbinom
#' @importFrom BiocParallel bplapply SerialParam bpnworkers
#' @noRd
.blockedInference <- function(fit, Y, block.size = NULL,
                              combine = c("cauchy", "brown"),
                              backend = c("cpu", "auto", "gpu"),
                              gpu.mem.budget = NULL,
                              BPPARAM = BiocParallel::SerialParam()) {
  combine <- match.arg(combine)
  backend <- match.arg(backend)
  W_full <- fit@W
  covtype <- as.character(fit@covtype)
  cols_gene <- grepl("Response", covtype) # Response + ResponseNiche columns
  Wsub <- W_full[, cols_gene, drop = FALSE]
  cov_niche <- covtype[cols_gene] == "ResponseNiche"
  coefmap_sub <- fit@coefmap[cols_gene, , drop = FALSE]
  index_ct <- coefmap_sub$index[cov_niche]
  uniq_index <- unique(index_ct)

  # mixed-effects fit: use the full penalised information, the working (Pearson)
  # dispersion, and a between-patient t reference
  has_re <- !is.null(fit@penalty)
  sel <- if (has_re) which(cols_gene) else NULL
  penalty <- if (has_re) fit@penalty else NULL
  re_df <- if (has_re) fit@df else NULL
  # residual df for the working-dispersion estimate (fixed columns only; the
  # penalised random columns contribute little effective df)
  disp_df <- if (has_re) max(nrow(W_full) - sum(!grepl("Random", covtype)), 1)

  alpha_full <- fit@alpha
  psi <- fit@psi
  ng <- nrow(alpha_full)

  # resolve GPU state once, in the parent process, before any dispatch. The
  # backend %in% c("gpu","auto") test short-circuits before checkGPU() probes
  # torch, so a forked worker never touches the GPU runtime.
  gpu_active <- backend %in% c("gpu", "auto") && SpaNorm::checkGPU()
  if (gpu_active && BiocParallel::bpnworkers(BPPARAM) > 1) {
    warning("GPU backend active (", SpaNorm::getBackendDevice(), "); forcing ",
            "BPPARAM = BiocParallel::SerialParam() to avoid multiple ",
            "processes contending for one GPU device. GPU inference is ",
            "already batched across genes and does not benefit from multiple ",
            "CPU workers.", call. = FALSE)
    BPPARAM <- BiocParallel::SerialParam()
  }

  # auto block size on the GPU path (bounded by the device memory budget);
  # NULL on the CPU path preserves the pre-existing single-block behaviour.
  if (is.null(block.size) && gpu_active) {
    p_eff <- if (has_re) ncol(W_full) else ncol(Wsub)
    block.size <- .inferenceBlockSize(ng, nrow(W_full), p_eff, backend,
                                      gpu.mem.budget)
  }
  blocks <- .chunkGenes(ng, block.size)

  # the design is shared across all genes, so it is built (and pushed to the
  # device) once, not per block. The Gram design is Wsub for a fixed-effects
  # fit, the full W for a mixed fit; it is only needed by the batched
  # "cauchy" path, since combine == "brown" builds its own covariance
  # per gene in .waldBrownGene().
  W_full_dev <- if (gpu_active) SpaNorm::toGPUMatrix(W_full, backend = backend)
                else W_full
  if (combine == "cauchy") {
    Wgram <- if (has_re) {
      W_full_dev
    } else if (gpu_active) {
      SpaNorm::toGPUMatrix(Wsub, backend = backend)
    } else {
      Wsub
    }
    # genes per covariance sub-batch -- bounds the (batch, p, p) stack
    # independently of block.size, on both backends
    cov_batch <- .covBatchSize(nrow(W_full),
                               if (has_re) ncol(W_full) else ncol(Wsub),
                               backend, gpu.mem.budget)
  }

  block_res <- BiocParallel::bplapply(blocks, function(gi) {
    Yb <- as.matrix(Y[gi, , drop = FALSE])
    alpha_block <- alpha_full[gi, , drop = FALSE]
    psib <- psi[gi]
    zero_gmean <- rep(0, length(gi))

    if (gpu_active) {
      Yb_dev <- SpaNorm::toGPUMatrix(Yb, backend = backend)
      mub <- SpaNorm::calculateMu(zero_gmean, alpha_block, W_full_dev,
                                  backend = backend) # tensor, block x ncells
      inv_mu <- torch::torch_reciprocal(mub)
      wtb <- torch::torch_reciprocal(
        SpaNorm::add_vec_mat_gpu(psib, inv_mu, backend = backend))
      loglikb <- as.numeric(SpaNorm::toRMatrix(SpaNorm::rowSums_gpu(
        SpaNorm::dnbinom_gpu(Yb_dev, mu = mub, size = 1 / psib, log = TRUE))))
      if (has_re) {
        num <- (Yb_dev - mub)^2
        denom <- mub + SpaNorm::mult_vec_mat_gpu(psib, mub * mub,
                                                 backend = backend)
        dispb <- as.numeric(SpaNorm::toRMatrix(
          SpaNorm::rowSums_gpu(num / denom))) / disp_df
      }
    } else {
      mub <- SpaNorm::calculateMu(zero_gmean, alpha_block, W_full)
      wtb <- 1 / (1 / mub + psib) # nblock x ncells
      loglikb <- rowSums(dnbinom(Yb, mu = mub, size = 1 / psib, log = TRUE))
      if (has_re) {
        dispb <- rowSums((Yb - mub)^2 / (mub + psib * mub^2)) / disp_df
      }
    }
    scale_block <- if (has_re) dispb else psib

    if (combine == "cauchy") {
      res <- .waldCauchyBlock(alpha_block[, cols_gene, drop = FALSE], Wsub,
                              wtb, scale_block, cov_niche, index_ct,
                              uniq_index, Wgram,
                              W_full = if (has_re) W_full else NULL,
                              penalty = penalty, sel = sel, df = re_df,
                              backend = backend, cov.batch = cov_batch)
      res$loglik <- loglikb
      return(res)
    }

    # Brown's method: per-gene, needs plain-R weights and design
    wtb_r <- if (gpu_active) SpaNorm::toRMatrix(wtb) else wtb
    per_gene <- lapply(seq_along(gi), function(i) {
      .waldBrownGene(alpha_block[i, cols_gene], Wsub, wtb_r[i, ],
                     scale_block[i], cov_niche, index_ct, uniq_index,
                     combine = "brown", W_full = if (has_re) W_full else NULL,
                     penalty = penalty, sel = sel, df = re_df)
    })
    list(
      t_stat = do.call(rbind, lapply(per_gene, `[[`, "t_stat")),
      se = do.call(rbind, lapply(per_gene, `[[`, "se")),
      p.pos = do.call(rbind, lapply(per_gene, `[[`, "p.pos")),
      p.neg = do.call(rbind, lapply(per_gene, `[[`, "p.neg")),
      loglik = loglikb
    )
  }, BPPARAM = BPPARAM)

  bind <- function(field) do.call(rbind, lapply(block_res, `[[`, field))
  t_stat <- bind("t_stat")
  se <- bind("se")
  p.pos <- bind("p.pos")
  p.neg <- bind("p.neg")
  loglik <- unlist(lapply(block_res, `[[`, "loglik"), use.names = FALSE)

  gnames <- rownames(alpha_full)
  rownames(t_stat) <- rownames(se) <- gnames
  rownames(p.pos) <- rownames(p.neg) <- gnames
  names(loglik) <- gnames

  fit@t_stat <- t_stat
  fit@se <- se
  fit@p.combined.pos <- p.pos
  fit@p.combined.neg <- p.neg
  fit@loglik <- loglik
  fit
}
