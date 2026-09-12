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

#' Should multi-worker dispatch run each worker's BLAS single-threaded?
#' @noRd
.singleBLAS <- function(BPPARAM) BiocParallel::bpnworkers(BPPARAM) > 1L

#' Run this worker's BLAS and OpenMP single-threaded
#'
#' Called inside a \code{bplapply()} worker of a multi-worker BPPARAM, never
#' in the parent.
#' @return the BLAS thread count before the call, so a persistent (snow)
#'   worker can be restored after the block.
#' @noRd
.workerBLAS <- function() {
  prev <- RhpcBLASctl::blas_get_num_procs()
  if (isTRUE(prev > 1L)) RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  invisible(prev)
}

#' bplapply() with single-threaded BLAS inside multi-worker dispatch
#'
#' A forked worker inherits the parent's OpenBLAS thread count, so N workers
#' on N cores run N x threads BLAS threads. Measured at the cohort's design
#' shape (77,454 cells, 345 dense + 660 nested columns; 2026-09-10): 4
#' workers x 4 threads on 4 cores take 2.3-2.6 s per Newton step against
#' 0.24-0.28 s at one thread each, while one worker gains only 1.5x from 4
#' threads -- the per-gene gram is memory-bound. The 0.99.19 cohort runs (4
#' workers x 8 threads on 8 cores) spent 245-395 min on the cold polish pass
#' this way. Both blocked stages -- the polish and the inference, whose
#' per-gene gram is the same BLAS work -- dispatch through this wrapper; the
#' niche construction and the gene-set stage are not gram-bound and do not.
#' A persistent worker gets its thread count back when the element is done;
#' a forked one dies with the call. Serial dispatch is plain
#' \code{bplapply()}.
#' @noRd
.bplapplySingleBLAS <- function(X, FUN, ..., BPPARAM = BiocParallel::SerialParam()) {
  if (!.singleBLAS(BPPARAM)) {
    return(BiocParallel::bplapply(X, FUN, ..., BPPARAM = BPPARAM))
  }
  BiocParallel::bplapply(X, function(x, ...) {
    prev <- .workerBLAS()
    on.exit(if (isTRUE(prev > 1L)) RhpcBLASctl::blas_set_num_threads(prev), add = TRUE)
    FUN(x, ...)
  }, ..., BPPARAM = BPPARAM)
}

#' Apply a t tail to a statistic matrix, with per-column df
#'
#' \code{df} may be a scalar (one df for all columns),
#' a length-\code{ncol(tmat)} vector (a df per tested column, broadcast down the
#' rows), or a matrix matching \code{tmat}. \code{tmat} may be a genes x k matrix
#' or a length-k vector (a single gene). Centralises the reference-distribution
#' logic so both the inference and FDR stages consume \code{@df} identically.
#'
#' There is no separate normal-reference branch: a fixed-effects fit carries
#' \code{df = Inf}, and \code{pt(x, Inf)} is exactly \code{pnorm(x)}. NULL is
#' still tolerated for objects serialised before \code{@df} was always
#' populated.
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

#' The reference df for a subset of the tested columns
#'
#' \code{@df} is either a single value shared by every tested column (a
#' fixed-effects \code{Inf}, or the scalar \code{S - 2} of
#' \code{df.method = "between"}) or one value per tested column
#' (\code{df.method = "satterthwaite"}). Consumers select a subset of the
#' tested columns and need the matching df; this is the one place that
#' distinction is made, rather than the same conditional at four call sites.
#'
#' Length 0 (a NULL \code{@df} on an object serialised before the slot was
#' always populated) passes through unchanged.
#'
#' @param f a SpiDEFit.
#' @param sel a logical or integer index into the tested columns.
#' @return the df for those columns.
#' @noRd
.dfFor <- function(f, sel) {
  if (length(f@df) <= 1L) f@df else f@df[sel]
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
.blockLoglik <- function(Y, alpha, W, psi, block.size = 2000L, winsor = 4) {
  ng <- nrow(alpha)
  blocks <- .chunkGenes(ng, block.size)
  ll <- numeric(ng)
  for (gi in blocks) {
    Yb <- as.matrix(Y[gi, , drop = FALSE])
    # ONE mean function per fit: fitNB's coefficients are read back through
    # calculateMu()'s 4-MAD clamp, exactly as inference reads them; a polished
    # fit is defined on the unclamped, floored mean and must be read back the
    # same way (winsor = Inf), or its log-likelihood is evaluated at a mean it
    # was never converged on -- the mismatch the inference review removed.
    mub <- SpaNorm::calculateMu(rep(0, length(gi)),
                                alpha[gi, , drop = FALSE], W, winsor = winsor)
    mub <- pmax(mub, .MU_FLOOR)
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
                           sel = NULL, df = NULL, w_rc = NULL) {
  combine <- match.arg(combine)
  # The information matrix is built from THIS gene's working weights, so it can
  # be singular for one gene and fine for the next. Unguarded, that aborts the
  # whole run: measured 2026-08-13, one gene of 13,348 (rcond 2.7e-20) killed an
  # 80-minute fit, while the same design on a 200-gene subset completed. A
  # singular gene must drop out (NA), not take the run with it.
  varcov <- tryCatch({
    if (is.null(W_full)) {
      # fixed-effects fit: covariance from the tested sub-design, normal ref
      SpaNorm::invert_mat(crossprod(Wsub * wt_g, Wsub))
    } else {
      # mixed-effects fit: full penalised information (X'WX + Lambda), then the
      # fixed-effect (Response/ResponseNiche) block -> larger between-sample SEs
      info <- crossprod(W_full * wt_g, W_full) + diag(penalty)
      SpaNorm::invert_mat(info)[sel, sel, drop = FALSE]
    }
  }, error = function(e) NULL)
  if (is.null(varcov)) {
    # shape MUST match the normal return exactly (t_stat, se, p.pos, p.neg,
    # se_pat), with p vectors named c("Gene", <index types>) as the live path
    # builds them -- a mismatched NA branch would fail only when a singular
    # gene appears, i.e. precisely when this guard is doing its job.
    nm <- colnames(Wsub)
    pnm <- c("Gene", uniq_index)
    return(list(
      t_stat = stats::setNames(rep(NA_real_, ncol(Wsub)), nm),
      se = stats::setNames(rep(NA_real_, ncol(Wsub)), nm),
      p.pos = stats::setNames(rep(NA_real_, length(pnm)), pnm),
      p.neg = stats::setNames(rep(NA_real_, length(pnm)), pnm),
      se_pat = NA_real_))
  }
  # NULL df -> normal reference (fixed fit); scalar/vector df -> t (mixed fit)
  ptail <- function(t, lower.tail) {
    .ptByCol(t, df, lower.tail)
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
    # .clampP() is the one place the ACAT floor is defined; it used to be
    # written out here and again in .waldCauchyBlock(), so the within-gene and
    # cross-bandwidth combinations could silently drift onto different floors.
    combineP <- function(xn) .cauchyCombine(matrix(.clampP(xn), nrow = 1))
    combinePsub <- function(xn, is_ct) {
      .cauchyCombine(matrix(.clampP(xn[is_ct]), nrow = 1))
    }
  }

  # Combination side. Brown/Fisher uses -2log(p), which is bounded at 0, so a
  # p-value near 1 contributes nothing and one-sided input is safe. ACAT uses
  # tan((0.5 - p)*pi), which diverges to -Inf as p -> 1, so a gene UP in one
  # niche and DOWN in another cancels exactly: the opposing column's one-sided
  # p is ~1 and its -3e14 term annihilates the signal. Verified: with a
  # bidirectional effect the one-sided ACAT gene p sits at ~0.5 at EVERY effect
  # size, while Brown's falls monotonically to 2e-11. So ACAT combines
  # TWO-SIDED p-values; Brown keeps one-sided.
  if (combine == "cauchy") {
    p_two <- 2 * pmin(ptail(t_stat, FALSE), ptail(t_stat, TRUE))
    p_two <- pmin(p_two, 1)
    xn <- p_two[cov_niche]
    comb <- c(Gene = combineP(xn),
              vapply(uniq_index, function(ct) combinePsub(xn, index_ct == ct),
                     numeric(1)))
    # the same two-sided combination is returned on both sides; the cascade is
    # told via two.sided that it must not apply the fdr/2 direction split, and
    # direction is taken from the niche-level coefficients instead.
    dirs <- list(p.pos = comb, p.neg = comb)
  } else {
    dirs <- lapply(c(p.pos = FALSE, p.neg = TRUE), function(lower.tail) {
      x <- ptail(t_stat, lower.tail)
      xn <- x[cov_niche]
      p_gene <- combineP(xn)
      p_ct <- vapply(uniq_index, function(ct) {
        combinePsub(xn, index_ct == ct)
      }, numeric(1))
      c(Gene = p_gene, p_ct)
    })
  }

  # patient-level contrast SE: sqrt(psi * w'Vw). beta = w'alpha is recovered
  # from the coefficients outside, so only the variance needs V.
  se_pat <- if (is.null(w_rc)) NA_real_ else {
    sqrt(psi_g * as.numeric(crossprod(w_rc, varcov %*% w_rc)))
  }
  list(t_stat = t_stat, se = se, p.pos = dirs$p.pos, p.neg = dirs$p.neg,
       se_pat = se_pat)
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
                             index_ct, uniq_index, Wgram, w_rc = NULL,
                             W_full = NULL,
                             penalty = NULL, sel = NULL, df = NULL,
                             backend = "cpu", cov.batch = NULL,
                             absorb = NULL) {
  b <- nrow(alpha_block)

  # The only per-gene quantity the Wald statistics need out of the (p x p)
  # covariance is the diagonal of its Response/ResponseNiche sub-block -- a
  # length-k vector, with k = ncol(Wsub) small regardless of how wide the
  # full design is. So the (batch, p, p) Gram/inverse stack is built, used
  # and discarded one sub-batch at a time, and only the (b, k) diagonals are
  # accumulated. This is what bounds peak memory for wide mixed-effects
  # designs independently of the caller's gene block.size.
  diagB <- matrix(0, nrow = b, ncol = ncol(Wsub))
  # patient-level contrast variance w'Vw, accumulated while V exists
  quadB <- rep(0, b)
  on_cpu <- !SpaNorm::is_torch_tensor(wtb)
  if (on_cpu && !is.null(absorb)) {
    # The nested (sample x cell type) block is a 0/1 partition of the cells, so
    # its blocks of the information matrix are a diagonal and a rowsum away and
    # can be absorbed by a Schur complement. The tested columns are never in
    # that block, and S^-1 IS the covariance restricted to the dense columns
    # (test-polish.R pins this against the dense inverse, including the
    # diagonal and the patient contrast). So the per-gene gram runs over the
    # dense columns alone instead of building and inverting the full design's
    # gram and then discarding the block.
    #
    # Measured at the real cohort shape (77,454 cells, 345 dense + 660 nested
    # columns): 2.16 s -> 0.35 s per gene, a 6.2x saving, or 8.0 h -> 1.3 h for
    # a 13,348-gene transcriptome, agreeing to 7e-21. The saving is NOT uniform
    # -- at a smaller shape (20,000 cells, 400 groups) the two are within 20%
    # and the absorption can even lose, because the Schur correction and the
    # group sums are then a bigger share of a smaller gram. It is kept
    # unconditional because at those sizes both paths cost well under a second
    # per gene, so only the large case is worth optimising for.
    for (g in seq_len(b)) {
      vcg <- absorb$solver$xcov(wtb[g, ])
      if (is.null(vcg)) {
        diagB[g, ] <- NA_real_
        if (!is.null(w_rc)) quadB[g] <- NA_real_
        next
      }
      vcg <- vcg[absorb$sel_x, absorb$sel_x, drop = FALSE]
      diagB[g, ] <- diag(vcg)
      if (!is.null(w_rc)) quadB[g] <- as.numeric(crossprod(w_rc, vcg %*% w_rc))
    }
  } else if (on_cpu && !is.null(cov.batch) && as.integer(cov.batch) == 1L) {
    # Very wide designs (large p) can only afford one gene per sub-batch, and
    # at that size the batched machinery is pure overhead: the (1, p, p)
    # array, invert_mat_batched()'s vapply/aperm round-trip and the
    # .subsetBatch() copy each duplicate a p x p matrix (192 MB at p = 4906)
    # to do the work of a single invert_mat() call. Go direct instead --
    # measurably faster, and identical by construction, which the
    # cov.batch-invariance test pins down.
    pen_mat <- if (is.null(penalty)) NULL else diag(penalty, nrow = ncol(Wgram))
    # NA (not error) for a gene whose weighted Gram is singular -- see the note
    # in .waldBrownGene(). NA propagates to se/t_stat and then to the p-values,
    # and p.adjust() already ignores NA, so such a gene is simply not tested.
    for (g in seq_len(b)) {
      info <- crossprod(Wgram * wtb[g, ], Wgram)
      if (!is.null(pen_mat)) info <- info + pen_mat
      vcg <- tryCatch(SpaNorm::invert_mat(info), error = function(e) NULL)
      if (is.null(vcg)) {
        diagB[g, ] <- NA_real_
        if (!is.null(w_rc)) quadB[g] <- NA_real_
        next
      }
      if (!is.null(W_full)) vcg <- vcg[sel, sel, drop = FALSE]
      diagB[g, ] <- diag(vcg)
      if (!is.null(w_rc)) quadB[g] <- as.numeric(crossprod(w_rc, vcg %*% w_rc))
    }
  } else {
    for (ii in .chunkGenes(b, cov.batch)) {
      info <- .gramBatch(Wgram, .rowsOf(wtb, ii), penalty_diag = penalty,
                         backend = backend)
      # Same singular-gene exposure as the direct path above, but batched: one
      # bad gene fails the whole Cholesky for its sub-batch. A per-gene fallback
      # here would have to reproduce the torch/base split of .gramBatch(), which
      # cannot be exercised without a GPU, so this reports a diagnosis with the
      # concrete lever instead of guessing at a recovery.
      vc <- tryCatch(SpaNorm::invert_mat_batched(info), error = function(e)
        stop("batched covariance inversion failed for a gene sub-batch: ",
             conditionMessage(e), "\n  A singular per-gene information matrix ",
             "is the usual cause (thin or collinear design: few cells relative ",
             "to design columns).\n  Levers: lower `cov.batch` to isolate the ",
             "offending genes, or set backend = 'cpu', which takes the ",
             "per-gene path that drops singular genes as NA instead of failing.",
             call. = FALSE))
      if (!is.null(W_full)) {
        vc <- .subsetBatch(vc, sel)
      }
      diagB[ii, ] <- SpaNorm::toRMatrix(.batchDiag(vc))
      if (!is.null(w_rc)) quadB[ii] <- .batchQuad(vc, w_rc)
    }
  }
  scale_block <- as.numeric(SpaNorm::toRMatrix(scale_block))
  se <- sqrt(scale_block * diagB)
  colnames(se) <- colnames(Wsub)
  t_stat <- SpaNorm::toRMatrix(alpha_block) / se
  colnames(t_stat) <- colnames(Wsub)

  ptail <- function(t, lower.tail) {
    .ptByCol(t, df, lower.tail)
  }

  # This is the ACAT path, so combine TWO-SIDED p-values (see .waldBrownGene for
  # why: one-sided ACAT cancels exactly on bidirectional genes). The identical
  # combination is returned on both sides; two.sided tells the cascade not to
  # apply the fdr/2 direction split.
  p_two_all <- pmin(2 * pmin(ptail(t_stat, FALSE), ptail(t_stat, TRUE)), 1)
  dirs <- lapply(c(p.pos = FALSE, p.neg = TRUE), function(lower.tail) {
    p_all <- p_two_all
    p_niche <- .clampP(p_all[, cov_niche, drop = FALSE])
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

  # patient-level contrast SE: sqrt(psi * w'Vw). beta = w'alpha is recovered
  # from the coefficients outside, so only the variance needs V.
  se_pat <- if (is.null(w_rc)) rep(NA_real_, b) else {
    sqrt(scale_block * quadB)
  }
  list(t_stat = t_stat, se = se, p.pos = dirs$p.pos, p.neg = dirs$p.neg,
       se_pat = se_pat)
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
                              BPPARAM = BiocParallel::SerialParam(),
                              dispersion = c("ql", "pearson")) {
  combine <- match.arg(combine)
  dispersion <- match.arg(dispersion)
  backend <- match.arg(backend)
  W_full <- fit@W
  covtype <- as.character(fit@covtype)
  mode <- .fitMode(fit)
  # Which columns are tested depends on the design mode: the response terms
  # under a condition, the two-way CellType:niche interactions without one.
  cols_gene <- .testedCols(covtype, mode)
  Wsub <- W_full[, cols_gene, drop = FALSE]
  cov_niche <- .nicheTestCols(covtype, mode)[cols_gene]
  coefmap_sub <- fit@coefmap[cols_gene, , drop = FALSE]
  index_ct <- coefmap_sub$index[cov_niche]
  uniq_index <- unique(index_ct)

  # Patient-level contrast weights. The CellType:condition block is cell-means
  # coded (one coefficient per cell type), so the tissue-level response effect
  # is the ABUNDANCE-WEIGHTED average of those coefficients: w_c proportional to
  # the number of cells of that type, summing to 1, and zero on every other
  # tested column. Absent that block (the niche-only design) w_rc stays NULL and the
  # patient-level outputs are NA.
  cov_rc <- covtype[cols_gene] == "ResponseCellType"
  w_rc <- NULL
  if (any(cov_rc)) {
    ct_of_col <- coefmap_sub$index[cov_rc]
    n_cells <- colSums(W_full[, paste0("CellType", ct_of_col), drop = FALSE] != 0)
    w_rc <- numeric(sum(cols_gene))
    w_rc[cov_rc] <- n_cells / sum(n_cells)
  }

  # A mixed fit needs the FULL penalised information (X'WX + Lambda)^-1 rather
  # than the reduced-design covariance formed from the tested columns alone,
  # plus the working (Pearson) dispersion and a between-patient t reference.
  # `full_cov` is named for what it selects rather than for random effects.
  full_cov <- !is.null(fit@penalty)
  sel <- if (full_cov) which(cols_gene) else NULL
  penalty <- if (full_cov) fit@penalty else NULL
  df_ref <- fit@df
  # Was this fit converged per gene? Two consequences.
  #
  # (1) The SE scale. A fixed-effects fit forms se = sqrt(psi * diag(V)), which
  # is safe only while psi is edgeR's cross-gene MODERATED dispersion. The
  # convergence stage replaces it with an unshrunk per-gene profile ML, so using
  # it as the SE scale rescales every gene by the square root of a noisy
  # estimate -- measured on the toy fixture, sd(t) 2.11 -> 2.45 and max|t|
  # 12.3 -> 14.2. The Pearson working dispersion is the robust scale the mixed
  # path already uses; a polished fit takes it too.
  #
  # (2) The mean function. calculateMu() winsorises each gene's log-mu at
  # rowmedian + 4*MAD, a robustness device for fitNB's shared-weight fit. A
  # converged fit maximises the UNCLAMPED likelihood, so evaluating it at the
  # clamped mean tests a point where the score is not zero. Measured on the toy,
  # the clamp bites 4 of 20 genes and moves their Pearson dispersion by up to
  # 0.67 on the log scale. Use one mean function for both stages.
  polished <- !is.null(fit@polish)
  use_pearson <- full_cov || polished
  # 4 is calculateMu()'s own default (SpaNorm's DEFAULT_WINSOR is internal, so
  # it cannot be referenced here); fitSpiDE()'s winsor default is the same value.
  winsor_use <- if (polished) Inf else 4
  # residual df for the working-dispersion estimate (fixed columns only; the
  # penalised random columns contribute little effective df)
  disp_df <- if (use_pearson) max(nrow(W_full) - sum(!grepl("Random", covtype)), 1)

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
    p_eff <- if (full_cov) ncol(W_full) else ncol(Wsub)
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
    Wgram <- if (full_cov) {
      W_full_dev
    } else if (gpu_active) {
      SpaNorm::toGPUMatrix(Wsub, backend = backend)
    } else {
      Wsub
    }
    # genes per covariance sub-batch -- bounds the (batch, p, p) stack
    # independently of block.size, on both backends
    cov_batch <- .covBatchSize(nrow(W_full),
                               if (full_cov) ncol(W_full) else ncol(Wsub),
                               backend, gpu.mem.budget)
  }

  # Absorb the nested indicator block when there is one and we are on the CPU
  # (the absorption uses rowsum(), which has no tensor equivalent here).
  absorb <- NULL
  if (full_cov && !gpu_active && !is.null(fit@re_group)) {
    nested_cols <- !is.na(fit@re_group) & fit@re_group == "SampleCellTypeInt"
    if (any(nested_cols)) {
      xi <- which(!nested_cols)
      # every tested column is a fixed effect, so it lies in the dense block
      stopifnot(all(sel %in% xi))
      absorb <- list(solver = .newtonSolver(W_full, penalty, nested_cols),
                     sel_x = match(sel, xi))
    }
  }

  # The quasi-likelihood scale (edgeR v4 form, SpaNorm::qlDispersion): per
  # gene, the adjusted NB deviance over its effective df, moderated ACROSS
  # genes by squeezeVar. It replaces the Pearson scale for the standard
  # errors; the reference df is unchanged. Its cross-gene moderation has to
  # see every gene before any block is scored, so it is a pre-pass over the
  # same blocks, on the same backend as the inference: the deviance moments
  # are one shared (log mu, log phi) table per block and the rest is
  # elementwise over genes x cells, which is what the device does.
  use_ql <- dispersion == "ql"
  ql_scale <- NULL
  if (use_ql && !use_pearson) {
    message("dispersion = \"ql\" needs a mixed or converged fit; this ",
            "fixed-effects, unconverged fit keeps its legacy NB-dispersion scale")
    use_ql <- FALSE
  }
  if (use_ql) {
    p_fixed <- sum(!grepl("Random", covtype))
    # ONE moments table for the whole call, over a range that does not depend
    # on how the genes are blocked: the floor the means are clamped to, the
    # largest count (a fitted mean beyond it sits on the table's edge, where
    # the moments are asymptotically flat), and the genes' dispersion range.
    # Results are then exactly invariant to block.size and workers, and the
    # same table serves both backends.
    ql_table <- SpaNorm::qlMomentTable(
      lmu_range = c(log(.MU_FLOOR), log(2 * max(Y) + 2)),
      lphi_range = log(range(psi[is.finite(psi) & psi > 0])))
    ql_parts <- .bplapplySingleBLAS(blocks, function(gi) {
      Yb <- as.matrix(Y[gi, , drop = FALSE])
      if (gpu_active) {
        Yb_q <- SpaNorm::toGPUMatrix(Yb, backend = backend)
        mub <- SpaNorm::calculateMu(rep(0, length(gi)), alpha_full[gi, , drop = FALSE],
                                    W_full_dev, winsor = winsor_use, backend = backend)
        mub <- torch::torch_clamp(mub, min = .MU_FLOOR)
      } else {
        Yb_q <- Yb
        mub <- SpaNorm::calculateMu(rep(0, length(gi)), alpha_full[gi, , drop = FALSE],
                                    W_full, winsor = winsor_use)
        mub <- pmax(mub, .MU_FLOOR)
      }
      q <- SpaNorm::qlDispersion(Yb_q, mub, psi[gi], p = p_fixed, table = ql_table)
      list(s2 = q$s2, df = q$df, ave = log2(rowMeans(Yb) + 0.5))
    }, BPPARAM = BPPARAM)
    s2 <- unlist(lapply(ql_parts, `[[`, "s2")); dfq <- unlist(lapply(ql_parts, `[[`, "df"))
    ave <- unlist(lapply(ql_parts, `[[`, "ave"))
    ok <- is.finite(s2) & s2 > 0 & dfq > 0
    sq <- limma::squeezeVar(s2[ok], dfq[ok], covariate = ave[ok], robust = TRUE)
    ql_scale <- rep(NA_real_, ng); ql_scale[ok] <- sq$var.post
  }

  block_res <- .bplapplySingleBLAS(blocks, function(gi) {
    Yb <- as.matrix(Y[gi, , drop = FALSE])
    alpha_block <- alpha_full[gi, , drop = FALSE]
    psib <- psi[gi]
    zero_gmean <- rep(0, length(gi))

    if (gpu_active) {
      Yb_dev <- SpaNorm::toGPUMatrix(Yb, backend = backend)
      mub <- SpaNorm::calculateMu(zero_gmean, alpha_block, W_full_dev,
                                  winsor = winsor_use,
                                  backend = backend) # tensor, block x ncells
      inv_mu <- torch::torch_reciprocal(mub)
      wtb <- torch::torch_reciprocal(
        SpaNorm::add_vec_mat_gpu(psib, inv_mu, backend = backend))
      loglikb <- as.numeric(SpaNorm::toRMatrix(SpaNorm::rowSums_gpu(
        SpaNorm::dnbinom_gpu(Yb_dev, mu = mub, size = 1 / psib, log = TRUE))))
      if (use_pearson) {
        num <- (Yb_dev - mub)^2
        denom <- mub + SpaNorm::mult_vec_mat_gpu(psib, mub * mub,
                                                 backend = backend)
        dispb <- as.numeric(SpaNorm::toRMatrix(
          SpaNorm::rowSums_gpu(num / denom))) / disp_df
      }
      rhob <- .rhoPartialGPU(Yb_dev, mub, psib, backend)
    } else {
      mub <- SpaNorm::calculateMu(zero_gmean, alpha_block, W_full,
                                  winsor = winsor_use)
      # same floor the polish uses: an unwinsorised linear predictor can
      # underflow exp() to exactly 0, and the Pearson dispersion below would
      # then be 0/0 = NaN for that gene
      mub <- pmax(mub, .MU_FLOOR)
      wtb <- 1 / (1 / mub + psib) # nblock x ncells
      loglikb <- rowSums(dnbinom(Yb, mu = mub, size = 1 / psib, log = TRUE))
      if (use_pearson) {
        dispb <- rowSums((Yb - mub)^2 / (mub + psib * mub^2)) / disp_df
      }
      rhob <- .rhoPartial(Yb, mub, psib)
    }
    scale_block <- if (use_ql) {
      # a gene the QL pre-pass could not score keeps its Pearson scale
      sb <- ql_scale[gi]; sb[is.na(sb)] <- dispb[is.na(sb)]; sb
    } else if (use_pearson) dispb else psib

    if (combine == "cauchy") {
      res <- .waldCauchyBlock(alpha_block[, cols_gene, drop = FALSE], Wsub,
                              wtb, scale_block, cov_niche, index_ct,
                              uniq_index, Wgram, w_rc = w_rc,
                              W_full = if (full_cov) W_full else NULL,
                              penalty = penalty, sel = sel, df = df_ref,
                              backend = backend, cov.batch = cov_batch,
                              absorb = absorb)
      res$loglik <- loglikb
      res$rho_s <- rhob$s
      res$rho_n <- rhob$n
      return(res)
    }

    # Brown's method: per-gene, needs plain-R weights and design
    wtb_r <- if (gpu_active) SpaNorm::toRMatrix(wtb) else wtb
    per_gene <- lapply(seq_along(gi), function(i) {
      .waldBrownGene(alpha_block[i, cols_gene], Wsub, wtb_r[i, ],
                     scale_block[i], cov_niche, index_ct, uniq_index,
                     combine = "brown", W_full = if (full_cov) W_full else NULL,
                     penalty = penalty, sel = sel, df = df_ref, w_rc = w_rc)
    })
    list(
      rho_s = rhob$s, rho_n = rhob$n,
      t_stat = do.call(rbind, lapply(per_gene, `[[`, "t_stat")),
      se = do.call(rbind, lapply(per_gene, `[[`, "se")),
      p.pos = do.call(rbind, lapply(per_gene, `[[`, "p.pos")),
      p.neg = do.call(rbind, lapply(per_gene, `[[`, "p.neg")),
      se_pat = vapply(per_gene, `[[`, numeric(1), "se_pat"),
      loglik = loglikb
    )
  }, BPPARAM = BPPARAM)

  bind <- function(field) do.call(rbind, lapply(block_res, `[[`, field))
  t_stat <- bind("t_stat")
  se <- bind("se")
  p.pos <- bind("p.pos")
  p.neg <- bind("p.neg")
  loglik <- unlist(lapply(block_res, `[[`, "loglik"), use.names = FALSE)
  # w_rc is NULL when the design carries no CellType:condition block (the
  # niche-only mode). The blocks then return NA for every gene; collapse that
  # to the empty vector the slot documents, so `length(@se_patient) > 0` stays
  # a reliable test for "this fit has a patient-level contrast" rather than
  # reporting a full-length vector of NAs.
  se_pat <- if (is.null(w_rc)) {
    numeric(0)
  } else {
    unlist(lapply(block_res, `[[`, "se_pat"), use.names = FALSE)
  }
  # Inter-gene correlation, reduced from the per-block partial sums. Blocks
  # contribute additively, so this is exact and independent of block size.
  rs <- lapply(block_res, `[[`, "rho_s")
  rs <- rs[!vapply(rs, is.null, logical(1))]
  if (length(rs)) {
    rn <- sum(unlist(lapply(block_res, `[[`, "rho_n"), use.names = FALSE))
    fit@rho <- .meanCorFromZ(Reduce(`+`, rs), nrow(fit@W), rn)
  }

  gnames <- rownames(alpha_full)
  rownames(t_stat) <- rownames(se) <- gnames
  rownames(p.pos) <- rownames(p.neg) <- gnames
  names(loglik) <- gnames
  if (length(se_pat) == length(gnames)) names(se_pat) <- gnames else se_pat <- numeric(0)

  fit@two.sided <- identical(combine, "cauchy")
  fit@se_patient <- se_pat
  fit@t_stat <- t_stat
  fit@se <- se
  fit@p.combined.pos <- p.pos
  fit@p.combined.neg <- p.neg
  fit@loglik <- loglik
  fit
}
