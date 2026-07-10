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
    ptail <- function(t, lower.tail) stats::pnorm(t, lower.tail = lower.tail)
  } else {
    # mixed-effects fit: full penalised information (X'WX + Lambda), then the
    # fixed-effect (Response/ResponseNiche) block -> larger, between-sample SEs
    info <- crossprod(W_full * wt_g, W_full) + diag(penalty)
    varcov <- SpaNorm::invert_mat(info)[sel, sel, drop = FALSE]
    ptail <- function(t, lower.tail) stats::pt(t, df = df, lower.tail = lower.tail)
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

#' Compute Wald inference + within-gene combination for a SpiDEFit, block-wise
#'
#' Fills the \code{t_stat}, \code{se}, \code{p.combined.pos}, \code{p.combined.neg}
#' and (recomputed) \code{loglik} slots. Genes are processed in blocks so the
#' genes x cells fitted-mean / weight matrices are never fully materialised, and
#' blocks are dispatched with BiocParallel.
#'
#' @param fit a SpiDEFit (with the NB fit populated).
#' @param Y the counts matrix (dense, sparse, or DelayedArray).
#' @param block.size a numeric, genes per block (NULL = a single block).
#' @param BPPARAM a BiocParallelParam.
#' @return the input \code{fit} with inference slots populated.
#' @importFrom stats dnbinom
#' @importFrom BiocParallel bplapply SerialParam
#' @noRd
.blockedInference <- function(fit, Y, block.size = NULL,
                              combine = c("cauchy", "brown"),
                              BPPARAM = BiocParallel::SerialParam()) {
  combine <- match.arg(combine)
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
  alpha_sub <- alpha_full[, cols_gene, drop = FALSE]
  psi <- fit@psi
  ng <- nrow(alpha_full)

  blocks <- .chunkGenes(ng, block.size)

  block_res <- BiocParallel::bplapply(blocks, function(gi) {
    Yb <- as.matrix(Y[gi, , drop = FALSE])
    mub <- SpaNorm::calculateMu(rep(0, length(gi)),
                                alpha_full[gi, , drop = FALSE], W_full)
    psib <- psi[gi]
    wtb <- 1 / (1 / mub + psib) # nblock x ncells
    loglikb <- rowSums(dnbinom(Yb, mu = mub, size = 1 / psib, log = TRUE))
    # per-gene working (Pearson) dispersion, used to scale the SEs when the
    # random effects are present (else the NB dispersion psi is used)
    if (has_re) {
      dispb <- rowSums((Yb - mub)^2 / (mub + psib * mub^2)) / disp_df
    }

    per_gene <- lapply(seq_along(gi), function(i) {
      scale_i <- if (has_re) dispb[i] else psib[i]
      .waldBrownGene(alpha_sub[gi[i], ], Wsub, wtb[i, ], scale_i,
                     cov_niche, index_ct, uniq_index, combine = combine,
                     W_full = if (has_re) W_full else NULL,
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
