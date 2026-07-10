# Wald + Brown inference. For each gene the working weights, Wald covariance,
# standard errors and t-statistics of the Response / ResponseNiche coefficients
# are computed, then per-covariate p-values are combined across correlated
# ResponseNiche covariates with Brown's method (poolr) into a gene-level and a
# per-index-cell-type p-value, separately for the up and down directions.
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

#' Wald + Brown inference for a single gene
#'
#' @param alpha_g numeric coefficients for the Response/ResponseNiche columns.
#' @param Wsub the Response/ResponseNiche sub-design (cells x k).
#' @param wt_g numeric working weights for this gene (length ncells).
#' @param psi_g the gene's dispersion.
#' @param cov_niche logical marking the ResponseNiche columns within the subset.
#' @param index_ct index cell type of each ResponseNiche column.
#' @param uniq_index the unique index cell types (column order of the output).
#' @return a list with t_stat, se (length k) and p.pos, p.neg (length
#'   1 + n_index: gene-level then per-index-cell-type).
#' @importFrom stats pnorm cov2cor setNames
#' @noRd
.waldBrownGene <- function(alpha_g, Wsub, wt_g, psi_g, cov_niche, index_ct, uniq_index) {
  varcov <- SpaNorm::invert_mat(crossprod(Wsub * wt_g, Wsub))
  se <- sqrt(psi_g * diag(varcov))
  t_stat <- alpha_g / se
  names(se) <- names(t_stat) <- colnames(Wsub)

  # correlation structure of the ResponseNiche coefficients -> Brown weights
  vc <- varcov[cov_niche, cov_niche, drop = FALSE]
  a <- diag(1 / diag(vc), nrow = nrow(vc))
  vc <- a %*% vc %*% t(a)
  vc <- stats::cov2cor(vc)
  vc <- poolr::mvnconv(vc, target = "m2lp", side = 1, cov2cor = FALSE)

  dirs <- lapply(c(p.pos = FALSE, p.neg = TRUE), function(lower.tail) {
    x <- stats::pnorm(t_stat, lower.tail = lower.tail)
    xn <- x[cov_niche]
    p_gene <- poolr::fisher(xn, side = 1, R = vc, adjust = "generalized")$p
    p_ct <- vapply(uniq_index, function(ct) {
      is_ct <- index_ct == ct
      poolr::fisher(xn[is_ct], side = 1, R = vc[is_ct, is_ct, drop = FALSE],
                    adjust = "generalized")$p
    }, numeric(1))
    c(Gene = p_gene, p_ct)
  })

  list(t_stat = t_stat, se = se, p.pos = dirs$p.pos, p.neg = dirs$p.neg)
}

#' Compute Wald + Brown inference for a SpiDEFit, block-wise
#'
#' Fills the \code{t_stat}, \code{se}, \code{p.brown.pos}, \code{p.brown.neg}
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
                              BPPARAM = BiocParallel::SerialParam()) {
  W_full <- fit@W
  covtype <- as.character(fit@covtype)
  cols_gene <- grepl("Response", covtype) # Response + ResponseNiche columns
  Wsub <- W_full[, cols_gene, drop = FALSE]
  cov_niche <- covtype[cols_gene] == "ResponseNiche"
  coefmap_sub <- fit@coefmap[cols_gene, , drop = FALSE]
  index_ct <- coefmap_sub$index[cov_niche]
  uniq_index <- unique(index_ct)

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

    per_gene <- lapply(seq_along(gi), function(i) {
      .waldBrownGene(alpha_sub[gi[i], ], Wsub, wtb[i, ], psib[i],
                     cov_niche, index_ct, uniq_index)
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
  fit@p.brown.pos <- p.pos
  fit@p.brown.neg <- p.neg
  fit@loglik <- loglik
  fit
}
