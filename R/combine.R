# Combine p-values across niche bandwidths with the Cauchy combination test,
# weighted by each gene's relative log-likelihood across bandwidths. Reproduces
# cauchy_combination() and the gene.w construction in YTMA_nicheDE_v9.md.

#' Weighted Cauchy combination of p-values across columns
#'
#' @param p a genes x k matrix of p-values (columns = bandwidths).
#' @param w a genes x k matrix of weights (default equal weights).
#' @return a numeric vector (length genes) of combined p-values.
#' @importFrom stats pcauchy weighted.mean
#' @noRd
.cauchyCombine <- function(p, w = NULL) {
  p <- as.matrix(p)
  if (is.null(w)) {
    w <- matrix(1, nrow(p), ncol(p))
  }
  tp <- tan((0.5 - p) * pi)
  stat <- vapply(seq_len(nrow(tp)), function(i) {
    stats::weighted.mean(tp[i, ], w[i, ], na.rm = TRUE)
  }, numeric(1))
  1 - stats::pcauchy(stat)
}

#' Per-gene, per-bandwidth Cauchy-combination weights
#'
#' Relative log-likelihood weights: \code{exp(loglik - max)} across bandwidths,
#' with small weights thresholded to zero.
#'
#' @param fits a list of SpiDEFit objects.
#' @param thresh a numeric, weights below this are set to 0.
#' @return a genes x bandwidth weight matrix.
#' @importFrom matrixStats rowMaxs
#' @noRd
.geneWeights <- function(fits, thresh = 0.1) {
  ll <- do.call(cbind, lapply(fits, function(f) f@loglik))
  w <- exp(ll - matrixStats::rowMaxs(ll))
  w[w < thresh] <- 0
  w
}

#' Combine each column of the per-bandwidth Brown p-values across bandwidths
#'
#' @param fits a list of SpiDEFit objects.
#' @param slot one of "p.brown.pos" or "p.brown.neg".
#' @param gene.w the gene weight matrix from .geneWeights().
#' @return a genes x (1 + n_index) matrix of combined p-values.
#' @noRd
.combineBandwidths <- function(fits, slot, gene.w) {
  ref <- methods::slot(fits[[1]], slot)
  out <- ref
  for (j in seq_len(ncol(ref))) {
    pj <- vapply(fits, function(f) methods::slot(f, slot)[, j], numeric(nrow(ref)))
    out[, j] <- .cauchyCombine(pj, gene.w)
  }
  out
}
