# Stage-1 machinery for the two-stage estimator: SpaNorm fit extraction, the
# one-step working response, the joint weighted slope fit, and the basis-R2
# diagnostic. See docs/superpowers/specs/2026-08-11-twostage-fixes-design.md.

#' Extract and validate the stored SpaNorm fit
#'
#' The fit supplies the biology/LS design split (`wtype`), the coefficients
#' and the dispersion that stage 1 linearises around. Kept as small pieces --
#' full genes x cells matrices are only materialised per (sample, index)
#' subset by .stage1Epsilon().
#' @param spe a SpatialExperiment previously normalised with SpaNorm().
#' @return list(alpha, gmean, psi, W, bio) with `bio` a logical over columns
#'   of `W` marking the biology block.
#' @noRd
.spanormComponents <- function(spe) {
  fit <- S4Vectors::metadata(spe)$SpaNorm
  if (is.null(fit)) {
    stop("no SpaNorm fit found in metadata(spe)$SpaNorm; run ",
         "SpaNorm::SpaNorm(spe) first (stage1 = \"spanorm\" needs it), or ",
         "use stage1 = \"ols\"")
  }
  if (fit$ncells != ncol(spe) || fit$ngenes != nrow(spe)) {
    stop("the stored SpaNorm fit (", fit$ngenes, " x ", fit$ncells,
         ") does not match spe (", nrow(spe), " x ", ncol(spe), "); ",
         "re-run SpaNorm::SpaNorm() on this object")
  }
  wt <- as.character(fit$wtype)
  list(alpha = fit$alpha, gmean = fit$gmean, psi = fit$psi, W = fit$W,
       bio = wt == "biology")
}

#' One-step working response and weights for a cell subset
#'
#' Linearises log counts at the FULL fitted mean (the correct expansion
#' point: handles zeros, no baseline leakage) and, under "addback", returns
#' the biology component so that, net, only the LS (and batch) effect is
#' removed. Materialises genes x subset matrices only.
#' @param Y counts matrix (genes x all cells; dense).
#' @param comp the .spanormComponents() list.
#' @param cells integer/logical index of the subset's cells.
#' @param epsilon "addback" (default; eta_bio + z) or "residual" (z).
#' @return list(eps, w): the response and NB working weights, genes x cells.
#' @noRd
.stage1Epsilon <- function(Y, comp, cells, epsilon = c("addback", "residual")) {
  epsilon <- match.arg(epsilon)
  Wc <- comp$W[cells, , drop = FALSE]
  eta <- comp$gmean + tcrossprod(comp$alpha, Wc)
  mu <- exp(eta)
  z <- (Y[, cells, drop = FALSE] - mu) / mu
  w <- mu / (1 + comp$psi * mu)
  eps <- if (epsilon == "addback") {
    comp$gmean + tcrossprod(comp$alpha[, comp$bio, drop = FALSE],
                            Wc[, comp$bio, drop = FALSE]) + z
  } else {
    z
  }
  list(eps = eps, w = w)
}

#' Joint WLS of the working response on all niche columns, per gene
#'
#' Design [1, X] (the intercept absorbs each gene's level, so no per-gene
#' weighted centring is needed). Per-gene weights preclude one shared Gram
#' matrix, so cross-products are batched as matrix products over the
#' upper-triangle column pairs, then each gene's small (k+1) system is
#' solved in a loop. Degenerate columns (weighted variance ~ 0) are dropped
#' up front and return NA. Variance is the usual WLS sandwich-free form
#' phi_g * (X'W_gX)^{-1} with phi_g the weighted RSS over n - rank.
#' @param eps,w genes x cells response and weights (Task 3).
#' @param X cells x k centred niche matrix.
#' @return list(beta, var): genes x k, NA for dropped columns.
#' @noRd
.jointSlopes <- function(eps, w, X) {
  G <- nrow(eps); n <- ncol(eps); k <- ncol(X)
  out_b <- out_v <- matrix(NA_real_, G, k,
                           dimnames = list(rownames(eps), colnames(X)))
  keep <- apply(X, 2, function(x) stats::var(x) > 1e-10)
  if (!any(keep)) return(list(beta = out_b, var = out_v))
  D <- cbind(`(Intercept)` = 1, X[, keep, drop = FALSE])
  p <- ncol(D)
  # per-gene Gram entries: (X'W_gX)_{ij} = sum_c w_gc D_ci D_cj
  ut <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  P <- D[, ut[, 1], drop = FALSE] * D[, ut[, 2], drop = FALSE]  # n x npairs
  Gm <- w %*% P                                                  # G x npairs
  Rhs <- (w * eps) %*% D                                         # G x p
  wrss_tot <- rowSums(w * eps^2)
  for (g in seq_len(G)) {
    A <- matrix(0, p, p)
    A[cbind(ut[, 1], ut[, 2])] <- Gm[g, ]
    A[cbind(ut[, 2], ut[, 1])] <- Gm[g, ]
    Ai <- tryCatch(solve(A), error = function(e) NULL)
    if (is.null(Ai)) next
    bg <- Ai %*% Rhs[g, ]
    df <- max(n - p, 1)
    phi <- max((wrss_tot[g] - sum(bg * Rhs[g, ])) / df, 1e-12)
    out_b[g, keep] <- bg[-1]
    out_v[g, keep] <- phi * diag(Ai)[-1]
  }
  list(beta = out_b, var = out_v)
}

#' R2 of each niche column on the biology basis, within a subset
#'
#' Measures the attenuation the "residual" response would suffer and the
#' smooth-trend overlap the "addback" response is exposed to. Reported per
#' (sample, index) subset in the diagnostics.
#' @noRd
.nicheBasisR2 <- function(X, B) {
  f <- stats::lm.fit(cbind(1, B), X)
  res <- as.matrix(f$residuals)
  tss <- colSums(sweep(X, 2, colMeans(X))^2)
  r2 <- 1 - colSums(res^2) / pmax(tss, 1e-12)
  stats::setNames(pmin(pmax(r2, 0), 1), colnames(X))
}
