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
