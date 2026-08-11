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
