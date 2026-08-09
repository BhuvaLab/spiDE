# Mixed-effects machinery: the PQL cell subsampling, the Schall
# variance-component loop, and the degrees-of-freedom calculations that depend
# on it. Split out of fitSpiDE.R, which keeps the fit orchestration -- these
# are the parts a reader only needs when `random != "none"`.
#
# Random effects are implemented via the ridge = random-effects equivalence:
# the random-effect columns are L2-penalised at fit time with
# `lambda.a = 1/tau2`, so no SpaNorm change is needed.

#' Draw a stratified cell subsample for the PQL loop
#'
#' Samples cells per \code{(cell_type, sample)} stratum so every sample and cell
#' type stays represented (the random-effect BLUPs need each sample populated).
#' For a stratum of size \code{n} the number kept is
#' \code{min(n, max(ceil(prop * n), min.cells))}: strata at or below
#' \code{min.cells} are taken whole, otherwise the larger of \code{prop} of the
#' cells and \code{min.cells}. No seed is set here; reproducibility is the
#' caller's responsibility via an external \code{set.seed()}.
#'
#' @param cell_type,sample per-cell cell-type and sample labels (length ncells).
#' @param prop the sampling proportion (in (0, 1]); \code{>= 1} keeps all cells.
#' @param min.cells the per-stratum floor.
#' @return a logical vector (length ncells) marking the sampled cells.
#' @noRd
.stratifiedCellIdx <- function(cell_type, sample, prop, min.cells = 100L) {
  n <- length(cell_type)
  if (is.null(prop) || prop >= 1) {
    return(rep(TRUE, n))
  }
  strata <- interaction(as.character(sample), as.character(cell_type),
                        drop = TRUE)
  idx <- logical(n)
  for (s in levels(strata)) {
    pos <- which(strata == s)
    ns <- length(pos)
    k <- min(ns, max(ceiling(prop * ns), min.cells))
    sel <- if (k >= ns) pos else sample(pos, k)
    idx[sel] <- TRUE
  }
  idx
}

#' Gene-averaged working weights over all cells, computed gene-block-wise
#'
#' The representative w-bar_c = mean_g mu_{c,g}/(1 + psi_g mu_{c,g}) used to form
#' the representative penalised inverse for the Satterthwaite df. Blocked over
#' genes so the genes x cells mean matrix is never fully materialised.
#' @noRd
.repWeights <- function(Y, alpha, W, psi, block.size = 2000L) {
  ng <- nrow(alpha)
  acc <- numeric(ncol(Y))
  for (gi in .chunkGenes(ng, block.size)) {
    mub <- SpaNorm::calculateMu(rep(0, length(gi)), alpha[gi, , drop = FALSE], W)
    acc <- acc + colSums(mub / (1 + psi[gi] * mub))
  }
  acc / ng
}

#' Estimate random-effect variance components (Schall / PQL, shared across genes)
#'
#' Implements the mixed model via ridge: random-effect columns are penalised by
#' \code{1/tau2}, and the variance components \code{tau2} are estimated by
#' alternating (i) a penalised NB-GLM fit (\code{SpaNorm::fitNB} with a per-column
#' \code{lambda.a}) and (ii) a Schall variance-component update using the pooled
#' BLUPs and the block effective degrees of freedom. A single \code{tau2} per
#' random-effect group is shared across genes (matching the dispersion
#' moderation), estimated with a representative (gene-averaged) working weight.
#'
#' For speed the inner loop fits on a stratified cell subsample (\code{idx}) with
#' a single dispersion iteration (\code{re.maxit.psi}); the variance components
#' being shared scalars do not need every cell. Once \code{tau2} converges a
#' single \strong{final fit} is run on all cells with full dispersion, and its
#' coefficients / dispersion are what inference uses.
#'
#' @param Y counts (genes x cells).
#' @param W the full design (fixed + random columns).
#' @param re_group per-column random-effect group (NA for fixed columns).
#' @param lambda.a base ridge penalty for the fixed columns (usually 0).
#' @param re.maxit,re.tol outer-loop iteration cap and relative tolerance (the
#'   tolerance is a relative change on \code{log(tau2)}).
#' @param tau2.init initial variance component.
#' @param idx logical over cells: the subsample used for the inner loop fits
#'   (NULL = all cells). The final fit always uses all cells.
#' @param re.maxit.psi dispersion iterations (\code{maxit.psi}) for the inner
#'   loop fits; the final fit uses full dispersion.
#' @param df.method one of "satterthwaite" (the default: a per-tested-column df
#'   vector via \code{.satterthwaiteDF()}) or "between" (a scalar \code{S - 2}
#'   residual df, the back-compatible behaviour).
#' @param cols_tested a logical over \code{colnames(W)} marking the
#'   Response/ResponseNiche columns needing a df (only used when
#'   \code{df.method == "satterthwaite"}).
#' @param mode the design mode, "condition" or "niche"; selects the
#'   \code{df.method = "between"} reference df (see the comment at its
#'   computation).
#' @return a list with the fitNB result plus \code{penalty} (per-column
#'   \code{lambda.a}), \code{tau2} (named variance components) and \code{df}
#'   (effective residual degrees of freedom: a scalar under "between", a
#'   named length-k vector aligned to \code{cols_tested} under
#'   "satterthwaite").
#' @importFrom stats setNames
#' @noRd
.fitNBmixed <- function(Y, W, re_group, lambda.a, winsor, backend, verbose,
                        re.maxit = 2L, re.tol = 1e-3, tau2.init = 1,
                        tau2.range = c(1e-8, 1e4), idx = NULL,
                        re.maxit.psi = 1L, df.method = "satterthwaite",
                        cols_tested = NULL, mode = "condition", ...) {
  p <- ncol(W)
  groups <- unique(re_group[!is.na(re_group)])
  base <- if (length(lambda.a) == 1) rep(lambda.a, p) else lambda.a
  base[!is.na(re_group)] <- 0
  tau2 <- stats::setNames(rep(tau2.init, length(groups)), groups)
  if (is.null(idx)) idx <- rep(TRUE, ncol(Y))
  Wi <- W[idx, , drop = FALSE]

  # the inner loop forces maxit.psi = re.maxit.psi; strip any user maxit.psi so
  # it does not collide, but forward it (if given) to the full final fit.
  dots <- list(...)
  dots_loop <- dots
  dots_loop$maxit.psi <- NULL

  pen <- base
  for (it in seq_len(re.maxit)) {
    pen <- base
    for (g in groups) pen[which(re_group == g)] <- 1 / tau2[[g]]
    tau2_fit <- tau2

    if (verbose) message(sprintf("  RE iteration %d: %s", it,
      paste(sprintf("%s=%.3g", groups, tau2), collapse = ", ")))
    fit <- do.call(SpaNorm::fitNB, c(
      list(Y, W, idx = idx, lambda.a = pen, winsor = winsor,
           backend = backend, maxit.psi = re.maxit.psi, verbose = FALSE),
      dots_loop))
    alpha <- fit$alpha

    # representative (gene-averaged) NB working weight per cell, shared info,
    # computed on the SAME sampled cells the fit used (Wi) so the Schall update
    # is internally consistent with the penalised fit.
    mu <- SpaNorm::calculateMu(rep(0, nrow(alpha)), alpha, Wi)
    wbar <- colMeans(mu / (1 + fit$psi * mu))
    info <- crossprod(Wi * sqrt(wbar))
    minv <- SpaNorm::invert_mat(info + diag(pen))

    tau2_new <- tau2
    ng <- nrow(alpha)
    for (g in groups) {
      gi <- which(re_group == g)
      trc <- sum(diag(minv[gi, gi, drop = FALSE]))
      edf <- max(length(gi) - (1 / tau2[[g]]) * trc, 1e-6)
      b2 <- sum(alpha[, gi, drop = FALSE]^2)
      # keep the random-effect columns at least weakly penalised: sample-level
      # covariates and per-sample columns are collinear in the tau2 -> Inf limit
      tau2_new[[g]] <- min(max(b2 / (ng * edf), tau2.range[1]), tau2.range[2])
    }
    # relative change on log(tau2): natural for a scale parameter and robust to
    # the slow monotone decay of the slope variance
    converged <- max(abs(log(tau2_new) - log(tau2))) < re.tol
    tau2 <- tau2_fit
    if (converged) break
    tau2 <- tau2_new
  }

  # final fit: ALL cells, FULL dispersion, at the converged penalty. This is the
  # fit inference uses (alpha / psi / gmean); the loop only supplied tau2.
  pen <- base
  for (g in groups) pen[which(re_group == g)] <- 1 / tau2[[g]]
  fit <- do.call(SpaNorm::fitNB, c(
    list(Y, W, lambda.a = pen, winsor = winsor, backend = backend,
         verbose = verbose),
    dots))

  # Reference df for the Wald t-test under df.method = "between".
  #
  # condition mode: the ResponseNiche coefficient is a DIFFERENCE BETWEEN
  #   within-sample niche slopes across conditions, so the replication unit for
  #   that comparison is the patient -> S - 2 (the condition has two levels).
  #   This, with the working (Pearson) dispersion used at inference time, is
  #   what corrects the cell-level pseudo-replication.
  # niche mode, random = "intercept": there is no between-condition comparison
  #   left. The CellType:niche slope is identified by niche density varying
  #   cell-to-cell INSIDE each sample, so cells are the replicates and the
  #   reference is the residual df ncells - p_fixed. The option name "between"
  #   is a misnomer in this case; it is kept for back-compatibility.
  # niche mode, random = "slope": the per-sample random slopes sit on the very
  #   columns being tested, moving that contrast back into the between-sample
  #   stratum -> S - 1 (no condition contrast spends a df here).
  #
  # df.method = "satterthwaite" instead derives a per-tested-column df from the
  # shared variance-component fit (see .satterthwaiteDF()), which computes this
  # same distinction rather than hard-coding it -- and is the default.
  n_samples <- sum(re_group == "SampleInt", na.rm = TRUE)
  has_slope <- any(re_group == "SampleSlope", na.rm = TRUE)
  df_between <- if (identical(mode, "niche")) {
    if (has_slope) {
      max(n_samples - 1, 1)
    } else {
      max(ncol(Y) - sum(is.na(re_group)), 1)
    }
  } else {
    max(n_samples - 2, 1)
  }
  if (df.method == "satterthwaite" && !is.null(cols_tested)) {
    wbar <- .repWeights(Y, fit$alpha, W, fit$psi)
    A <- crossprod(W * sqrt(wbar))
    minv <- SpaNorm::invert_mat(A + diag(pen))
    tested <- which(cols_tested)
    df <- .satterthwaiteDF(A, minv, pen, re_group, tau2, tested, ncol(Y),
                          colnames(W)[tested])
  } else {
    df <- df_between
  }

  c(fit, list(penalty = pen, tau2 = tau2, df = df))
}

#' Covariance of the working-model variance parameters (phi, tau2_groups)
#'
#' Reduced-form REML expected (Fisher) information at phi = 1 for the working
#' linear mixed model V = W-bar^{-1} + Z G Z', built entirely from p x p M^{-1}
#' sub-blocks (never an n x n matrix). Parameter order is (phi, groups...). The
#' phi-phi entry reduces to 0.5*(ncells - p + tr((M^{-1} Lambda)^2)); the
#' group entries use `B_{a,b} = A[ca,cb] - A[ca,] M^{-1} A[,cb] = Z_a' P Z_b`.
#' @noRd
.varParamCov <- function(A, minv, pen, re_group, ncells) {
  p <- ncol(A)
  groups <- unique(re_group[!is.na(re_group)])
  Mn <- length(groups)
  gcols <- lapply(groups, function(gr) which(re_group == gr))
  I <- matrix(0, Mn + 1L, Mn + 1L)
  ML <- minv * rep(pen, each = p)          # M^{-1} Lambda (scale columns by pen)
  I[1, 1] <- 0.5 * (ncells - p + sum(ML * t(ML)))   # tr((M^{-1}Lambda)^2)
  minvLminv <- ML %*% minv                 # M^{-1} Lambda M^{-1}
  for (a in seq_len(Mn)) {
    ca <- gcols[[a]]
    Aa <- A[, ca, drop = FALSE]            # U_a = A[, group a]
    MinvAa <- minv %*% Aa
    for (b in seq_len(a)) {
      cb <- gcols[[b]]
      Bab <- A[ca, cb, drop = FALSE] - crossprod(Aa, minv %*% A[, cb, drop = FALSE])
      I[a + 1L, b + 1L] <- I[b + 1L, a + 1L] <- 0.5 * sum(Bab * Bab)  # 0.5||B||_F^2
    }
    tr_Baa <- sum(diag(A)[ca]) - sum(MinvAa * Aa)
    tr_UMU <- sum((minvLminv %*% Aa) * Aa)
    I[1, a + 1L] <- I[a + 1L, 1] <- 0.5 * (tr_Baa - tr_UMU)
  }
  cov <- tryCatch(solve(I),
                  error = function(e) SpaNorm::invert_mat(I))
  list(cov = cov, groups = groups, gcols = gcols)
}

#' Per-coefficient Satterthwaite degrees of freedom (shared across genes)
#'
#' df_j = 2 v_jj^2 / (d_j' Cov(theta-hat) d_j), with v_jj = (M^{-1})_jj, gradient
#' d_j = (v_jj, g_j1, ..., g_jM) and g_jm = (pen_m / tau2_m) * sum_{k in m} M^{-1}_{jk}^2.
#' Invariant to the per-gene dispersion phi (spec), hence one shared vector.
#' @noRd
.satterthwaiteDF <- function(A, minv, pen, re_group, tau2, tested, ncells,
                             tested_names) {
  vp <- .varParamCov(A, minv, pen, re_group, ncells)
  Vt  <- minv[tested, , drop = FALSE]      # k x p
  vjj <- diag(minv)[tested]                # k
  grad <- matrix(0, length(tested), length(vp$groups) + 1L)
  grad[, 1] <- vjj
  for (a in seq_along(vp$groups)) {
    gc <- vp$gcols[[a]]
    pen_m <- pen[gc[1]]                     # constant within a group
    Sm <- rowSums(Vt[, gc, drop = FALSE]^2)
    grad[, a + 1L] <- (pen_m / tau2[[vp$groups[a]]]) * Sm
  }
  varse2 <- rowSums((grad %*% vp$cov) * grad)
  df <- 2 * vjj^2 / varse2
  df <- pmin(pmax(df, 1), ncells)
  stats::setNames(df, tested_names)
}
