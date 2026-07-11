# Per-bandwidth model fitting. For each niche bandwidth a design matrix is
# built and a per-gene negative binomial GLM is fit over ALL genes at once using
# SpaNorm::fitNB (dispersion is moderated across genes by edgeR, so genes must
# not be blocked at fit time). The gene-independent Wald/Brown inference is
# added by testSpiDE()/.blockedInference() (see inference.R).

# resolve the bandwidth grid from stored niche reducedDims when not supplied
.detectSigma <- function(spe, name = "Niche") {
  rdn <- SingleCellExperiment::reducedDimNames(spe)
  hits <- grep(sprintf("^%s[0-9.]+$", name), rdn, value = TRUE)
  if (length(hits) == 0) {
    stop("no niche reducedDims found; run buildNiches() or supply 'sigma'")
  }
  sort(as.numeric(sub(sprintf("^%s", name), "", hits)))
}

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
#' @return a list with the fitNB result plus \code{penalty} (per-column
#'   \code{lambda.a}), \code{tau2} (named variance components) and \code{df}
#'   (effective residual degrees of freedom).
#' @importFrom stats setNames
#' @noRd
.fitNBmixed <- function(Y, W, re_group, lambda.a, winsor, backend, verbose,
                        re.maxit = 10L, re.tol = 1e-3, tau2.init = 1,
                        tau2.range = c(1e-8, 1e4), idx = NULL,
                        re.maxit.psi = 1L, ...) {
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

  # between-patient reference df for the Wald t-test: the response contrasts live
  # in the between-sample error stratum, so they are tested against the number of
  # samples, not cells (the condition has two levels -> S - 2). This, together
  # with the working (Pearson) dispersion used at inference time, is what
  # corrects the cell-level pseudo-replication.
  n_samples <- sum(re_group == "SampleInt", na.rm = TRUE)
  df <- max(n_samples - 2, 1)

  c(fit, list(penalty = pen, tau2 = tau2, df = df))
}

#' Fit the spiDE negative binomial GLM for one bandwidth
#'
#' @return a SpiDEFit with the fit populated and inference slots empty.
#' @importFrom SummarizedExperiment assay
#' @importFrom stats dnbinom
#' @noRd
.fitOneBandwidth <- function(Y, spe, condition, sigma, index, niche, covariates,
                             cell_type, winsor, lambda.a, backend, name,
                             verbose, sample_id = "sample_id", random = "none",
                             re.maxit = 10L, re.tol = 1e-3, tau2.init = 1,
                             re.prop = 0.1, re.maxit.psi = 1L,
                             re.min.cells = 100L, ...) {
  des <- .buildNicheDesign(spe, condition, sigma, index, niche, covariates,
                           cell_type, name, sample_id, random)
  W <- des$W

  # fit all genes at once (dispersion moderated across genes). With random
  # effects, an outer Schall/PQL loop estimates the variance components.
  if (random == "none") {
    fit <- SpaNorm::fitNB(Y, W, lambda.a = lambda.a, winsor = winsor,
                          backend = backend, verbose = verbose, ...)
    penalty <- NULL
    tau2 <- NULL
    df <- NULL
  } else {
    cd <- SummarizedExperiment::colData(spe)
    idx <- .stratifiedCellIdx(cd[[cell_type]], cd[[sample_id]], re.prop,
                              re.min.cells)
    fit <- .fitNBmixed(Y, W, des$re_group, lambda.a, winsor, backend, verbose,
                       re.maxit = re.maxit, re.tol = re.tol,
                       tau2.init = tau2.init, idx = idx,
                       re.maxit.psi = re.maxit.psi, ...)
    penalty <- fit$penalty
    tau2 <- fit$tau2
    df <- fit$df
  }
  alpha <- fit$alpha
  rownames(alpha) <- rownames(Y)
  colnames(alpha) <- colnames(W)

  # per-gene log-likelihood for Cauchy weighting (recomputed from the fit; the
  # fitNB $loglik is per-iteration, not per-gene). Computed gene-block-wise so
  # the whole counts matrix is never densified (the invariant); .blockedInference
  # recomputes this too, so this value is only used if inference is skipped.
  loglik <- .blockLoglik(Y, alpha, W, fit$psi)

  new(
    "SpiDEFit",
    sigma = sigma,
    ngenes = nrow(alpha),
    ncells = ncol(Y),
    W = W,
    covtype = des$covtype,
    coefmap = S4Vectors::DataFrame(des$coefmap),
    alpha = alpha,
    gmean = as.numeric(fit$gmean),
    psi = as.numeric(fit$psi),
    loglik = as.numeric(loglik),
    re_group = if (random == "none") NULL else des$re_group,
    tau2 = tau2,
    penalty = penalty,
    df = df,
    t_stat = NULL,
    se = NULL,
    p.combined.pos = NULL,
    p.combined.neg = NULL,
    sampling = fit$sampling
  )
}

#' Fit the spiDE negative binomial model
#'
#' Builds the neighbourhood-interaction design for each niche bandwidth and fits
#' a per-gene negative binomial GLM over all genes using the SpaNorm
#' \code{\link[SpaNorm]{fitNB}} engine. The heavy fit uses fitNB's compute
#' \code{backend} (CPU or GPU); the counts \code{Y} may be dense, sparse, or a
#' \code{DelayedArray} and is not densified up front. The returned object holds
#' one \code{SpiDEFit} per bandwidth. Neighbourhood effects are tested
#' separately with [testSpiDE()].
#'
#' @param spe a SpatialExperiment with niche reducedDims (see [buildNiches()]).
#' @param condition a character, the colData column of the tested condition
#'   (must have exactly two levels).
#' @param index,niche character vectors restricting the index / niche cell
#'   types considered (NULL = all).
#' @param covariates a character vector of nuisance colData columns to adjust
#'   for (e.g. library size, age, sex).
#' @param sigma a numeric vector of bandwidths to fit; NULL (default) uses every
#'   \code{Niche<sigma>} reducedDim present.
#' @param assay a character, the counts assay to model.
#' @param cell_type a character, the colData column of cell type labels.
#' @param sample_id a character, the colData column identifying samples
#'   (patients); used only when \code{random != "none"}.
#' @param random one of "none" (fixed-effects fit, the default), "intercept" or
#'   "slope". Adds patient-level random effects (implemented as ridge-penalised
#'   design columns) to correct anti-conservative inference caused by
#'   cell-level pseudo-replication. "intercept" adds a per-sample random
#'   intercept; "slope" additionally adds per-sample random slopes on the niche
#'   covariates (recommended, since it protects the response x niche tests). See
#'   the mixed-effects vignette.
#' @param winsor,lambda.a fitting parameters forwarded to
#'   \code{\link[SpaNorm]{fitNB}} (coefficient winsorisation and the base ridge
#'   penalty on the fixed columns).
#' @param backend a character, the fitNB compute backend
#'   ("auto", "cpu", or "gpu").
#' @param name a character, the niche reducedDim prefix.
#' @param re.maxit,re.tol iteration cap and relative tolerance for the
#'   random-effect variance-component (Schall/PQL) loop (used when
#'   \code{random != "none"}). The tolerance is a relative change on
#'   \code{log(tau2)}.
#' @param tau2.init initial random-effect variance component.
#' @param re.prop the cell-subsampling proportion used to speed up the
#'   variance-component (PQL) loop, sampled per cell type within each sample
#'   (\code{random != "none"}). For a stratum of \code{n} cells,
#'   \code{min(n, max(ceil(re.prop * n), re.min.cells))} are used. \code{1}
#'   disables subsampling (all cells, fully reproducible). The final fit that
#'   feeds inference always uses all cells; only the shared \code{tau2} estimate
#'   is affected. No seed is set internally — set one externally for
#'   reproducibility.
#' @param re.maxit.psi dispersion iterations for the inner PQL loop fits (the
#'   final all-cell fit always uses full dispersion). \code{1} (default) skips
#'   the redundant re-estimation of the barely-moving dispersion each iteration.
#' @param re.min.cells the per-stratum floor for \code{re.prop} subsampling.
#' @param BPPARAM a BiocParallelParam (reserved for the inference stage).
#' @param verbose a logical, whether to print fitting progress.
#' @param ... further arguments forwarded to \code{\link[SpaNorm]{fitNB}}.
#'
#' @return a [SpiDEResults] object (inference not yet computed).
#'
#' @examples
#' data(toySpiDE)
#' spe <- toySpiDE
#' spe <- buildNiches(spe, sigma = 20)
#' fit <- fitSpiDE(spe, condition = "condition", sigma = 20, verbose = FALSE)
#' fit
#'
#' @rdname fitSpiDE
#' @importFrom BiocParallel SerialParam
#' @export
setMethod(
  "fitSpiDE",
  signature = "ANY",
  definition = function(spe, condition, index = NULL, niche = NULL,
                        covariates = character(), sigma = NULL, assay = "counts",
                        cell_type = "cell_type", sample_id = "sample_id",
                        random = c("none", "intercept", "slope"),
                        winsor = 4, lambda.a = 0,
                        backend = c("auto", "cpu", "gpu"), name = "Niche",
                        re.maxit = 10L, re.tol = 1e-3, tau2.init = 1,
                        re.prop = 0.1, re.maxit.psi = 1L, re.min.cells = 100L,
                        BPPARAM = BiocParallel::SerialParam(), verbose = TRUE, ...) {
    backend <- match.arg(backend)
    random <- match.arg(random)
    checkSPE(spe, assay = assay, cell_type = cell_type, sample_id = sample_id)
    checkCondition(spe, condition)
    checkCovariates(spe, covariates)
    if (random != "none") {
      checkSample(spe, condition, sample_id, covariates)
      if (!is.numeric(re.prop) || length(re.prop) != 1 || re.prop <= 0 ||
          re.prop > 1) {
        stop("'re.prop' must be a single number in (0, 1]")
      }
    }

    if (is.null(sigma)) {
      sigma <- .detectSigma(spe, name)
    }
    checkNiche(spe, sigma, name = name)

    Y <- SummarizedExperiment::assay(spe, assay)
    checkCounts(Y)

    fits <- lapply(sigma, function(sg) {
      if (verbose) message(sprintf("Fitting bandwidth sigma = %s", sg))
      .fitOneBandwidth(Y, spe, condition, sg, index, niche, covariates,
                       cell_type, winsor, lambda.a, backend, name, verbose,
                       sample_id = sample_id, random = random,
                       re.maxit = re.maxit, re.tol = re.tol,
                       tau2.init = tau2.init, re.prop = re.prop,
                       re.maxit.psi = re.maxit.psi,
                       re.min.cells = re.min.cells, ...)
    })
    names(fits) <- paste0(name, sigma)

    # resolve the index / niche cell type sets actually used
    coefmap0 <- fits[[1]]@coefmap
    rn <- coefmap0$type == "ResponseNiche"
    index_used <- sort(unique(coefmap0$index[rn]))
    niche_used <- sort(unique(coefmap0$niche[rn]))

    cd <- SummarizedExperiment::colData(spe)

    new(
      "SpiDEResults",
      fits = fits,
      sigma = sigma,
      condition = condition,
      index = index_used,
      niche = niche_used,
      covariates = covariates,
      coldata = cd,
      gene.weights = NULL,
      p.cauchy.pos = NULL,
      p.cauchy.neg = NULL,
      results = data.frame(),
      fdr = NA_real_,
      call = match.call()
    )
  }
)
