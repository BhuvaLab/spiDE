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
                        re.maxit = 10L, re.tol = 1e-3, tau2.init = 1,
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
#' group entries use B_{a,b} = A[ca,cb] - A[ca,] M^{-1} A[,cb] = Z_a' P Z_b.
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
                             re.prop = 1, re.maxit.psi = 1L,
                             re.min.cells = 100L, df.method = "satterthwaite", ...) {
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
                       re.maxit.psi = re.maxit.psi,
                       df.method = df.method,
                       cols_tested = .testedCols(des$covtype, des$mode),
                       mode = des$mode,
                       ...)
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
    mode = des$mode,
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
#'   (must have exactly two levels), or \code{NULL} for a condition-free
#'   (niche-only) analysis. With \code{NULL} the condition terms are dropped
#'   from the design and the two-way \code{CellType:niche} interactions become
#'   the tested effects: within index cell type \emph{c}, how expression
#'   changes with the local density of niche cell type \emph{n}. The
#'   \code{results(type = "celltype")} and \code{results(type = "patient")}
#'   tables are empty in that mode, there being no condition to contrast.
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
#'   covariates, which protects the response x niche tests when the niche-slope
#'   varies between samples. Note that "slope" estimates an extra variance
#'   component (\code{tau2} for \code{SampleSlope}) that is collinear with the
#'   tested fixed effect and, with few samples, is less stable than the
#'   intercept-only variance; prefer "intercept" at small \code{S}, and "slope"
#'   when between-sample niche-slope variation is expected. The fit is also
#'   stochastic (\code{fitNB} subsamples cells for the dispersion estimate), so
#'   set a seed for reproducible variance components. See the mixed-effects and
#'   simulation vignettes.
#'
#'   In a condition-free analysis (\code{condition = NULL}) the random-slope
#'   block sits on exactly the \code{CellType:niche} columns being tested, so
#'   \code{"slope"} is the natural correction there when between-sample
#'   variation in niche slopes is plausible. Be aware of a limitation specific
#'   to that mode: the tested slope is a \emph{within-sample} contrast on a
#'   \emph{spatially autocorrelated} covariate, and spiDE does not model
#'   spatial autocorrelation. Neighbouring cells are therefore not independent
#'   replicates of the slope, and neither random-effect structure can recover
#'   that -- on a null fixture with a planted per-sample intercept, the
#'   fixed-effects fit made 37 calls, \code{"intercept"} 5 and \code{"slope"}
#'   5 (see \code{longtests/testthat/test-nicheOnly-mixed.R}). Random effects
#'   remove most of the inflation but niche mode remains mildly
#'   anti-conservative; treat borderline calls with corresponding caution.
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
#'   (the default) disables subsampling (all cells, fully reproducible). The
#'   final fit that feeds inference always uses all cells; only the shared
#'   \code{tau2} estimate is affected. No seed is set internally — set one
#'   externally for reproducibility of \code{re.prop < 1} runs. Lowering
#'   \code{re.prop} trades accuracy for speed, and the trade is worse than it
#'   looks: a replicate study on real data (\code{vignettes/spiDE-mixed-
#'   benchmark.Rmd}) found that subsampling noise in \code{tau2} does not
#'   shrink as \code{re.prop} rises from 0.2 to 0.8 (it stays comparable to or
#'   larger than genuine between-patient variation), and — more importantly —
#'   \code{tau2} is systematically \emph{biased downward} at every
#'   \code{re.prop < 1} tested, an attenuation that averaging replicates
#'   cannot fix, only shrinking as \code{re.prop} approaches 1. Using
#'   \code{re.prop < 1} is therefore rarely advised; only do so when the
#'   variance component's absolute scale doesn't matter (e.g. a quick
#'   feasibility check) and treat its \code{tau2} as a lower bound, not a
#'   point estimate.
#' @param re.maxit.psi dispersion iterations for the inner PQL loop fits (the
#'   final all-cell fit always uses full dispersion). \code{1} (default) skips
#'   the redundant re-estimation of the barely-moving dispersion each iteration.
#' @param re.min.cells the per-stratum floor for \code{re.prop} subsampling.
#' @param df.method one of "satterthwaite" (default) or "between"; only used
#'   when \code{random != "none"}. "satterthwaite" derives a separate df per
#'   tested column from the shared variance-component fit, distinguishing
#'   between-sample contrasts (Response: small df, close to "between") from
#'   within-sample contrasts (ResponseNiche: larger df, more power) rather than
#'   applying \code{S - 2} to both; \code{@df} is then a named per-column
#'   vector. "between" tests every Response/ResponseNiche coefficient against
#'   the same scalar between-sample reference df (\code{S - 2}), the original
#'   back-compatible behaviour, and \code{@df} is a scalar.
#'
#'   The default changed to "satterthwaite" after the benchmark study
#'   (\code{research/}) measured both arms on identically seeded data: "between"
#'   is severely over-conservative when samples are few (null type-I
#'   \eqn{\approx 0.001} at \eqn{S = 4} against a nominal 0.05, with
#'   correspondingly near-zero power), while "satterthwaite" holds type-I in
#'   \eqn{0.042}-\eqn{0.065} across the whole sampled range and gains
#'   \eqn{\approx 0.10} mean TPR. The trade is a mild liberal drift at larger
#'   \eqn{S} (worst measured \eqn{\approx 0.065}); use "between" when strict
#'   conservatism matters more than power, or for back-compatibility.
#'   Ignored when \code{random == "none"}.
#'
#'   In a condition-free analysis (\code{condition = NULL}) the
#'   \code{"between"} reference df changes, because the tested
#'   \code{CellType:niche} slope is a within-sample contrast rather than a
#'   between-condition one: it is \code{ncells - p_fixed} under
#'   \code{random = "intercept"} (cells are the replicates) and \code{S - 1}
#'   under \code{random = "slope"} (the per-sample random slopes sit on the
#'   tested columns, moving the contrast into the between-sample stratum). The
#'   name \code{"between"} is therefore a misnomer in the intercept case; it is
#'   retained for back-compatibility. \code{"satterthwaite"} computes this
#'   distinction from the fitted variance components and is preferred.
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
#' fit0 <- fitSpiDE(spe, condition = NULL, sigma = 20, verbose = FALSE)
#' fit0
#'
#' @rdname fitSpiDE
#' @importFrom BiocParallel SerialParam
#' @export
setMethod(
  "fitSpiDE",
  signature = "ANY",
  definition = function(spe, condition = NULL, index = NULL, niche = NULL,
                        covariates = character(), sigma = NULL, assay = "counts",
                        cell_type = "cell_type", sample_id = "sample_id",
                        random = c("none", "intercept", "slope"),
                        winsor = 4, lambda.a = 0,
                        backend = c("auto", "cpu", "gpu"), name = "Niche",
                        re.maxit = 10L, re.tol = 1e-3, tau2.init = 1,
                        re.prop = 1, re.maxit.psi = 1L, re.min.cells = 100L,
                        df.method = c("satterthwaite", "between"),
                        BPPARAM = BiocParallel::SerialParam(), verbose = TRUE, ...) {
    backend <- match.arg(backend)
    random <- match.arg(random)
    df.method <- match.arg(df.method)
    checkSPE(spe, assay = assay, cell_type = cell_type, sample_id = sample_id)
    # condition = NULL selects the condition-free (niche-only) design
    if (!is.null(condition)) checkCondition(spe, condition)
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
                       re.min.cells = re.min.cells, df.method = df.method, ...)
    })
    names(fits) <- paste0(name, sigma)

    mode <- if (is.null(condition)) "niche" else "condition"

    # Resolve the index / niche cell type sets actually used. This must go
    # through the mode predicate like every other tested-column lookup: a bare
    # == "ResponseNiche" is all-FALSE in niche mode, which would leave @index
    # and @niche empty on every condition-free result.
    coefmap0 <- fits[[1]]@coefmap
    rn <- .nicheTestCols(coefmap0$type, mode)
    index_used <- sort(unique(coefmap0$index[rn]))
    niche_used <- sort(unique(coefmap0$niche[rn]))

    cd <- SummarizedExperiment::colData(spe)

    new(
      "SpiDEResults",
      fits = fits,
      sigma = sigma,
      condition = if (is.null(condition)) NA_character_ else condition,
      mode = mode,
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
