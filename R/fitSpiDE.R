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

#' Fit the spiDE negative binomial GLM for one bandwidth
#'
#' @return a SpiDEFit with the fit populated and inference slots empty.
#' @importFrom SummarizedExperiment assay
#' @importFrom stats dnbinom
#' @noRd
.fitOneBandwidth <- function(Y, spe, condition, sigma, index, niche, covariates,
                             cell_type, winsor, lambda.a, backend, name,
                             verbose, sample_id = "sample_id", random = "none",
                             re.maxit = 2L, re.tol = 1e-3, tau2.init = 1,
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
    # Inf, not NULL: a t reference with infinite df IS the normal reference
    # (pt(x, Inf) == pnorm(x)), so @df is always a populated numeric and no
    # consumer needs an is.null() branch.
    df <- Inf
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
#'
#'   The default was lowered from 10 to 2 on measurement: for
#'   \code{random = "intercept"} one iteration is indistinguishable from ten on
#'   null type-I error (to three decimal places) and on \code{tau2} (to two),
#'   because the loop typically converges in a couple of steps. It also largely
#'   dissolves a hazard of the larger cap -- \code{tau2} can enter a 2-cycle,
#'   making the result depend on the parity of \code{re.maxit}.
#'
#'   \strong{This evidence covers the intercept model only.} Under
#'   \code{random = "slope"} the slope variance component decays monotonically
#'   across all ten iterations without meeting \code{re.tol}, so a cap of 2
#'   leaves it far from where 10 leaves it. Pass \code{re.maxit = 10} for
#'   slope fits.
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
                        re.maxit = 2L, re.tol = 1e-3, tau2.init = 1,
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
