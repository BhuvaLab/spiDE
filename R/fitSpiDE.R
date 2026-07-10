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
                             verbose, ...) {
  des <- .buildNicheDesign(spe, condition, sigma, index, niche, covariates,
                           cell_type, name)
  W <- des$W

  # fit all genes at once (dispersion moderated across genes)
  fit <- SpaNorm::fitNB(Y, W, lambda.a = lambda.a, winsor = winsor,
                        backend = backend, verbose = verbose, ...)
  alpha <- fit$alpha
  rownames(alpha) <- rownames(Y)
  colnames(alpha) <- colnames(W)

  # per-gene log-likelihood for Cauchy weighting (recomputed from the fit; the
  # fitNB $loglik is per-iteration, not per-gene). Computed on the fitted-mean
  # matrix; .blockedInference() recomputes this block-wise for large data.
  mu <- SpaNorm::calculateMu(rep(0, nrow(alpha)), alpha, W)
  loglik <- rowSums(dnbinom(as.matrix(Y), mu = mu, size = 1 / fit$psi, log = TRUE))

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
    t_stat = NULL,
    se = NULL,
    p.brown.pos = NULL,
    p.brown.neg = NULL,
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
#' @param winsor,lambda.a fitting parameters forwarded to
#'   \code{\link[SpaNorm]{fitNB}} (coefficient winsorisation and ridge penalty).
#' @param backend a character, the fitNB compute backend
#'   ("auto", "cpu", or "gpu").
#' @param name a character, the niche reducedDim prefix.
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
                        cell_type = "cell_type", winsor = 4, lambda.a = 0,
                        backend = c("auto", "cpu", "gpu"), name = "Niche",
                        BPPARAM = BiocParallel::SerialParam(), verbose = TRUE, ...) {
    backend <- match.arg(backend)
    checkSPE(spe, assay = assay, cell_type = cell_type)
    checkCondition(spe, condition)
    checkCovariates(spe, covariates)

    if (is.null(sigma)) {
      sigma <- .detectSigma(spe, name)
    }
    checkNiche(spe, sigma, name = name)

    Y <- SummarizedExperiment::assay(spe, assay)
    checkCounts(Y)

    fits <- lapply(sigma, function(sg) {
      if (verbose) message(sprintf("Fitting bandwidth sigma = %s", sg))
      .fitOneBandwidth(Y, spe, condition, sg, index, niche, covariates,
                       cell_type, winsor, lambda.a, backend, name, verbose, ...)
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
