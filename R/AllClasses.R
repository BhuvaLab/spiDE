#' @name SpiDEFit
#'
#' @title An S4 class to store a single-bandwidth spiDE model fit
#'
#' @description Stores the negative binomial GLM fit and Wald/Brown inference
#'   for a single niche bandwidth. Several `SpiDEFit` objects (one per
#'   bandwidth) are combined in a [SpiDEResults] container.
#'
#' @slot sigma a numeric, the niche bandwidth (kernel standard deviation) used.
#' @slot ngenes a numeric, the number of genes.
#' @slot ncells a numeric, the number of cells/spots.
#' @slot W a matrix, the design matrix (cells x covariates).
#' @slot covtype a factor, the covariate type of each column of `W`, one of
#'   "CellType", "Niche", "Response", "ResponseNiche", "Other", or "Random"
#'   (patient random-effect columns for the mixed-effects fit).
#' @slot coefmap a DataFrame mapping each covariate to its index cell type,
#'   niche cell type, and type.
#' @slot alpha a matrix, the per-gene coefficients (genes x covariates).
#' @slot gmean a numeric, the per-gene intercept (zero for the generic fit).
#' @slot psi a numeric, the per-gene negative binomial dispersion.
#' @slot loglik a numeric, the per-gene log-likelihood (used for Cauchy weights).
#' @slot re_group a character (or NULL), the random-effect group of each column
#'   of `W` (`NA` for fixed columns); NULL for a fixed-effects fit.
#' @slot tau2 a numeric (or NULL), the fitted random-effect variance components
#'   (one per random-effect group); NULL for a fixed-effects fit.
#' @slot penalty a numeric (or NULL), the per-column ridge penalty (`lambda.a`)
#'   used at fit time (0 on fixed columns, `1/tau2` on random columns).
#' @slot df a numeric (or NULL), the Wald reference degrees of freedom. NULL for
#'   a fixed-effects fit (normal reference); a scalar between-patient `S - 2`
#'   under `df.method = "between"`; or a named per-tested-coefficient vector
#'   (aligned to the columns of `t_stat`/`se`) under `df.method = "satterthwaite"`
#'   — the Response main effect stays ~ `S - 2` while niche interactions, informed
#'   by within-sample variation, get a larger df.
#' @slot t_stat a matrix, per-gene Wald t-statistics (genes x covariates).
#' @slot se a matrix, per-gene coefficient standard errors (genes x covariates).
#' @slot p.combined.pos a matrix, combined p-values for up-regulation
#'   (genes x (gene-level plus per-index-cell-type)). The within-gene combiner
#'   is Brown's method or the Cauchy combination test (see \code{combine}).
#' @slot p.combined.neg a matrix, combined p-values for down-regulation.
#' @slot sampling a factor, cells used for GLM/dispersion estimation (from
#'   [SpaNorm::fitNB()]).
#'
#' @param x an object of class SpiDEFit.
#' @param name a character, the slot to retrieve.
#' @return Return value varies depending on method.
NULL

#' @rdname SpiDEFit
#' @export
#' @import methods
setClass(
  Class = "SpiDEFit",
  slots = c(
    sigma = "numeric",
    ngenes = "numeric",
    ncells = "numeric",
    W = "matrix",
    covtype = "factor",
    coefmap = "ANY",
    alpha = "matrix",
    gmean = "numeric",
    psi = "numeric",
    loglik = "numeric",
    re_group = "ANY",
    tau2 = "ANY",
    penalty = "ANY",
    df = "ANY",
    t_stat = "ANY",
    se = "ANY",
    p.combined.pos = "ANY",
    p.combined.neg = "ANY",
    sampling = "ANY"
  ),
  prototype = list(
    re_group = NULL, tau2 = NULL, penalty = NULL, df = NULL
  )
)

#' @rdname SpiDEFit
#' @export
setMethod(
  f = "$",
  signature = "SpiDEFit",
  definition = function(x, name) {
    return(slot(x, name))
  }
)

setMethod(
  f = "show",
  signature = "SpiDEFit",
  definition = function(object) {
    cat(
      is(object)[[1]],
      sprintf("Bandwidth (sigma): %s", object@sigma),
      sprintf("Data: %d genes, %d cells/spots", object@ngenes, object@ncells),
      sprintf("Covariates: %d (%s)", ncol(object@W), paste(sprintf("%s=%d", levels(object@covtype), table(object@covtype)), collapse = ", ")),
      sprintf("alpha: %s", utils::capture.output(utils::str(object@alpha))),
      sprintf("psi: %s", utils::capture.output(utils::str(object@psi))),
      sprintf("Inference computed: %s", !is.null(object@t_stat)),
      sep = "\n"
    )
  }
)

validSpiDEFit <- function(object) {
  # dimension checks
  if (ncol(object@alpha) != ncol(object@W)) {
    stop("ncol of 'alpha' does not match ncol of 'W'")
  }
  if (length(object@covtype) != ncol(object@W)) {
    stop("length of 'covtype' does not match ncol of 'W'")
  }
  if (nrow(object@alpha) != object@ngenes) {
    stop("nrow of 'alpha' does not match 'ngenes'")
  }
  if (nrow(object@W) != object@ncells) {
    stop("nrow of 'W' does not match 'ncells'")
  }
  if (length(object@psi) != object@ngenes) {
    stop("length of 'psi' does not match 'ngenes'")
  }
  # NA checks
  if (any(is.na(object@W))) {
    stop("'W' cannot have missing values")
  }
  if (any(is.na(object@alpha))) {
    stop("'alpha' cannot have missing values")
  }
  # covtype levels
  valid_levels <- c("CellType", "Niche", "Response", "ResponseNiche", "Other",
                    "Random")
  if (!all(levels(object@covtype) %in% valid_levels)) {
    stop(sprintf("'covtype' levels should be a subset of: %s", paste(valid_levels, collapse = ", ")))
  }
  TRUE
}

setValidity("SpiDEFit", validSpiDEFit)


#' @name SpiDEResults
#'
#' @title An S4 class to store multi-bandwidth spiDE results
#'
#' @description Container for one or more [SpiDEFit] objects (one per niche
#'   bandwidth) plus the cross-bandwidth combined inference and the tidy
#'   results table produced by [testSpiDE()].
#'
#' @slot fits a list of [SpiDEFit] objects, one per bandwidth.
#' @slot sigma a numeric, the bandwidth grid.
#' @slot condition a character, the condition (colData column) tested.
#' @slot index a character, the index cell types considered.
#' @slot niche a character, the niche cell types considered.
#' @slot covariates a character, the nuisance covariates included.
#' @slot coldata a DataFrame of sample-level metadata.
#' @slot gene.weights a matrix of per-gene Cauchy-combination weights
#'   (genes x bandwidth).
#' @slot p.cauchy.pos a matrix, Cauchy-combined up-regulation p-values.
#' @slot p.cauchy.neg a matrix, Cauchy-combined down-regulation p-values.
#' @slot results a data.frame, the tidy results table (empty until
#'   [testSpiDE()] is run), keyed by (gene, ct_index, ct_niche, bandwidth).
#' @slot fdr a numeric, the FDR threshold used.
#' @slot call the matched call that produced the object.
#'
#' @param x an object of class SpiDEResults.
#' @param name a character, the slot to retrieve.
#' @return Return value varies depending on method.
NULL

#' @rdname SpiDEResults
#' @export
setClass(
  Class = "SpiDEResults",
  slots = c(
    fits = "list",
    sigma = "numeric",
    condition = "character",
    index = "character",
    niche = "character",
    covariates = "character",
    coldata = "ANY",
    gene.weights = "ANY",
    p.cauchy.pos = "ANY",
    p.cauchy.neg = "ANY",
    results = "data.frame",
    fdr = "numeric",
    call = "ANY"
  )
)

#' @rdname SpiDEResults
#' @export
setMethod(
  f = "$",
  signature = "SpiDEResults",
  definition = function(x, name) {
    return(slot(x, name))
  }
)

setMethod(
  f = "show",
  signature = "SpiDEResults",
  definition = function(object) {
    ng <- if (length(object@fits) > 0) object@fits[[1]]@ngenes else NA_integer_
    cat(
      is(object)[[1]],
      sprintf("Bandwidths (sigma): %s", paste(object@sigma, collapse = ", ")),
      sprintf("Genes: %s", ng),
      sprintf("Condition: %s", object@condition),
      sprintf("Index cell types: %s", paste(object@index, collapse = ", ")),
      sprintf("Niche cell types: %s", paste(object@niche, collapse = ", ")),
      sprintf("Tested: %s (%d rows in results table)", nrow(object@results) > 0, nrow(object@results)),
      sep = "\n"
    )
  }
)

validSpiDEResults <- function(object) {
  if (length(object@fits) > 0 && !all(vapply(object@fits, is, logical(1), "SpiDEFit"))) {
    stop("'fits' should be a list of SpiDEFit objects")
  }
  if (length(object@fits) != length(object@sigma)) {
    stop("length of 'fits' does not match length of 'sigma'")
  }
  TRUE
}

setValidity("SpiDEResults", validSpiDEResults)
