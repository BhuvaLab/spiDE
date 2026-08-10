#' @name SpiDEFit
#'
#' @title An S4 class to store a single-bandwidth spiDE model fit
#'
#' @description Stores the negative binomial GLM fit and Wald/Brown inference
#'   for a single niche bandwidth. Several `SpiDEFit` objects (one per
#'   bandwidth) are combined in a [SpiDEResults] container.
#'
#' @slot sigma a numeric, the niche bandwidth (kernel standard deviation) used.
#' @slot mode a character, the design mode: "condition" (the default — a
#'   `CellType:condition:niche` design) or "niche" (a condition-free design in
#'   which the two-way `CellType:niche` interactions are the tested effects).
#' @slot ngenes a numeric, the number of genes.
#' @slot ncells a numeric, the number of cells/spots.
#' @slot W a matrix, the design matrix (cells x covariates).
#' @slot two.sided logical, TRUE when the within-gene combination used
#'   two-sided p-values (Cauchy/ACAT); FALSE for Brown's method.
#' @slot se_patient numeric, per-gene standard error of the abundance-weighted
#'   patient-level response contrast; length 0 when the design has no
#'   CellType:condition block.
#' @slot covtype a factor, the covariate type of each column of `W`, one of
#'   "CellType", "Niche", "Response", "ResponseNiche", "ResponseCellType",
#'   "Other", or "Random"
#'   (patient random-effect columns for the mixed-effects fit).
#' @slot coefmap a DataFrame mapping each covariate to its index cell type,
#'   niche cell type, and type.
#' @slot alpha a matrix, the per-gene coefficients (genes x covariates).
#' @slot gmean a numeric, the per-gene intercept (zero for the generic fit).
#' @slot psi a numeric, the per-gene negative binomial dispersion.
#' @slot rho a numeric, the average inter-gene correlation of the Pearson
#'   residuals, accumulated by \code{testSpiDE()} from the gene blocks it
#'   already loads. Used by \code{spiGSEA()} as the variance-inflation term.
#'   Empty until inference has run.
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
    mode = "character",
    ngenes = "numeric",
    ncells = "numeric",
    W = "matrix",
    covtype = "factor",
    se_patient = "numeric",
    two.sided = "logical",
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
    rho = "numeric",
    sampling = "ANY"
  ),
  prototype = list(
    re_group = NULL, tau2 = NULL, penalty = NULL, df = NULL,
    # "condition" is also what .fillSlots() gives objects serialised before
    # this slot existed -- every one of those is a condition-mode fit.
    mode = "condition",
    # Average inter-gene correlation of the Pearson residuals, accumulated by
    # .blockedInference() as a by-product of the blocks it already loads. Empty
    # until inference has run; spiGSEA() uses it as the variance-inflation term
    # so it needs no second pass over the counts.
    rho = numeric(0),
    # E2: empty unless the design carries a CellType:condition block
    se_patient = numeric(0),
    # TRUE when the within-gene combination used TWO-SIDED p-values
    # (Cauchy/ACAT). Brown's method keeps one-sided p-values, where the
    # -2log(p) transform is bounded at 0 and cannot cancel.
    two.sided = FALSE
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
  # E2 fix: "ResponseCellType" tags the CellType:condition block. The substring
  # "Response" puts it in cols_tested (so it receives SEs and t-statistics),
  # while the exact matches in inference.R and fdr.R keep it out of the
  # within-gene Cauchy combination and the triplet FDR cascade.
  valid_levels <- c("CellType", "Niche", "Response", "ResponseNiche",
                    "ResponseCellType", "Other",
                    "Random")
  if (!all(levels(object@covtype) %in% valid_levels)) {
    stop(sprintf("'covtype' levels should be a subset of: %s", paste(valid_levels, collapse = ", ")))
  }
  .checkMode(object@mode)
  TRUE
}

# Shared by both classes' validity methods.
.checkMode <- function(mode) {
  if (length(mode) != 1L || !mode %in% c("condition", "niche")) {
    stop("'mode' should be a single value, either \"condition\" or \"niche\"")
  }
  invisible(TRUE)
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
#' @slot condition a character, the condition (colData column) tested;
#'   `NA_character_` in "niche" mode.
#' @slot mode a character, the design mode: "condition" or "niche". In "niche"
#'   mode `condition` is `NA_character_` and the `results.celltype` /
#'   `results.patient` tables are empty.
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
#' @slot results.celltype a data.frame of cell-type-specific response calls
#'   keyed by (gene, ct_index), from the `CellType:condition` block. Empty
#'   unless the design carries that block. Retrieved with
#'   `results(object, type = "celltype")`.
#' @slot results.patient a data.frame of patient-level response calls, one
#'   abundance-weighted contrast per gene. Empty unless the design carries the
#'   `CellType:condition` block. Retrieved with
#'   `results(object, type = "patient")`.
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
    mode = "character",
    index = "character",
    niche = "character",
    covariates = "character",
    coldata = "ANY",
    gene.weights = "ANY",
    p.cauchy.pos = "ANY",
    p.cauchy.neg = "ANY",
    results = "data.frame",
    # E2: response results that are not niche-dependent. Empty unless the
    # design carries a CellType:condition block.
    results.celltype = "data.frame",
    results.patient = "data.frame",
    fdr = "numeric",
    call = "ANY"
  ),
  prototype = list(
    results.celltype = data.frame(),
    results.patient = data.frame(),
    mode = "condition"
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
      if (identical(.fitMode(object), "niche")) {
        "Mode: niche-only (no condition)"
      } else {
        sprintf("Condition: %s", object@condition)
      },
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
  # The fits/sigma pairing is an invariant of the GLM path, where there is one
  # SpiDEFit per bandwidth. twoStageSpiDE() has no per-gene GLM fit at all -- it
  # estimates per-patient slopes and contrasts them -- so it carries its
  # bandwidth in `sigma` with an EMPTY `fits`. Keyed on the empty list rather
  # than on a new `mode` value, because `mode` drives .testedCols() /
  # .nicheTestCols() and a third level there would silently change which columns
  # downstream code treats as tested.
  if (length(object@fits) > 0 &&
      length(object@fits) != length(object@sigma)) {
    stop("length of 'fits' does not match length of 'sigma'")
  }
  .checkMode(object@mode)
  TRUE
}

setValidity("SpiDEResults", validSpiDEResults)

# ---------------------------------------------------------------------------
# Migration of objects serialised by earlier versions.
#
# Slots added to a class do not appear in objects pickled before they existed.
# Reading such an object still works, but anything that triggers validity --
# notably initialize(), which is the documented idiom for re-combining fits
# across bandwidths -- fails with "slots in class definition but not in object".
#
# These objects can represent many hours of cluster time, so the fix is the
# Bioconductor-standard one: fill absent slots from the prototype and leave
# every present slot untouched.
# ---------------------------------------------------------------------------

#' Update a serialised spiDE object to the current class definition
#'
#' Slots added to a class do not appear in objects serialised before they
#' existed. Reading such an object still works, but anything that triggers
#' validity — notably `initialize()`, the documented idiom for re-combining
#' fits across bandwidths — fails. These methods fill any absent slot from the
#' class prototype and leave every present slot untouched.
#'
#' @param object a SpiDEFit or SpiDEResults, possibly from an older version.
#' @param ... ignored.
#' @param verbose report which slots were filled.
#' @return the object, with any slots missing since serialisation filled from
#'   the class prototype.
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = 20)
#' fit <- fitSpiDE(spe, condition = "condition", sigma = 20, verbose = FALSE)
#' updateObject(fit)
#' @importFrom BiocGenerics updateObject
#' @rdname updateObject
#' @export
setMethod("updateObject", "SpiDEFit", function(object, ..., verbose = FALSE) {
  .fillSlots(object, "SpiDEFit", verbose)
})

#' @rdname updateObject
#' @export
setMethod("updateObject", "SpiDEResults", function(object, ..., verbose = FALSE) {
  object <- .fillSlots(object, "SpiDEResults", verbose)
  if (length(object@fits)) {
    object@fits <- lapply(object@fits, function(f) {
      if (is(f, "SpiDEFit")) .fillSlots(f, "SpiDEFit", verbose) else f
    })
  }
  object
})

#' Fill slots absent from a serialised object using the class prototype
#' @noRd
.fillSlots <- function(object, Class, verbose = FALSE) {
  proto <- methods::new(Class)
  missing <- setdiff(methods::slotNames(Class), names(attributes(object)))
  for (s in missing) attr(object, s) <- methods::slot(proto, s)
  if (verbose && length(missing)) {
    message("updateObject(", Class, "): filled ", paste(missing, collapse = ", "))
  }
  object
}
