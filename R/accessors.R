# Accessors for SpiDEResults.

#' Extract the tidy spiDE results table
#'
#' @param object a [SpiDEResults] object.
#' @param type which response result to return. "niche" (default) gives the
#'   neighbourhood-dependent calls keyed by (gene, ct_index, ct_niche,
#'   bandwidth), gated by the 3-level gene -> index -> niche cascade.
#'   "celltype" gives cell-type-specific response calls keyed by
#'   (gene, ct_index), gated by a 2-level gene -> cell type cascade.
#'   "patient" gives one abundance-weighted response contrast per gene.
#'   The latter two are empty unless the design carries a CellType:condition
#'   block.
#' @return a data.frame of significant calls; see \code{type}. Empty until
#'   [testSpiDE()] is run.
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = 20)
#' res <- spiDE(spe, condition = "condition", sigma = 20, verbose = FALSE)
#' head(results(res))
#' @rdname results
#' @export
setMethod("results", "SpiDEResults",
          function(object, type = c("niche", "celltype", "patient"), ...) {
  type <- match.arg(type)
  switch(type,
    niche    = object@results,
    celltype = object@results.celltype,
    patient  = object@results.patient)
})

#' Extract per-bandwidth fits
#'
#' @param object a [SpiDEResults] object.
#' @return a named list of [SpiDEFit] objects, one per bandwidth.
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = 20)
#' res <- fitSpiDE(spe, condition = "condition", sigma = 20, verbose = FALSE)
#' fits(res)
#' @rdname fits
#' @export
setMethod("fits", "SpiDEResults", function(object) object@fits)

#' Bandwidths used in a spiDE analysis
#'
#' @param object a [SpiDEResults] object.
#' @return a numeric vector of niche bandwidths.
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = c(10, 20))
#' res <- fitSpiDE(spe, condition = "condition", verbose = FALSE)
#' bandwidths(res)
#' @rdname bandwidths
#' @export
setMethod("bandwidths", "SpiDEResults", function(object) object@sigma)
