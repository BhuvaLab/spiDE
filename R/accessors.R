# Accessors for SpiDEResults.

#' Extract the tidy spiDE results table
#'
#' @param object a [SpiDEResults] object.
#' @return a data.frame of significant neighbourhood-dependent DE calls, keyed
#'   by (gene, ct_index, ct_niche, bandwidth). Empty until [testSpiDE()] is run.
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = 20)
#' res <- spiDE(spe, condition = "condition", sigma = 20, verbose = FALSE)
#' head(results(res))
#' @rdname results
#' @export
setMethod("results", "SpiDEResults", function(object) object@results)

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
