# Generics for the spiDE user-facing API. Methods are defined in the
# corresponding implementation files (buildNiches.R, fitSpiDE.R, etc.).

#' @rdname buildNiches
#' @export
setGeneric("buildNiches", function(spe, ...) standardGeneric("buildNiches"))

#' @rdname mergeNiches
#' @export
setGeneric("mergeNiches", function(spe, groups, ...) standardGeneric("mergeNiches"))

#' @rdname computeSizeFactors
#' @export
setGeneric("computeSizeFactors", function(spe, ...) standardGeneric("computeSizeFactors"))

#' @rdname fitSpiDE
#' @export
setGeneric("fitSpiDE", function(spe, condition, ...) standardGeneric("fitSpiDE"))

#' @rdname testSpiDE
#' @export
setGeneric("testSpiDE", function(object, ...) standardGeneric("testSpiDE"))

#' @rdname spiDE
#' @export
setGeneric("spiDE", function(spe, condition, ...) standardGeneric("spiDE"))

#' @rdname results
#' @export
setGeneric("results", function(object, ...) standardGeneric("results"))

#' @rdname fits
#' @export
setGeneric("fits", function(object) standardGeneric("fits"))

#' @rdname bandwidths
#' @export
setGeneric("bandwidths", function(object) standardGeneric("bandwidths"))

#' @rdname spiGSEA
#' @export
# `spe` defaults to NULL in the GENERIC as well as the method: testSpiDE()
# stores the inter-gene correlation on each fit, so the counts are usually not
# needed, and a generic without the default makes the argument mandatory no
# matter what the method signature says.
setGeneric("spiGSEA", function(object, spe = NULL, genesets, ...)
  standardGeneric("spiGSEA"))