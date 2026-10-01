# Generics for the spiDE user-facing API. Methods are defined in the
# corresponding implementation files (buildNiches.R, fitSpiDE.R, testSpiDE.R,
# ...). Add a new exported generic here, not inline in its implementation file.

#' @rdname buildNiches
#' @export
setGeneric("buildNiches", function(spe, ...) standardGeneric("buildNiches"))

#' @rdname mergeNiches
#' @export
setGeneric("mergeNiches", function(spe, groups, ...) standardGeneric("mergeNiches"))

#' @rdname fitSpiDE
#' @export
setGeneric("fitSpiDE", function(spe, condition = NULL, ...) standardGeneric("fitSpiDE"))

#' @rdname testSpiDE
#' @export
setGeneric("testSpiDE", function(object, condition = NULL, ...) standardGeneric("testSpiDE"))

#' @rdname spiDE
#' @export
setGeneric("spiDE", function(spe, condition = NULL, ...) standardGeneric("spiDE"))

#' @rdname results
#' @export
setGeneric("results", function(object, ...) standardGeneric("results"))

#' @rdname patientSlopes
#' @export
setGeneric("patientSlopes", function(object, ...) standardGeneric("patientSlopes"))

#' @rdname testNicheAbundance
#' @export
setGeneric("testNicheAbundance", function(spe, ...) standardGeneric("testNicheAbundance"))

#' @rdname spiGSEA
#' @export
setGeneric("spiGSEA", function(object, genesets, ...) standardGeneric("spiGSEA"))

#' @rdname patientIntercepts
#' @export
setGeneric("patientIntercepts", function(object, ...) standardGeneric("patientIntercepts"))
