#' Run spiDE end to end
#'
#' Builds the niche covariates at bandwidth \code{sigma} (unless \code{spe}
#' already carries them), fits the per-patient niche models ([fitSpiDE()]) and
#' tests them ([testSpiDE()]): the pooled niche test always, and the
#' condition-specific test when a \code{condition} is given.
#'
#' @param spe a SpatialExperiment with raw counts.
#' @param condition \code{NULL} or the colData column of a two-level,
#'   patient-level condition.
#' @param sigma the niche bandwidth.
#' @param engine \code{"slopes"} or \code{"sandwich"}, see [fitSpiDE()].
#' @param procedure,fdr passed to [testSpiDE()].
#' @param cell_type,sample_id the colData columns of cell type and patient.
#' @param BPPARAM a BiocParallelParam.
#' @param verbose report progress.
#' @param ... further arguments to [fitSpiDE()].
#' @return a [SpiDEResults-class].
#' @examples
#' data(toySpiDE)
#' res <- spiDE(toySpiDE, condition = "condition", sigma = 20, index = "A",
#'              min.patients = 6, procedure = "all")
#' head(results(res))
#' @rdname spiDE
#' @export
setMethod(
  "spiDE", "SpatialExperiment",
  function(spe, condition = NULL, sigma, engine = c("slopes", "sandwich"),
           procedure = c("filtered", "all"), fdr = 0.05, cell_type = "cell_type",
           sample_id = "sample_id", BPPARAM = BiocParallel::SerialParam(), verbose = TRUE, ...) {
    engine <- match.arg(engine)
    procedure <- match.arg(procedure)
    if (missing(sigma) || length(sigma) != 1L) stop("'sigma' must be a single bandwidth", call. = FALSE)
    if (!paste0("Niche", sigma) %in% SingleCellExperiment::reducedDimNames(spe)) {
      spe <- buildNiches(spe, sigma = sigma, cell_type = cell_type, sample_id = sample_id,
                         BPPARAM = BPPARAM, verbose = verbose)
    }
    fit <- fitSpiDE(spe, condition = condition, engine = engine, sigma = sigma,
                    cell_type = cell_type, sample_id = sample_id, BPPARAM = BPPARAM,
                    verbose = verbose, ...)
    testSpiDE(fit, condition = condition, procedure = procedure, fdr = fdr)
  }
)
