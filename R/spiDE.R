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
#' @param procedure,fdr,pooled.df passed to [testSpiDE()]; \code{procedure}
#'   defaults to the engine's default there.
#' @param cell_type,sample_id the colData columns of cell type and patient.
#' @param BPPARAM a BiocParallelParam.
#' @param verbose report progress.
#' @param ... further arguments to [fitSpiDE()].
#' @return a [SpiDEResults-class].
#' @examples
#' data(toySpiDE)
#' res <- spiDE(toySpiDE, condition = "condition", sigma = 30, index = "A",
#'              procedure = "all")
#' head(results(res))
#' @rdname spiDE
#' @export
setMethod(
  "spiDE", "SpatialExperiment",
  function(spe, condition = NULL, sigma, engine = c("slopes", "sandwich"),
           procedure = c("heterogeneity", "filtered", "all"), fdr = 0.05,
           pooled.df = c("proportional", "capped"), cell_type = "cell_type",
           sample_id = "sample_id", BPPARAM = BiocParallel::SerialParam(), verbose = TRUE, ...) {
    engine <- match.arg(engine)
    procedure <- if (missing(procedure)) NULL else match.arg(procedure)
    pooled.df <- match.arg(pooled.df)
    # refuse before the fit, not after it
    if (identical(procedure, "heterogeneity") && engine != "slopes") {
      stop("procedure = \"heterogeneity\" needs the slopes engine's per-patient slopes; ",
           "use procedure = \"filtered\" or \"all\" with the sandwich engine", call. = FALSE)
    }
    if (missing(sigma) || length(sigma) != 1L) stop("'sigma' must be a single bandwidth", call. = FALSE)
    if (!paste0("Niche", sigma) %in% SingleCellExperiment::reducedDimNames(spe)) {
      spe <- buildNiches(spe, sigma = sigma, cell_type = cell_type, sample_id = sample_id,
                         BPPARAM = BPPARAM, verbose = verbose)
    }
    fit <- fitSpiDE(spe, condition = condition, engine = engine, sigma = sigma,
                    cell_type = cell_type, sample_id = sample_id, BPPARAM = BPPARAM,
                    verbose = verbose, ...)
    testSpiDE(fit, condition = condition, procedure = procedure, fdr = fdr, pooled.df = pooled.df,
              BPPARAM = BPPARAM)
  }
)
