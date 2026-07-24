# Convenience wrapper for the full spiDE workflow.

#' Run the full spiDE workflow
#'
#' Convenience wrapper that (optionally) builds niche covariates, fits the
#' per-gene negative binomial model over the neighbourhood-interaction design,
#' and tests for neighbourhood-dependent differential expression. Equivalent to
#' calling [buildNiches()] (if the niches are not already present),
#' [fitSpiDE()], and [testSpiDE()] in sequence.
#'
#' @inheritParams fitSpiDE
#' @param sample_id a character, the colData column identifying samples (used
#'   when niches must be built and for the random-effects fit).
#' @param backend a character, the compute backend for **both** the model fit
#'   ("auto", "cpu", or "gpu", forwarded to \code{\link[SpaNorm]{fitNB}}) and
#'   the inference stage (where the GPU backend batches the per-gene Wald
#'   covariance and negative-binomial working weights across each gene-block
#'   on the accelerator, forcing a serial \code{BPPARAM} in the process).
#' @param fdr a numeric, the target false discovery rate.
#' @param combine one of "cauchy" (default) or "brown", the within-gene combiner
#'   for the correlated niche p-values (passed to [testSpiDE()]).
#' @param block.size a numeric, genes per inference block (NULL = a single
#'   block on the CPU backend, or a memory-bounded auto-selected size on the
#'   GPU backend).
#' @param gpu.mem.budget a numeric, the GPU memory budget in bytes used to
#'   size inference blocks (NULL auto-detects; only relevant for the GPU
#'   backend).
#' @param BPPARAM a BiocParallelParam for niche construction and inference.
#'
#' @return a [SpiDEResults] object with the tidy results table populated (see
#'   [results()]).
#'
#' @examples
#' data(toySpiDE)
#' spe <- toySpiDE
#' res <- spiDE(spe, condition = "condition", sigma = c(10, 20), verbose = FALSE)
#' head(results(res))
#'
#' @rdname spiDE
#' @importFrom SingleCellExperiment reducedDimNames
#' @importFrom BiocParallel SerialParam
#' @export
setMethod(
  "spiDE",
  signature = "ANY",
  definition = function(spe, condition, index = NULL, niche = NULL,
                        covariates = character(), sigma = c(10, 30, 50, 70),
                        assay = "counts", cell_type = "cell_type",
                        sample_id = "sample_id",
                        random = c("none", "intercept", "slope"),
                        winsor = 4, lambda.a = 0,
                        backend = c("auto", "cpu", "gpu"), name = "Niche",
                        fdr = 0.05, combine = c("cauchy", "brown"),
                        block.size = NULL, gpu.mem.budget = NULL,
                        BPPARAM = BiocParallel::SerialParam(), verbose = TRUE, ...) {
    backend <- match.arg(backend)
    random <- match.arg(random)
    combine <- match.arg(combine)

    # build niches if the requested bandwidths are not already present
    need <- paste0(name, sigma)
    have <- SingleCellExperiment::reducedDimNames(spe)
    if (!all(need %in% have)) {
      if (verbose) message("Building niche covariates")
      spe <- buildNiches(spe, sigma = sigma, cell_type = cell_type,
                         sample_id = sample_id, name = name, BPPARAM = BPPARAM)
    }

    res <- fitSpiDE(spe, condition = condition, index = index, niche = niche,
                    covariates = covariates, sigma = sigma, assay = assay,
                    cell_type = cell_type, sample_id = sample_id,
                    random = random, winsor = winsor, lambda.a = lambda.a,
                    backend = backend, name = name, BPPARAM = BPPARAM,
                    verbose = verbose, ...)

    testSpiDE(res, spe = spe, assay = assay, fdr = fdr, combine = combine,
              block.size = block.size, backend = backend,
              gpu.mem.budget = gpu.mem.budget, BPPARAM = BPPARAM)
  }
)
