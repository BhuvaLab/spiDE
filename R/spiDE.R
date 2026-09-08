# Convenience wrapper for the full spiDE workflow.

#' Run the full spiDE workflow
#'
#' Convenience wrapper that (optionally) builds niche covariates, fits the
#' per-gene negative binomial model over the neighbourhood-interaction design,
#' and tests for neighbourhood-dependent differential expression. Equivalent to
#' calling [buildNiches()] (if the niches are not already present),
#' [fitSpiDE()], and [testSpiDE()] in sequence.
#'
#' @param condition a character, the colData column of the tested condition, or
#'   \code{NULL} for a condition-free (niche-only) analysis. See [fitSpiDE()].
#' @inheritParams fitSpiDE
#' @param sample_id a character, the colData column identifying samples (used
#'   when niches must be built and for the random-effects fit).
#' @param backend a character, the compute backend for **both** the model fit
#'   ("auto", "cpu", or "gpu", forwarded to \code{\link[SpaNorm]{fitNB}}) and
#'   the inference stage (where the GPU backend batches the per-gene Wald
#'   covariance and negative-binomial working weights across each gene-block
#'   on the accelerator, forcing a serial \code{BPPARAM} in the process).
#' @param re.celltype logical; add a nested (sample x cell type) random
#'   intercept so the tested niche slopes are within-group. Default
#'   \code{TRUE}. See [fitSpiDE()].
#' @param polish logical; run the polish stage ([polishSpiDE()]) between the
#'   fit and the test: converge each gene, set its dispersion and re-estimate
#'   the variance components from the converged fit. Default \code{TRUE};
#'   \code{FALSE} tests the shared fit as \code{fitNB} returned it. Run the
#'   stages by hand to control the polish's own settings.
#' @param fdr a numeric, the target false discovery rate.
#' @param combine one of "cauchy" (default) or "brown", the within-gene combiner
#'   for the correlated niche p-values (passed to [testSpiDE()]).
#' @param dispersion the standard-error scale; see [testSpiDE()].
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
#' res <- spiDE(spe, condition = "condition", sigma = 20, verbose = FALSE)
#' head(results(res))
#'
#' @rdname spiDE
#' @importFrom SingleCellExperiment reducedDimNames
#' @importFrom BiocParallel SerialParam
#' @export
setMethod(
  "spiDE",
  signature = "ANY",
  definition = function(spe, condition = NULL, index = NULL, niche = NULL,
                        covariates = character(), sigma = c(10, 30, 50, 70),
                        assay = "counts", cell_type = "cell_type",
                        sample_id = "sample_id",
                        random = c("intercept", "none", "slope"),
                        winsor = 4, lambda.a = 0,
                        backend = c("auto", "cpu", "gpu"), name = "Niche",
                        fdr = 0.05, combine = c("cauchy", "brown"),
                        df.method = c("satterthwaite", "between"),
                        re.celltype = TRUE, polish = TRUE,
                        dispersion = c("ql", "pearson"),
                        block.size = NULL, gpu.mem.budget = NULL,
                        BPPARAM = BiocParallel::SerialParam(), verbose = TRUE, ...) {
    backend <- match.arg(backend)
    random <- match.arg(random)
    combine <- match.arg(combine)
    dispersion <- match.arg(dispersion)
    df.method <- match.arg(df.method)

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
                    backend = backend, name = name, df.method = df.method,
                    re.celltype = re.celltype, block.size = block.size,
                    BPPARAM = BPPARAM, verbose = verbose, ...)
    if (polish) {
      res <- polishSpiDE(res, spe, assay = assay, block.size = block.size,
                         BPPARAM = BPPARAM, verbose = verbose)
    }

    testSpiDE(res, spe = spe, assay = assay, fdr = fdr, combine = combine,
              block.size = block.size, backend = backend,
              gpu.mem.budget = gpu.mem.budget, BPPARAM = BPPARAM,
                     dispersion = dispersion)
  }
)
