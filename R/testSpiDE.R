# Cross-bandwidth combination and hierarchical FDR. Given per-bandwidth fits
# with Wald/Brown inference, this combines p-values across bandwidths (Cauchy)
# and applies the hierarchical FDR to produce the tidy results table.

#' Test neighbourhood-dependent differential expression
#'
#' Completes a spiDE analysis: for each bandwidth, per-gene Wald/Brown inference
#' is computed (block-wise, in parallel) if not already present; p-values are
#' combined across bandwidths by a log-likelihood-weighted Cauchy combination;
#' and a hierarchical (gene -> index cell type -> niche cell type) Benjamini-
#' Hochberg procedure identifies significant neighbourhood-dependent DE. Results
#' are returned as a tidy table via [results()].
#'
#' @param object a [SpiDEResults] from [fitSpiDE()].
#' @param spe the SpatialExperiment used for fitting (required to compute the
#'   Wald/Brown inference if it is not already present on the fits).
#' @param assay a character, the counts assay (used only if inference must be
#'   computed).
#' @param fdr a numeric in (0, 1], the target false discovery rate (default
#'   0.05). Up/Down association directions are tested separately and combined,
#'   which is mathematically equivalent to gating the reported, combined
#'   `fdr.gene`/`fdr.index`/`fdr.niche` values at `fdr` directly - so
#'   `fdr = 1` returns every (gene, index, niche) result, unfiltered.
#' @param weight.thresh a numeric, Cauchy weights below this are set to zero.
#' @param block.size a numeric, genes per inference block (NULL = single block).
#' @param BPPARAM a BiocParallelParam for the inference stage.
#' @param ... ignored.
#'
#' @return the input \code{object} with combined p-values and the tidy results
#'   table populated.
#'
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = c(10, 20))
#' res <- fitSpiDE(spe, condition = "condition", verbose = FALSE)
#' res <- testSpiDE(res, spe = spe)
#' head(results(res))
#'
#' @rdname testSpiDE
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocParallel SerialParam
#' @export
setMethod(
  "testSpiDE",
  signature = "SpiDEResults",
  definition = function(object, spe = NULL, assay = "counts", fdr = 0.05,
                        weight.thresh = 0.1, block.size = NULL,
                        BPPARAM = BiocParallel::SerialParam(), ...) {
    checkFdr(fdr)
    fits <- object@fits

    # compute Wald/Brown inference if missing
    needs_inf <- any(vapply(fits, function(f) is.null(f@t_stat), logical(1)))
    if (needs_inf) {
      if (is.null(spe)) {
        stop("inference not present on the fits; supply 'spe' to compute it")
      }
      Y <- SummarizedExperiment::assay(spe, assay)
      fits <- lapply(fits, function(f) {
        .blockedInference(f, Y, block.size = block.size, BPPARAM = BPPARAM)
      })
      object@fits <- fits
    }

    # combine across bandwidths
    gene.w <- .geneWeights(fits, weight.thresh)
    p.pos <- .combineBandwidths(fits, "p.brown.pos", gene.w)
    p.neg <- .combineBandwidths(fits, "p.brown.neg", gene.w)

    object@gene.weights <- gene.w
    object@p.cauchy.pos <- p.pos
    object@p.cauchy.neg <- p.neg

    # hierarchical FDR -> tidy results
    object@results <- .hierarchicalFDR(fits, p.pos, p.neg, gene.w, fdr)
    object@fdr <- fdr
    object
  }
)
