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
#' @param combine one of "cauchy" or "brown", the within-gene combiner for the
#'   correlated niche p-values (gene-level and per-index-cell-type). "cauchy"
#'   (default) is the correlation-agnostic Cauchy combination test (ACAT), which
#'   controls type-I error under correlation without estimating a correlation
#'   matrix; "brown" is the correlation-aware Brown's method. Used only if the
#'   Wald inference is not already present on the fits.
#' @param block.size a numeric, genes per inference block (NULL = a single
#'   block on the CPU backend, or a memory-bounded auto-selected size on the
#'   GPU backend).
#' @param backend a character, the compute backend for the inference stage
#'   ("auto", "cpu", or "gpu"). The GPU backend batches the per-gene Wald
#'   covariance and NB math across each block via \code{SpaNorm}'s tensor
#'   engine; it also forces a serial \code{BPPARAM} (with a warning) to avoid
#'   multiple processes contending for one GPU device. Used only if the Wald
#'   inference is not already present on the fits.
#' @param gpu.mem.budget a numeric, the GPU memory budget in bytes used to
#'   size inference blocks (NULL auto-detects; only relevant for the GPU
#'   backend). The batched Wald covariance is bounded separately, by a gene
#'   sub-batch sized from this budget on the GPU path and from
#'   \code{getOption("spiDE.cov.mem.budget", 2e9)} on the CPU path -- raise
#'   the latter to trade memory for speed on wide (random-slope) designs.
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
                        weight.thresh = 0.1, combine = c("cauchy", "brown"),
                        block.size = NULL,
                        backend = c("auto", "cpu", "gpu"),
                        gpu.mem.budget = NULL,
                        BPPARAM = BiocParallel::SerialParam(), ...) {
    checkFdr(fdr)
    combine <- match.arg(combine)
    backend <- match.arg(backend)
    if (combine == "brown" && !requireNamespace("poolr", quietly = TRUE)) {
        stop("combine = \"brown\" requires the 'poolr' package; ",
             "install it or use combine = \"cauchy\".", call. = FALSE)
    }
    fits <- object@fits

    # compute Wald + within-gene combination inference if missing
    needs_inf <- any(vapply(fits, function(f) is.null(f@t_stat), logical(1)))
    if (needs_inf) {
      if (is.null(spe)) {
        stop("inference not present on the fits; supply 'spe' to compute it")
      }
      Y <- SummarizedExperiment::assay(spe, assay)
      fits <- lapply(fits, function(f) {
        .blockedInference(f, Y, block.size = block.size, combine = combine,
                          backend = backend, gpu.mem.budget = gpu.mem.budget,
                          BPPARAM = BPPARAM)
      })
      object@fits <- fits
    }

    # combine across bandwidths
    gene.w <- .geneWeights(fits, weight.thresh)
    p.pos <- .combineBandwidths(fits, "p.combined.pos", gene.w)
    p.neg <- .combineBandwidths(fits, "p.combined.neg", gene.w)

    object@gene.weights <- gene.w
    object@p.cauchy.pos <- p.pos
    object@p.cauchy.neg <- p.neg

    # hierarchical FDR -> tidy results
    two_sided <- isTRUE(fits[[1]]@two.sided)
    object@results <- .hierarchicalFDR(fits, p.pos, p.neg, gene.w, fdr,
                                      two.sided = two_sided)

    # E2: cell-type-specific and patient-level response results. Like the
    # niche layer, both are combined ACROSS BANDWIDTHS with the same weighted
    # Cauchy combination -- every bandwidth contributes evidence, weighted by
    # that gene's relative log-likelihood. Reading them off a single fit would
    # discard most of the evidence and make the answer depend on which
    # bandwidths happened to be requested.
    #
    # These terms are niche-independent as model terms, but they are NOT
    # identical between bandwidths: each is estimated alongside that
    # bandwidth's niche columns, so combining is the right operation rather
    # than a convenience.
    .rcCols <- function(f) {
      ct <- as.character(f@covtype)
      tested <- colnames(f@t_stat)
      tt <- ct[match(tested, as.character(f@coefmap$covariate))]
      which(tt == "ResponseCellType")
    }
    rc1 <- .rcCols(fits[[1]])
    if (length(rc1)) {
      f1 <- fits[[1]]
      tested1 <- colnames(f1@t_stat)[rc1]
      cts <- f1@coefmap$index[match(tested1, f1@coefmap$covariate)]

      # per-bandwidth two-sided p on each CellType:condition column
      pcube <- lapply(fits, function(f) {
        rc <- .rcCols(f)
        dfv <- if (is.null(f@df) || length(f@df) == 1L) f@df else f@df[rc]
        tm <- f@t_stat[, rc, drop = FALSE]
        pmin(2 * pmin(.ptByCol(tm, dfv, lower.tail = FALSE),
                      .ptByCol(tm, dfv, lower.tail = TRUE)), 1)
      })
      # combine across bandwidths, column by column, with the gene weights
      p_ct <- pcube[[1]]
      for (j in seq_len(ncol(p_ct))) {
        pj <- vapply(pcube, function(P) P[, j], numeric(nrow(p_ct)))
        p_ct[, j] <- .cauchyCombine(pj, gene.w)
      }
      object@results.celltype <- .cellTypeFDRp(p_ct, cts, fdr)

      # patient level: one contrast per gene per bandwidth, combined the same way
      have_se <- vapply(fits, function(f) length(f@se_patient) > 0L, logical(1))
      if (all(have_se)) {
        pp <- vapply(fits, function(f) {
          rc <- .rcCols(f)
          tn <- colnames(f@t_stat)[rc]
          cc <- f@coefmap$index[match(tn, f@coefmap$covariate)]
          n_cells <- colSums(f@W[, paste0("CellType", cc), drop = FALSE] != 0)
          w <- n_cells / sum(n_cells)
          b <- as.numeric(f@alpha[, tn, drop = FALSE] %*% w)
          dfp <- if (is.null(f@df) || length(f@df) == 1L) f@df else
            stats::median(f@df[rc])
          tt <- b / f@se_patient[rownames(f@alpha)]
          tm <- matrix(tt, ncol = 1)
          pmin(2 * pmin(.ptByCol(tm, dfp, lower.tail = FALSE),
                        .ptByCol(tm, dfp, lower.tail = TRUE)), 1)[, 1]
        }, numeric(nrow(fits[[1]]@alpha)))
        p_pat <- .cauchyCombine(pp, gene.w)
        names(p_pat) <- rownames(fits[[1]]@alpha)
        # direction and effect size from the best-supported bandwidth
        best <- max.col(gene.w, ties.method = "first")
        bw_beta <- vapply(seq_along(fits), function(i) {
          f <- fits[[i]]; rc <- .rcCols(f)
          tn <- colnames(f@t_stat)[rc]
          cc <- f@coefmap$index[match(tn, f@coefmap$covariate)]
          n_cells <- colSums(f@W[, paste0("CellType", cc), drop = FALSE] != 0)
          as.numeric(f@alpha[, tn, drop = FALSE] %*% (n_cells / sum(n_cells)))
        }, numeric(nrow(fits[[1]]@alpha)))
        beta <- bw_beta[cbind(seq_len(nrow(bw_beta)), best)]
        names(beta) <- names(p_pat)
        object@results.patient <- .patientFDRp(p_pat, beta, fdr)
      }
    }
    object@fdr <- fdr
    object
  }
)
