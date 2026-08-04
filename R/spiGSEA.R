# Gene-set inference on a fitted spiDE model.
#
# The per-gene niche tests are individually under-powered on sparse spatial
# data: a single gene's response-niche coefficient is estimated from few
# effective replicates, so the hit list skews toward whatever is most highly
# expressed. Averaging a t-statistic over the tens-to-hundreds of genes in a
# pathway recovers power, because the gene-level noise partly averages out
# while a coherent pathway-level shift does not.
#
# Both available tests share a variance inflation factor: for a set of m genes
# with average inter-gene correlation rho, the mean of m unit-variance
# statistics has
#
#     sd = sqrt((1 + rho (m - 1)) / m)
#
# rather than sqrt(1/m). This term is not optional. With rho = 0 a co-regulated
# pathway looks enormously significant purely because its genes are not
# independent -- the effective sample size is a fraction of m.
#
# What the two tests differ on is the null they compare against, and they answer
# genuinely different questions:
#
#   "self-contained" (default) tests the set's mean z against ZERO. Null: no
#       gene in the set responds. This is the test the flat-script pipeline this
#       replaces used, and it is what its `fry_res` object holds.
#   "competitive" tests the set's mean z against the mean of the genes OUTSIDE
#       the set, as limma::camera does. Null: genes in the set respond no more
#       than the rest of the assayed genes.
#
# The distinction matters here more than usual. spiDE's per-gene statistics are
# known to track expression, so a set of abundant genes can carry a non-zero
# mean shift without carrying niche-specific signal -- which the self-contained
# test will call and the competitive test largely will not, because the
# background it subtracts is subject to the same bias.
#
# Two departures from the flat-script version this replaces:
#
#   1. Statistics are converted to z BEFORE averaging. The set statistic assumes
#      unit variance, but a t with v df has variance v/(v-2). For the niche
#      coefficients v is large and the distinction is cosmetic; for the
#      CellType:condition coefficients v is the between-patient S-2 (single
#      digits), where t has variance ~1.3 and treating it as unit-variance makes
#      the test anti-conservative by ~15%.
#   2. Bandwidths are combined on TWO-SIDED p-values. Under one-sided input the
#      Cauchy transform tan((0.5 - p)pi) diverges to -Inf as p -> 1, so a set
#      that shifts up at one bandwidth and down at another cancels exactly
#      instead of registering as heterogeneous. This is the same correction
#      already applied to the gene-level combination (see R/combine.R).

#' Signed z-scores from t-statistics
#'
#' The tail is evaluated on the log scale and inverted there, because
#' \code{pt()} underflows to 0 for |t| beyond about 9 at large df -- which
#' \code{qnorm()} would then map to \code{-Inf}, silently turning the most
#' significant genes in a set into missing values.
#'
#' @param tmat a genes x k matrix of t-statistics.
#' @param df NULL (normal reference), a scalar, or a per-column vector.
#' @return a matrix of z-scores with the same shape and signs as \code{tmat}.
#' @importFrom stats qnorm
#' @noRd
.tToZ <- function(tmat, df) {
  if (is.null(df)) {
    return(tmat)
  }
  lp <- .ptByCol(-abs(tmat), df, lower.tail = TRUE)
  # .ptByCol returns probabilities, not logs; guard the underflow explicitly
  lp[lp < .Machine$double.xmin] <- .Machine$double.xmin
  z <- -stats::qnorm(lp)
  z[!is.finite(z)] <- 0
  sign(tmat) * z
}

#' Mean off-diagonal correlation from accumulated per-cell z sums
#'
#' The identity the streaming estimator rests on, isolated so it can be checked
#' against \code{cor()} directly rather than only through a fitted model.
#'
#' @param s per-cell sum, over genes, of standardised residuals.
#' @param n number of cells.
#' @param p number of genes contributing to \code{s}.
#' @noRd
.meanCorFromZ <- function(s, n, p) {
  if (p < 2L) {
    return(0)
  }
  (sum(s^2) / (n - 1) - p) / (p * (p - 1))
}

#' Mean pairwise correlation of Pearson residuals, without forming the matrix
#'
#' The average inter-gene correlation is a single scalar, but the obvious route
#' to it -- correlate every gene against every other -- is a p x p matrix. At
#' 13,000 genes that is 1.4 GB purely to take one mean of it. Standardise each
#' gene's residuals and note that \code{sum(cor)} over ALL entries equals
#' \code{sum_cells (sum_genes z)^2 / (n - 1)}, so the per-cell sum over genes is
#' accumulable one block at a time: exact, in O(cells) memory, one pass.
#'
#' The estimator streams the counts one gene block at a time, which is the only
#' heavy step in the gene-set layer -- everything downstream operates on a
#' sets x contrasts matrix that is orders of magnitude smaller. It therefore
#' carries the same backend and parallel controls as \code{.blockedInference()}:
#' blocks are independent (each contributes an additive per-cell partial sum),
#' so they distribute over CPU workers exactly, and the per-block arithmetic is
#' dense elementwise work on a block x cells matrix, which is what the
#' accelerator is for.
#'
#' @param Y counts (genes x cells; dense, sparse or DelayedArray).
#' @param fit a SpiDEFit.
#' @param genes integer row indices to use.
#' @param block.size genes per block; auto-sized on the GPU backend.
#' @param backend "auto", "cpu" or "gpu".
#' @param gpu.mem.budget device memory budget, passed to the block sizer.
#' @param BPPARAM a BiocParallelParam for the CPU block loop.
#' @return a single numeric, the mean pairwise correlation.
#' @noRd
.interGeneCor <- function(Y, fit, genes, block.size = 500L,
                          backend = c("auto", "cpu", "gpu"),
                          gpu.mem.budget = NULL,
                          BPPARAM = BiocParallel::SerialParam()) {
  backend <- match.arg(backend)
  n <- ncol(Y)

  # Resolve GPU state once in the parent, before any dispatch: the short-circuit
  # keeps a forked worker from ever probing the torch runtime.
  gpu_active <- backend %in% c("gpu", "auto") && SpaNorm::checkGPU()
  if (gpu_active && BiocParallel::bpnworkers(BPPARAM) > 1) {
    warning("GPU backend active (", SpaNorm::getBackendDevice(), "); forcing ",
            "BPPARAM = BiocParallel::SerialParam() to avoid multiple ",
            "processes contending for one GPU device.", call. = FALSE)
    BPPARAM <- BiocParallel::SerialParam()
  }
  if (is.null(block.size) && gpu_active) {
    block.size <- .inferenceBlockSize(length(genes), n, ncol(fit@W), backend,
                                      gpu.mem.budget)
  }
  if (is.null(block.size)) block.size <- 500L

  # the design is shared by every block, so it crosses to the device once
  W_dev <- if (gpu_active) {
    SpaNorm::toGPUMatrix(fit@W, backend = backend)
  } else {
    fit@W
  }

  idx <- split(genes, ceiling(seq_along(genes) / block.size))
  parts <- BiocParallel::bplapply(idx, function(g) {
    gmb <- fit@gmean[g]
    ab <- fit@alpha[g, , drop = FALSE]
    psib <- fit@psi[g]

    if (gpu_active) {
      Yb <- SpaNorm::toGPUMatrix(as.matrix(Y[g, , drop = FALSE]),
                                 backend = backend)
      mub <- SpaNorm::calculateMu(gmb, ab, W_dev, backend = backend)
      # NB variance mu + psi mu^2 == mu (1 + psi mu); the vec-mat helper does
      # the per-gene broadcast that plain recycling cannot on a tensor
      v <- mub * (1 + SpaNorm::mult_vec_mat_gpu(psib, mub, backend = backend))
      r <- (Yb - mub) / torch::torch_sqrt(v)
      rmean <- as.numeric(SpaNorm::toRMatrix(
        SpaNorm::rowSums_gpu(r))) / n
      rc <- SpaNorm::add_vec_mat_gpu(-rmean, r, backend = backend)
      ss <- as.numeric(SpaNorm::toRMatrix(SpaNorm::rowSums_gpu(rc * rc)))
      rsd <- sqrt(ss / (n - 1))
      ok <- is.finite(rsd) & rsd > 0
      # zeroing the scale drops a degenerate gene's contribution entirely,
      # which is the same thing as excluding its row but needs no masking
      inv <- ifelse(ok, 1 / rsd, 0)
      z <- SpaNorm::mult_vec_mat_gpu(inv, rc, backend = backend)
      sblk <- as.numeric(SpaNorm::toRMatrix(torch::torch_sum(z, dim = 1L)))
    } else {
      Yb <- as.matrix(Y[g, , drop = FALSE])
      mub <- as.matrix(SpaNorm::toRMatrix(
        SpaNorm::calculateMu(gmb, ab, W_dev)))
      r <- (Yb - mub) / sqrt(mub + psib * mub^2)
      rmean <- rowMeans(r)
      rc <- r - rmean
      rsd <- sqrt(rowSums(rc^2) / (n - 1))
      # A gene with no residual variation has an undefined correlation with
      # anything; dropping it is the only defensible choice, since keeping it
      # propagates NaN through the whole estimate.
      ok <- is.finite(rsd) & rsd > 0
      z <- rc * ifelse(ok, 1 / rsd, 0)
      z[!is.finite(z)] <- 0
      sblk <- colSums(z)
    }
    list(s = sblk, kept = sum(ok))
  }, BPPARAM = BPPARAM)

  s <- Reduce(`+`, lapply(parts, `[[`, "s"))
  kept <- sum(vapply(parts, function(x) as.integer(x$kept), integer(1)))
  .meanCorFromZ(s, n, kept)
}

#' Per-set sums of a statistic matrix
#'
#' \code{rowsum()} does the grouped accumulation at C level in one pass, which
#' matters because the naive version is a loop over tens of thousands of sets.
#'
#' @param mat a genes x k matrix.
#' @param sets a list of integer row-index vectors.
#' @return a length(sets) x k matrix of column means per set.
#' @noRd
.setColMeans <- function(mat, sets) {
  gi <- unlist(sets, use.names = FALSE)
  grp <- rep(seq_along(sets), lengths(sets))
  sums <- rowsum(mat[gi, , drop = FALSE], group = grp, reorder = FALSE)
  sums / lengths(sets)
}

#' Relative log-likelihood weights across bandwidths
#'
#' Shared with \code{.geneWeights()}; a gene set's weight is built from the
#' summed log-likelihood of its member genes, so the same transform has to be
#' applied to a matrix that is not per-gene.
#'
#' @param ll a rows x bandwidth matrix of log-likelihoods.
#' @param thresh weights below this are set to 0.
#' @importFrom matrixStats rowMaxs
#' @noRd
.relWeights <- function(ll, thresh = 0.1) {
  w <- exp(ll - matrixStats::rowMaxs(ll))
  w[w < thresh] <- 0
  w
}

#' Nested BH cascade over gene sets
#'
#' Mirrors \code{.hierarchicalFDR()}: each level is tested only within the
#' survivors of the level above, so the reported q-values are conditional on
#' that gating rather than being three independent marginal corrections.
#'
#' @param tab a data.frame with geneset, collection, ct_index, ct_niche, p.
#' @param fdr the FDR threshold.
#' @param nested TRUE for the 3-level (niche) cascade, FALSE for 2-level.
#' @importFrom stats p.adjust
#' @noRd
.gseaCascade <- function(tab, fdr, nested = TRUE) {
  cauchy_by <- function(p, key) {
    sp <- split(p, key)
    cb <- vapply(sp, function(x) .cauchyCombine(matrix(x, nrow = 1L)), numeric(1))
    unname(cb[key])
  }
  bh_by <- function(p, key) {
    unsplit(lapply(split(p, key), stats::p.adjust, method = "BH"), key)
  }

  # level 1: the set as a whole, corrected within its collection
  # (each level re-derives its keys from the surviving rows, because the
  # previous level's filtering invalidates any index computed before it)
  k1 <- paste(tab$collection, tab$geneset, sep = "\r")
  tab$p.geneset <- cauchy_by(tab$p, k1)
  lv1 <- !duplicated(k1)
  q1 <- stats::p.adjust(tab$p.geneset[lv1], "BH")
  tab$fdr.geneset <- q1[match(k1, k1[lv1])]
  tab <- tab[tab$fdr.geneset < fdr, , drop = FALSE]
  if (!nrow(tab)) return(tab)

  # level 2: index cell type, within surviving sets
  k2 <- paste(tab$collection, tab$geneset, tab$ct_index, sep = "\r")
  tab$p.index <- cauchy_by(tab$p, k2)
  lv2 <- !duplicated(k2)
  sub <- tab[lv2, , drop = FALSE]
  q2 <- unsplit(lapply(split(sub$p.index, sub$geneset), stats::p.adjust,
                       method = "BH"), sub$geneset)
  tab$fdr.index <- q2[match(k2, k2[lv2])]
  tab <- tab[tab$fdr.index < fdr, , drop = FALSE]
  if (!nrow(tab) || !nested) return(tab)

  # level 3: niche cell type, within surviving (set, index) pairs
  k3 <- paste(tab$collection, tab$geneset, tab$ct_index, sep = "\r")
  tab$fdr.niche <- bh_by(tab$p, k3)
  tab[tab$fdr.niche < fdr, , drop = FALSE]
}

#' Gene-set enrichment on a fitted spiDE model
#'
#' Tests whether the genes of a set shift coherently in their response-niche
#' (or cell-type response) statistics. Gene-set statistics average over tens to
#' hundreds of genes, which recovers the power the per-gene niche tests lose on
#' sparse spatial data, with the inter-gene correlation carried explicitly so
#' that co-regulation is not mistaken for evidence.
#'
#' Evidence is combined across bandwidths with the same log-likelihood-weighted
#' Cauchy combination used for the gene-level results, then gated by a nested
#' Benjamini-Hochberg cascade: set level, then index cell type within surviving
#' sets, then (for \code{type = "niche"}) niche cell type within those.
#'
#' @section Which test, and what it licenses you to say:
#' \code{test = "self-contained"} (the default) asks whether the set's mean
#' statistic differs from **zero** -- null: no gene in the set responds. It is
#' the more powerful of the two and reproduces the established pipeline, but a
#' significant result does **not** mean the set is special: if a global effect
#' shifts most genes, most sets become significant, correctly but uselessly.
#'
#' \code{test = "competitive"} asks whether the set's mean differs from the
#' mean of the genes **outside** it, as \code{limma::camera} does -- null: the
#' set responds no more than the assayed background. Use it when the claim is
#' "this pathway specifically", not merely "this pathway responds".
#'
#' This matters concretely here: spiDE's per-gene statistics track expression,
#' so a set of abundant genes can carry a real non-zero mean without carrying
#' niche-specific signal. The competitive test absorbs much of that, because
#' its background is subject to the same bias. Where the distinction is load
#' bearing, also compare each set against expression-matched random sets of
#' the same size rather than reading the q-value alone.
#'
#' @param object a [SpiDEResults] from [testSpiDE()] or [spiDE()].
#' @param spe the [SpatialExperiment::SpatialExperiment] the model was fitted
#'   on; needed for the counts the inter-gene correlation is estimated from.
#' @param genesets a named list of gene identifiers (character, matched against
#'   \code{rownames}) or integer row indices.
#' @param type "niche" (default) tests the three-way celltype:condition:niche
#'   statistics; "celltype" tests the CellType:condition statistics, which are
#'   empty unless the design carries that block.
#' @param test "self-contained" (default) compares the set's mean statistic
#'   with zero; "competitive" compares it with the genes outside the set. See
#'   the section below -- they answer different questions and the default is
#'   the more permissive one.
#' @param fdr the FDR threshold applied at every level of the cascade.
#' @param min.size,max.size sets outside this size range (after intersecting
#'   with the fitted genes) are dropped. The lower bound keeps the mean from
#'   being dominated by one gene; the upper bound drops sets so broad that the
#'   competitive null they are tested against is largely themselves.
#' @param rho the average inter-gene correlation: \code{NULL} to estimate it
#'   from Pearson residuals separately for each bandwidth (each is a different
#'   design and leaves different residuals), a scalar to use one value
#'   throughout, or one value per bandwidth. Pass 0 only if you are certain the
#'   genes are independent -- they are not, and the test is severely
#'   anti-conservative without this term.
#' @param rho.genes number of genes to subsample when estimating \code{rho}
#'   (\code{NULL} uses all). A random subset of genes gives a random subset of
#'   gene PAIRS, so the estimate is unbiased; the default trades a little
#'   precision on a nuisance scalar for a large saving on the counts pass.
#' @param block.size genes per block when streaming the counts. \code{NULL}
#'   auto-sizes against the device memory budget on the GPU backend.
#' @param backend the compute backend for the inter-gene correlation pass
#'   ("auto", "cpu" or "gpu"), which is the only step that touches the counts
#'   and so the only one worth accelerating -- the gene-set arithmetic
#'   downstream is a sets x contrasts matrix and negligible by comparison.
#'   Ignored when \code{rho} is supplied, since no pass is then needed.
#' @param gpu.mem.budget device memory budget in bytes for the GPU block sizer.
#' @param BPPARAM a BiocParallelParam. Gene blocks contribute additive partial
#'   sums, so they parallelise exactly; a multi-worker param is downgraded to
#'   serial on the GPU backend to stop workers contending for one device.
#' @param verbose report progress.
#' @return a data.frame of significant sets, keyed by (geneset, ct_index) and
#'   additionally ct_niche when \code{type = "niche"}, with the set size, the
#'   z-statistic at its most informative bandwidth, the direction, and the
#'   q-value from each level of the cascade. Empty if nothing survives.
#' @references Wu D, Smyth GK (2012). "Camera: a competitive gene set test
#'   accounting for inter-gene correlation." \emph{Nucleic Acids Research}
#'   40(17):e133. (Source of the variance inflation factor, and of the
#'   two-sample form used by \code{test = "competitive"}.)
#' @references Wu D, Lim E, Vaillant F, Asselin-Labat ML, Visvader JE, Smyth GK
#'   (2010). "ROAST: rotation gene set tests for complex microarray
#'   experiments." \emph{Bioinformatics} 26(17):2176-2182. (Self-contained
#'   gene-set null.)
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = 20)
#' res <- spiDE(spe, condition = "condition", sigma = 20, verbose = FALSE)
#' gs <- list(setA = rownames(spe)[1:5], setB = rownames(spe)[6:10])
#' spiGSEA(res, spe, gs, min.size = 3, fdr = 1)
#' @importFrom stats pnorm p.adjust
#' @importFrom matrixStats colVars
#' @rdname spiGSEA
#' @export
setMethod(
  "spiGSEA", "SpiDEResults",
  function(object, spe, genesets, type = c("niche", "celltype"),
           test = c("self-contained", "competitive"),
           fdr = 0.05, min.size = 5L, max.size = 500L, rho = NULL,
           rho.genes = 2000L, block.size = NULL,
           backend = c("auto", "cpu", "gpu"), gpu.mem.budget = NULL,
           BPPARAM = BiocParallel::SerialParam(), verbose = TRUE) {
    type <- match.arg(type)
    test <- match.arg(test)
    backend <- match.arg(backend)
    fl <- fits(object)
    if (!length(fl)) stop("`object` holds no fits")
    f1 <- fl[[1]]
    if (is.null(f1@t_stat) || !length(f1@t_stat)) {
      stop("no inference on `object`; run testSpiDE() first")
    }
    gnames <- rownames(f1@t_stat)
    Y <- SummarizedExperiment::assay(spe, "counts")
    if (!identical(nrow(Y), f1@ngenes) && !all(gnames %in% rownames(Y))) {
      stop("`spe` does not match the fitted genes")
    }

    # --- map sets onto fitted rows, then size-filter ------------------------
    sets <- lapply(genesets, function(g) {
      if (is.character(g)) which(gnames %in% g) else as.integer(g)
    })
    sizes <- lengths(sets)
    keep <- sizes >= min.size & sizes <= max.size
    if (!any(keep)) {
      stop("no gene set has between ", min.size, " and ", max.size,
           " of the fitted genes")
    }
    sets <- sets[keep]
    if (verbose) {
      message(sprintf("spiGSEA: %d of %d sets within [%d, %d] genes",
                      length(sets), length(genesets), min.size, max.size))
    }

    # --- inter-gene correlation, PER BANDWIDTH -----------------------------
    # Each bandwidth is a different design, so it leaves different residuals
    # and a different average inter-gene correlation. Reusing one bandwidth's
    # rho for all of them mis-scales every other bandwidth's set statistics,
    # in a direction that depends on how the niche columns happen to soak up
    # correlated structure -- so it is estimated per fit unless supplied.
    gsub_ <- seq_along(gnames)
    if (!is.null(rho.genes) && rho.genes < length(gsub_)) {
      gsub_ <- sort(sample(gsub_, rho.genes))
    }
    rhov <- if (is.null(rho)) NULL else rep_len(rho, length(fl))
    if (verbose && is.null(rhov)) {
      message(sprintf("spiGSEA: estimating inter-gene correlation on %d genes x %d bandwidths",
                      length(gsub_), length(fl)))
    }

    # --- per-bandwidth set statistics --------------------------------------
    want <- if (type == "niche") "ResponseNiche" else "ResponseCellType"
    per <- lapply(seq_along(fl), function(i) {
      f <- fl[[i]]
      ct <- as.character(f@covtype)
      respcols <- grepl("Response", ct)
      sel <- ct[respcols] == want
      if (!any(sel)) {
        return(NULL)
      }
      tm <- f@t_stat[, sel, drop = FALSE]
      dfn <- if (is.null(f@df) || length(f@df) == 1L) f@df else f@df[sel]
      zm <- .tToZ(tm, dfn)
      zbar <- .setColMeans(zm, sets)
      m <- lengths(sets)
      rho_i <- if (is.null(rhov)) {
        .interGeneCor(Y, f, gsub_, block.size, backend = backend,
                      gpu.mem.budget = gpu.mem.budget, BPPARAM = BPPARAM)
      } else {
        rhov[i]
      }
      if (verbose) {
        message(sprintf("spiGSEA: sigma %g  rho = %.4f", f@sigma, rho_i))
      }
      # equicorrelation variance inflation: the mean of m unit-variance
      # statistics with average correlation rho
      vif <- 1 + rho_i * (m - 1)
      stat <- if (test == "self-contained") {
        zbar / sqrt(vif / m)
      } else {
        # camera's two-sample form: the set's mean against the mean of
        # everything outside it, on the observed spread of the statistics
        # rather than an assumed unit variance. Standardising by what the data
        # actually show is what makes this robust to a global shift affecting
        # every gene -- which a self-contained test would report as every set
        # being significant.
        G <- nrow(zm)
        mOut <- G - m
        # mean of the genes NOT in the set, from the column total minus the
        # set's own sum -- avoids materialising a complement per set
        meanOut <- sweep(-(zbar * m), 2L, colSums(zm), "+") / mOut
        vpool <- matrixStats::colVars(zm)
        sweep(zbar - meanOut, 2L, sqrt(vpool), "/") /
          sqrt(vif / m + 1 / mOut)
      }
      list(z = stat, rho = rho_i,
           ll = vapply(sets, function(g) sum(f@loglik[g]), numeric(1)),
           cm = f@coefmap[respcols, , drop = FALSE][sel, , drop = FALSE],
           sigma = f@sigma)
    })
    per <- Filter(Negate(is.null), per)
    if (!length(per)) {
      stop("the fit carries no ", want, " columns; ",
           if (type == "celltype")
             "the design has no CellType:condition block" else
             "check the design")
    }

    zl <- lapply(per, `[[`, "z")
    sig <- vapply(per, `[[`, numeric(1), "sigma")
    gene.w <- .relWeights(do.call(cbind, lapply(per, `[[`, "ll")))

    # --- combine bandwidths on two-sided p ---------------------------------
    ns <- nrow(zl[[1]]); nc <- ncol(zl[[1]])
    pcomb <- matrix(NA_real_, ns, nc)
    for (j in seq_len(nc)) {
      pj <- vapply(zl, function(z) 2 * stats::pnorm(-abs(z[, j])), numeric(ns))
      # With a SINGLE gene set vapply returns a length-k vector rather than a
      # 1 x k matrix, and as.matrix() inside .cauchyCombine would then stand it
      # up as k x 1 -- the transpose of what it expects, against a 1 x k weight
      # matrix, which errors outright. Testing one set is an ordinary thing to
      # do, so pin the shape rather than rely on vapply's simplification.
      if (!is.matrix(pj)) pj <- matrix(pj, nrow = ns)
      pcomb[, j] <- .cauchyCombine(pj, gene.w)
    }
    # Report the statistic at the bandwidth that saw it most strongly, and take
    # the direction from there: a set can genuinely shift in opposite directions
    # at different scales, and averaging the signed z would hide that.
    az <- abs(zl[[1]]); zbest <- zl[[1]]
    bw <- matrix(sig[1], ns, nc)
    for (i in seq_along(zl)[-1]) {
      hit <- abs(zl[[i]]) > az
      az[hit] <- abs(zl[[i]])[hit]
      zbest[hit] <- zl[[i]][hit]
      bw[hit] <- sig[i]
    }

    cm <- per[[1]]$cm
    coll <- attr(genesets, "collection")
    collv <- if (!is.null(coll)) coll[names(sets)] else rep("all", length(sets))
    collv[is.na(collv)] <- "all"

    tab <- data.frame(
      geneset = rep(names(sets), times = nc),
      collection = rep(unname(collv), times = nc),
      size = rep(lengths(sets), times = nc),
      ct_index = rep(as.character(cm$index), each = ns),
      ct_niche = rep(as.character(cm$niche), each = ns),
      bandwidth.max = as.numeric(bw),
      z = as.numeric(zbest),
      p = as.numeric(pcomb),
      stringsAsFactors = FALSE
    )
    tab$Direction <- ifelse(tab$z > 0, "Up", "Down")

    out <- .gseaCascade(tab, fdr, nested = (type == "niche"))
    if (!nrow(out)) {
      if (verbose) message("spiGSEA: no set survived the cascade")
      return(out)
    }
    cols <- c("geneset", "collection", "size", "ct_index",
              if (type == "niche") "ct_niche", "bandwidth.max", "z",
              "Direction", "fdr.geneset", "fdr.index",
              if (type == "niche") "fdr.niche")
    out <- out[, intersect(cols, names(out)), drop = FALSE]
    key <- if (type == "niche") out$fdr.niche else out$fdr.index
    out <- out[order(key, -abs(out$z)), , drop = FALSE]
    rownames(out) <- NULL
    # the correlation actually used scales every statistic here, so it travels
    # with the result rather than only appearing in a log
    attr(out, "rho") <- stats::setNames(
      vapply(per, `[[`, numeric(1), "rho"), paste0("sigma", sig))
    attr(out, "test") <- test
    if (verbose) {
      message(sprintf("spiGSEA: %d significant rows over %d sets",
                      nrow(out), length(unique(out$geneset))))
    }
    out
  }
)

#' Average inter-gene correlation of a fitted model's residuals
#'
#' The mean pairwise correlation of Pearson residuals, per bandwidth. This is
#' the \code{rho} that [spiGSEA()] uses to inflate the variance of a gene-set
#' mean, and it is worth having in its own right as a model diagnostic: it says
#' how much of the residual variation is shared across genes, and therefore how
#' far a set of \code{m} genes falls short of carrying \code{m} genes' worth of
#' independent evidence.
#'
#' Compute it once and pass it to [spiGSEA()] as \code{rho} when testing
#' several gene-set collections against the same fit. It depends only on the
#' fit, not on the sets, so re-estimating it per collection repeats identical
#' work. It is deliberately **not** cached on the fitted object: R's copy
#' semantics mean a cache would have to be written to a copy and reassigned by
#' the caller anyway, and a silently stale correlation would mis-scale every
#' gene-set p-value derived from it.
#'
#' Estimated without ever forming the gene-by-gene correlation matrix, so it
#' costs one streaming pass over the counts (about 5 seconds per bandwidth on
#' 13,000 genes x 77,000 cells with 16 workers) and memory linear in cells.
#'
#' @param object a [SpiDEResults].
#' @param spe the [SpatialExperiment::SpatialExperiment] the model was fitted
#'   on.
#' @param rho.genes number of genes to subsample (\code{NULL} uses all). A
#'   random subset of genes gives a random subset of gene PAIRS, so the
#'   estimate is unbiased; the default trades a little precision for a large
#'   saving on the counts pass.
#' @param block.size genes per block; \code{NULL} auto-sizes on the GPU
#'   backend.
#' @param backend "auto", "cpu" or "gpu".
#' @param gpu.mem.budget device memory budget in bytes for the GPU block sizer.
#' @param BPPARAM a BiocParallelParam. Blocks contribute additive partial sums,
#'   so they parallelise exactly.
#' @param verbose report progress.
#' @return a named numeric, one value per bandwidth.
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = 20)
#' res <- spiDE(spe, condition = "condition", sigma = 20, verbose = FALSE)
#' interGeneCor(res, spe)
#' @rdname interGeneCor
#' @export
setMethod(
  "interGeneCor", "SpiDEResults",
  function(object, spe, rho.genes = 2000L, block.size = NULL,
           backend = c("auto", "cpu", "gpu"), gpu.mem.budget = NULL,
           BPPARAM = BiocParallel::SerialParam(), verbose = FALSE) {
    backend <- match.arg(backend)
    fl <- fits(object)
    if (!length(fl)) stop("`object` holds no fits")
    Y <- SummarizedExperiment::assay(spe, "counts")
    g <- seq_len(nrow(Y))
    if (!is.null(rho.genes) && rho.genes < length(g)) {
      g <- sort(sample(g, rho.genes))
    }
    out <- vapply(fl, function(f) {
      r <- .interGeneCor(Y, f, g, block.size, backend = backend,
                         gpu.mem.budget = gpu.mem.budget, BPPARAM = BPPARAM)
      if (verbose) message(sprintf("sigma %g: rho = %.4f", f@sigma, r))
      r
    }, numeric(1))
    stats::setNames(out, paste0("sigma", vapply(fl, function(f) f@sigma,
                                                numeric(1))))
  }
)
