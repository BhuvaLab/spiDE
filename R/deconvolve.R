# Undo transcript spillover by inverting the mixing operator.
#
# Segmentation bleed and ambient RNA move transcripts from a cell into its
# spatial neighbours. Model it as the linear operator every decontamination
# method assumes: a cell keeps a fraction (1 - kappa) of its own transcripts and
# donates kappa, spread over the cells within `radius`,
#
#     Y_obs = X A'    with    A[i, j] = (1 - kappa) delta_ij + kappa w_ij
#
# where w_ij is the row-normalised neighbour indicator (w_ii = 0, rows summing
# to 1; isolated cells keep everything). A is sparse, banded by the radius, and
# strictly diagonally dominant for kappa < 1/2 -- hence invertible, and solvable
# by a sparse LU per sample rather than a dense inverse.
#
# Recovering X is then one solve:
#
#     X = Y_obs (A')^{-1}   <=>   A X' = Y_obs'
#
# WHY THE EXACT SOLVE AND NOT A FIRST-ORDER CORRECTION. The obvious cheap
# version, X ~ (Y - kappa W Y)/(1 - kappa), is the first Neumann term and is
# accurate only to O(kappa^2). On simulated data at kappa = 0.3 that residual
# matters a great deal: the first-order correction took a contaminated
# niche-only analysis from 68 false calls to 24, while the exact solve took the
# same data from 70 to 1 -- the clean-data count -- and recovered all three
# planted effects that the uncontaminated analysis had missed. Do not replace
# this with the cheap form.
#
# WHAT COMES BACK IS NOT COUNT DATA, STRICTLY. (A')^{-1} is a sharpening filter
# with negative off-diagonal entries, so the solve returns reals and some are
# negative. They are rounded and floored at zero. Measured on the simulation
# (kappa = 0.3, 60 genes, 900 cells): 8.0% of entries went negative, the most
# negative was -1.04, and flooring destroyed 0.46% of the total count mass. On a
# sparser fixture 35% of entries floored -- but 97.7% of those had an OBSERVED
# count of zero, where the inverse merely overshoots a true zero and clipping is
# correct. Hence two QC numbers: `floored` tracks sparsity, `floored_mass`
# tracks damage, and only the latter is worth warning on.
# Reassuringly the mean-variance relationship is RESTORED rather than damaged --
# neighbour mixing is a smoother, so contamination suppressed the dispersion
# (var/mean 4.64 clean -> 3.70 contaminated) and inverting it brings it back
# (-> 4.54). The floored fraction is returned as a QC value precisely because
# this argument stops holding as kappa grows and the operator becomes less well
# conditioned.

#' Build the row-normalised neighbour operator for one sample
#'
#' Sparse by construction: only pairs within \code{radius} are stored. Cells
#' with no neighbour inside the radius get an all-zero row, so \code{A} keeps a
#' 1 on their diagonal and they pass through untouched.
#'
#' @param xy a cells x 2 matrix of coordinates.
#' @param radius the leak radius, in the units of \code{xy}.
#' @return a sparse \code{dgCMatrix}, rows summing to 1 or 0.
#' @importFrom Matrix sparseMatrix rowSums Diagonal
#' @noRd
.neighbourOperator <- function(xy, radius) {
  n <- nrow(xy)
  # spatstat's closepairs is the sparse neighbour query the package already
  # depends on; it returns only pairs within the radius, never an n x n matrix.
  pp <- spatstat.geom::ppp(xy[, 1], xy[, 2],
                           window = spatstat.geom::owin(range(xy[, 1]),
                                                        range(xy[, 2])),
                           check = FALSE)
  cp <- spatstat.geom::closepairs(pp, radius, what = "indices")
  if (!length(cp$i)) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                dims = c(n, n)))
  }
  # closepairs returns each pair in both directions, so tabulating the first
  # index gives the degree directly -- no need to materialise W to count.
  deg <- tabulate(cp$i, nbins = n)
  Matrix::sparseMatrix(i = cp$i, j = cp$j, x = 1 / deg[cp$i],
                       dims = c(n, n))
}

#' Remove transcript spillover by inverting the mixing operator
#'
#' Reverses segmentation bleed and ambient contamination, which move transcripts
#' from a cell into its spatial neighbours. This matters most for a
#' condition-free ([fitSpiDE()] with `condition = NULL`) analysis: contamination
#' is delivered in proportion to the local density of each neighbouring cell
#' type, which is exactly the shape of the `CellType:niche` covariate being
#' tested, so it is confounded with the biology at the point of estimation. A
#' conditioned analysis is largely protected — spillover carries no condition
#' label — so this is optional there.
#'
#' @section Model:
#' A cell keeps \eqn{1 - \kappa} of its transcripts and donates \eqn{\kappa},
#' spread evenly over the cells within `radius`:
#' \deqn{Y_{obs} = X A', \quad A_{ij} = (1-\kappa)\delta_{ij} + \kappa w_{ij}}
#' with \eqn{w} the row-normalised neighbour indicator. Recovering \eqn{X} is
#' one sparse solve per sample. \eqn{A} is strictly diagonally dominant for
#' \eqn{\kappa < 1/2}, so the solve is stable in that range.
#'
#' @section Choosing kappa:
#' There is no estimator in this function — `kappa` is yours to supply. The
#' error is strongly asymmetric, so **err high**: on simulated data with a true
#' \eqn{\kappa} of 0.30, correcting at 0.40 left 3 false calls (and recovered
#' every planted effect) while correcting at 0.20 left 47. Overshooting costs
#' little; undershooting leaves most of the artefact. [spilloverScore()]
#' evaluated over a grid of `kappa` is a usable objective — it should be
#' minimised near the truth.
#'
#' @param spe a SpatialExperiment with spatial coordinates and a counts assay.
#' @param kappa the leaked fraction, in \[0, 0.5). Zero returns `spe` unchanged.
#' @param radius the leak radius in the units of `spatialCoords(spe)`, normally
#'   about one cell diameter.
#' @param assay a character, the assay to correct.
#' @param name a character, the assay to write. Deliberately not `assay`: the
#'   operation is lossy and the original should stay available for comparison.
#' @param sample_id a character, the colData column identifying samples. Cells
#'   are only ever mixed within a sample.
#' @param det.drop.max the per-gene detection loss above which a gene is
#'   flagged as bimodalised by the filter (default 0.02). The cohort's median
#'   loss is 0.000 and its 99th percentile 0.006, so this sits well clear of
#'   normal behaviour while catching the technical outliers.
#' @param BPPARAM a BiocParallelParam; samples are independent.
#' @param verbose a logical, report the per-sample floored fraction.
#'
#' @return `spe` with the deconvolved counts in `assay(spe, name)` and a QC
#'   summary in `metadata(spe)$spiDE_deconvolution`: the `kappa` and `radius`
#'   used, plus two flooring diagnostics.
#'
#'   `floored` is the fraction of *entries* that came back negative and were
#'   truncated at zero. Read it as a measure of sparsity, not of damage: on a
#'   toy fixture 35% of entries floored, but **97.7% of those had an observed
#'   count of zero** — the inverse overshoots slightly negative where a gene was
#'   simply not detected, and clipping to zero is the right answer there. Expect
#'   this number to be large on sparse panels.
#'
#'   `floored_mass` is the fraction of total count *mass* destroyed by that
#'   truncation, and is the number that actually bounds the damage (0.46% on the
#'   simulation study, 6.9% on a low-count toy). Above ~10% the operator is
#'   outside the range where this approximation was validated and a warning is
#'   raised — treat the result with suspicion rather than confidence.
#'
#'   `det_drop` is the per-gene loss of detection (fraction of cells with a
#'   non-zero count) caused by the correction, and `amplified` names the genes
#'   exceeding `det.drop.max`. The
#'   operator is the inverse of a smoother, hence a sharpener, so it amplifies
#'   heavy tails while flooring the middle — which shows up precisely as lost
#'   detection. Real panels carry technical outliers (one gene in the YTMA
#'   cohort has a mean of 8 and a maximum of 7711) that this turns into bimodal
#'   distributions no negative-binomial fit can handle; left in, they dominate
#'   every downstream test. **Filter them before fitting.** Simulated NB counts
#'   have no such tail, which is why this failure mode appears only on real
#'   data.
#'
#' @examples
#' data(toySpiDE)
#' spe <- deconvolveSpillover(toySpiDE, kappa = 0.2, radius = 20)
#' SummarizedExperiment::assayNames(spe)
#' S4Vectors::metadata(spe)$spiDE_deconvolution$floored
#'
#' @rdname deconvolveSpillover
#' @importFrom Matrix solve t
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom S4Vectors metadata metadata<-
#' @export
setMethod(
  "deconvolveSpillover", "ANY",
  function(spe, kappa, radius, assay = "counts", name = "counts_deconv",
           sample_id = "sample_id", det.drop.max = 0.02,
           BPPARAM = BiocParallel::SerialParam(), verbose = TRUE, ...) {
    checkSPE(spe, assay = assay, sample_id = sample_id)
    if (!is.numeric(kappa) || length(kappa) != 1L || is.na(kappa) ||
        kappa < 0 || kappa >= 0.5) {
      stop("'kappa' should be a single number in [0, 0.5); the mixing operator ",
           "stops being diagonally dominant at 0.5 and the solve is no longer ",
           "stable.", call. = FALSE)
    }
    if (!is.numeric(radius) || length(radius) != 1L || is.na(radius) ||
        radius <= 0) {
      stop("'radius' should be a single positive number", call. = FALSE)
    }

    Y <- SummarizedExperiment::assay(spe, assay)
    if (kappa == 0) {
      SummarizedExperiment::assay(spe, name) <- Y
      S4Vectors::metadata(spe)[["spiDE_deconvolution"]] <-
        list(kappa = 0, radius = radius, floored = 0, floored_mass = 0)
      return(spe)
    }

    xy <- SpatialExperiment::spatialCoords(spe)
    smp <- as.character(SummarizedExperiment::colData(spe)[[sample_id]])
    groups <- split(seq_along(smp), smp)

    res <- BiocParallel::bplapply(names(groups), function(s) {
      j <- groups[[s]]
      W <- .neighbourOperator(xy[j, , drop = FALSE], radius)
      deg <- Matrix::rowSums(W)
      # isolated cells keep everything: their row of W is empty, so the
      # diagonal must be 1 rather than 1 - kappa
      dg <- ifelse(deg > 0, 1 - kappa, 1)
      A <- Matrix::Diagonal(x = dg) + kappa * W
      # A X' = Y'  ->  one sparse LU, applied to every gene at once
      Xt <- Matrix::solve(A, Matrix::t(as.matrix(Y[, j, drop = FALSE])))
      X <- t(as.matrix(Xt))
      Xr <- pmax(round(X), 0)
      # Per-gene DETECTION, before and after. The operator is a sharpener (the
      # inverse of a smoother): on a heavy-tailed gene it pushes the tail
      # further out while flooring the middle to zero, which is exactly a drop
      # in the fraction of cells with a non-zero count. That turned out to be a
      # far better discriminator than amplification of the maximum -- see the
      # note above the function.
      Yj <- as.matrix(Y[, j, drop = FALSE])
      list(j = j, X = Xr, nneg = sum(X < 0), n = length(X),
           negmass = -sum(X[X < 0]), tot = sum(Yj),
           det_r = rowSums(Yj > 0), det_z = rowSums(Xr > 0), ncell = length(j))
    }, BPPARAM = BPPARAM)

    # Assemble sparse. Densifying the whole assay is not an option at realistic
    # size: 13k genes x 77k cells is ~8 GB as doubles, while the deconvolved
    # counts are as sparse as the originals.
    ord <- unlist(groups, use.names = FALSE)
    out <- do.call(cbind, lapply(res, function(r) {
      methods::as(Matrix::Matrix(r$X, sparse = TRUE), "dgCMatrix")
    }))
    out <- out[, order(ord), drop = FALSE]
    dimnames(out) <- dimnames(Y)
    gt <- function(f) sum(vapply(res, `[[`, numeric(1), f))
    floored <- gt("nneg") / gt("n")
    floored_mass <- gt("negmass") / gt("tot")

    # Genes the sharpening filter has bimodalised. Real panels carry technical
    # outliers (one gene in the YTMA cohort has a mean of 8 and a maximum of
    # 7711); the inverse pushes their tail further out while flooring their
    # middle, and the result is a distribution no negative-binomial fit can
    # handle, which then dominates every downstream test. Simulated NB counts
    # have no such tail, which is why this appears only on real data.
    #
    # DETECTION DROP is the diagnostic, chosen by measurement. On the cohort it
    # ranked the known-bad genes 1st, 2nd and 7th of 10,422 against a median of
    # 0.000 and a 99th percentile of 0.006. Inflation of the per-gene MAXIMUM
    # was tried first and is useless here: deconvolution raises every gene's
    # maximum by roughly the same 1/(1-kappa), so the worst offender sat well
    # inside the range of well-behaved genes.
    nc <- sum(vapply(res, `[[`, numeric(1), "ncell"))
    det_r <- Reduce(`+`, lapply(res, `[[`, "det_r")) / nc
    det_z <- Reduce(`+`, lapply(res, `[[`, "det_z")) / nc
    det_drop <- det_r - det_z
    names(det_drop) <- rownames(Y)
    amplified <- det_drop > det.drop.max

    if (verbose) {
      message(sprintf(
        "deconvolveSpillover: kappa = %.3f, radius = %g; %.1f%% of entries floored (%.2f%% of count mass)",
        kappa, radius, 100 * floored, 100 * floored_mass))
      if (any(amplified)) {
        worst <- names(sort(det_drop[amplified], decreasing = TRUE))
        warning(sum(amplified), " gene(s) lost more than ",
                round(100 * det.drop.max, 1), "% of their detection to the ",
                "correction (worst: ", paste(utils::head(worst, 5), collapse = ", "),
                "). These are heavy-tailed, usually technical, genes that the ",
                "inverse filter bimodalises; drop them before fitting or the ",
                "corrected data will be worse than the uncorrected.",
                call. = FALSE)
      }
      if (floored_mass > 0.10) {
        warning("flooring destroyed more than 10% of the count mass; this is ",
                "outside the range where the approximation was validated -- ",
                "check 'kappa' and 'radius'", call. = FALSE)
      }
    }

    SummarizedExperiment::assay(spe, name) <- out
    S4Vectors::metadata(spe)[["spiDE_deconvolution"]] <-
      list(kappa = kappa, radius = radius, floored = floored,
           floored_mass = floored_mass, det_drop = det_drop,
           amplified = names(det_drop)[amplified])
    spe
  }
)
