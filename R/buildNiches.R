# Niche covariate construction. For each biological sample and each bandwidth,
# a Gaussian kernel density estimate of every cell type is evaluated at each
# cell's location, giving a cells x cellTypes matrix of "effective niche"
# densities. Reimplements effective_niche()/speapply()/fill_missing_dims() from
# the original spiDE analysis scripts, standardised on BiocParallel.

#' Compute the effective niche of a single sample
#'
#' Gaussian kernel density of each cell type evaluated at every cell, scaled by
#' \code{pi * sigma^2} so values are on the scale of expected neighbour counts.
#'
#' @param x a SpatialExperiment for a single sample.
#' @param sigma a numeric, the kernel bandwidth (standard deviation).
#' @param cell_type a character, the colData column holding cell type labels.
#' @param edge,diggle logicals passed to \code{spatstat.explore::densityfun}.
#' @return a cells x cellTypes matrix of niche densities.
#' @importFrom SpatialExperiment spatialCoords
#' @noRd
.effectiveNiche <- function(x, sigma, cell_type = "cell_type", edge = TRUE, diggle = TRUE) {
  coords <- SpatialExperiment::spatialCoords(x)
  cts <- as.character(SummarizedExperiment::colData(x)[[cell_type]])

  xrange <- range(coords[, 1])
  yrange <- range(coords[, 2])
  window <- spatstat.geom::owin(xrange = xrange, yrange = yrange)
  ppp_obj <- spatstat.geom::ppp(
    coords[, 1], coords[, 2],
    window = window, marks = factor(cts)
  )

  uct <- unique(cts)
  res <- lapply(uct, function(ct) {
    dfun <- spatstat.explore::densityfun(
      ppp_obj[ppp_obj$marks == ct],
      sigma = sigma, edge = edge, diggle = diggle
    )
    dfun(coords[, 1], coords[, 2]) * pi * sigma^2
  })
  res <- do.call(cbind, res)
  colnames(res) <- uct
  rownames(res) <- colnames(x)
  return(res)
}

#' Pad a matrix so its columns (or rows) match a reference set
#'
#' Absent columns/rows are added (filled with \code{NA}) and the result is
#' reordered to the reference. Used to reconcile per-sample niche matrices that
#' may be missing cell types absent from that sample.
#'
#' @param mat a matrix.
#' @param ref_names a character vector of reference names.
#' @param dim one of "column" or "row".
#' @return a matrix with dimnames matching \code{ref_names} along \code{dim}.
#' @noRd
.fillMissingDims <- function(mat, ref_names, dim = c("column", "row")) {
  dim <- match.arg(dim)
  if (dim == "column") {
    missing_names <- setdiff(ref_names, colnames(mat))
    if (length(missing_names) > 0) {
      pad <- matrix(
        NA_real_, nrow = nrow(mat), ncol = length(missing_names),
        dimnames = list(rownames(mat), missing_names)
      )
      mat <- cbind(mat, pad)
    }
    mat <- mat[, ref_names, drop = FALSE]
  } else {
    missing_names <- setdiff(ref_names, rownames(mat))
    if (length(missing_names) > 0) {
      pad <- matrix(
        NA_real_, nrow = length(missing_names), ncol = ncol(mat),
        dimnames = list(missing_names, colnames(mat))
      )
      mat <- rbind(mat, pad)
    }
    mat <- mat[ref_names, , drop = FALSE]
  }
  return(mat)
}

#' Split a SpatialExperiment by sample and apply a function per sample
#'
#' @param spe a SpatialExperiment.
#' @param FUN a function applied to each single-sample subset.
#' @param ... further arguments to FUN.
#' @param sample_id a character, the colData column identifying samples.
#' @param BPPARAM a BiocParallelParam.
#' @return a named list of per-sample results.
#' @importFrom BiocParallel bplapply SerialParam
#' @noRd
.speApply <- function(spe, FUN, ..., sample_id = "sample_id", BPPARAM = BiocParallel::SerialParam()) {
  sids <- unique(as.character(SummarizedExperiment::colData(spe)[[sample_id]]))
  spl <- lapply(sids, function(s) spe[, as.character(SummarizedExperiment::colData(spe)[[sample_id]]) == s])
  names(spl) <- sids
  BiocParallel::bplapply(spl, FUN, ..., BPPARAM = BPPARAM)
}

#' Build niche covariate matrices
#'
#' Computes, for each biological sample and each bandwidth \code{sigma}, a
#' Gaussian kernel density estimate of every cell type evaluated at each cell's
#' location. The resulting cells x cellTypes matrices are stored as
#' \code{reducedDim(spe, "Niche<sigma>")}, one per bandwidth, with a consistent
#' set of cell-type columns (cell types absent from a sample are set to 0).
#'
#' @param spe a SpatialExperiment with spatial coordinates, cell type labels,
#'   and sample identifiers.
#' @param sigma a numeric vector of kernel bandwidths (in the units of
#'   \code{spatialCoords}), default \code{c(10, 30, 50, 70)}.
#' @param cell_type a character, the colData column holding cell type labels.
#' @param sample_id a character, the colData column identifying samples.
#' @param edge,diggle logicals passed to \code{spatstat.explore::densityfun}
#'   for edge correction (default TRUE).
#' @param name a character, the prefix for the stored reducedDims (default
#'   "Niche").
#' @param BPPARAM a BiocParallelParam for parallelising over samples.
#' @param ... ignored.
#'
#' @return the input \code{spe} with one \code{reducedDim} per bandwidth added.
#'
#' @examples
#' data(toySpiDE)
#' spe <- toySpiDE
#' spe <- buildNiches(spe, sigma = c(10, 20))
#' SingleCellExperiment::reducedDimNames(spe)
#'
#' @rdname buildNiches
#' @importFrom SingleCellExperiment reducedDim reducedDim<- reducedDimNames
#' @importFrom BiocParallel SerialParam
#' @export
setMethod(
  "buildNiches",
  signature = "ANY",
  definition = function(spe, sigma = c(10, 30, 50, 70), cell_type = "cell_type",
                        sample_id = "sample_id", edge = TRUE, diggle = TRUE,
                        name = "Niche", BPPARAM = BiocParallel::SerialParam(), ...) {
    checkSPE(spe, cell_type = cell_type, sample_id = sample_id)

    all_cell_types <- sort(unique(as.character(SummarizedExperiment::colData(spe)[[cell_type]])))

    # per-sample list, each element a list of per-sigma matrices
    per_sample <- .speApply(
      spe,
      function(x) {
        lapply(sigma, function(sg) {
          m <- .effectiveNiche(x, sigma = sg, cell_type = cell_type, edge = edge, diggle = diggle)
          m <- .fillMissingDims(m, all_cell_types, "column")
          m[is.na(m)] <- 0
          m
        })
      },
      sample_id = sample_id, BPPARAM = BPPARAM
    )

    # .speApply splits by unique(sample_id) in first-appearance order and each
    # subset preserves cell order, so rbind lands cells in `perm` order. Invert
    # that permutation to restore the original cell order without relying on
    # colnames(spe) (SPEs frequently have none).
    sids <- as.character(SummarizedExperiment::colData(spe)[[sample_id]])
    perm <- order(match(sids, unique(sids)))
    inv <- integer(length(perm))
    inv[perm] <- seq_along(perm)

    for (i in seq_along(sigma)) {
      mat <- do.call(rbind, lapply(per_sample, function(s) s[[i]]))
      mat <- mat[inv, , drop = FALSE]
      rownames(mat) <- colnames(spe)
      SingleCellExperiment::reducedDim(spe, paste0(name, sigma[i])) <- mat
    }

    return(spe)
  }
)
