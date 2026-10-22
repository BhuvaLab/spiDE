# Niche covariate construction. For each biological sample and each bandwidth,
# a Gaussian kernel density estimate of every cell type is evaluated at each
# cell's location, giving a cells x cellTypes matrix of "effective niche"
# densities. Reimplements effective_niche()/fill_missing_dims() from
# the original spiDE analysis scripts, standardised on BiocParallel.

#' Compute the effective niche of a single sample
#'
#' Gaussian kernel density of each cell type evaluated at every cell, scaled by
#' \code{pi * sigma^2} so values are on the scale of expected neighbour counts.
#'
#' @param coords a cells x 2 coordinate matrix of one sample.
#' @param cts a character vector, the cell type of each cell of the sample.
#' @param sigma a numeric, the kernel bandwidth (standard deviation).
#' @param edge,diggle logicals passed to \code{spatstat.explore::densityfun}.
#' @return a cells x cellTypes matrix of niche densities (columns in order of
#'   first appearance in \code{cts}).
#' @noRd
.effectiveNiche <- function(coords, cts, sigma, edge = TRUE, diggle = TRUE) {
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

    cts <- as.character(SummarizedExperiment::colData(spe)[[cell_type]])
    all_cell_types <- sort(unique(cts))
    coords <- SpatialExperiment::spatialCoords(spe)
    sids <- as.character(SummarizedExperiment::colData(spe)[[sample_id]])
    # cell indices per sample; only coordinates and cell types go to the workers
    idx <- split(seq_len(ncol(spe)), factor(sids, levels = unique(sids)))

    per_sample <- BiocParallel::bplapply(
      idx,
      function(i, coords, cts, sigma, edge, diggle, all_cell_types) {
        lapply(sigma, function(sg) {
          m <- .effectiveNiche(coords[i, , drop = FALSE], cts[i], sigma = sg, edge = edge, diggle = diggle)
          m <- .fillMissingDims(m, all_cell_types, "column")
          m[is.na(m)] <- 0
          m
        })
      },
      coords = coords, cts = cts, sigma = sigma, edge = edge, diggle = diggle,
      all_cell_types = all_cell_types, BPPARAM = BPPARAM
    )

    for (j in seq_along(sigma)) {
      mat <- matrix(0, ncol(spe), length(all_cell_types),
                    dimnames = list(colnames(spe), all_cell_types))
      for (s in names(idx)) mat[idx[[s]], ] <- per_sample[[s]][[j]]
      SingleCellExperiment::reducedDim(spe, paste0(name, sigma[j])) <- mat
    }

    return(spe)
  }
)
