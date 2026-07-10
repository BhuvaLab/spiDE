# Per-cell-type sample size factors. For each (sample, cell type) a log size
# factor is computed from total counts / total area, centred per cell type, then
# the median across cell types is taken as the per-sample library-size offset
# (LS). Reimplements the size-factor block in batch_nichede_v9.R.

#' Compute the per-(sample, cell-type) log size factor matrix
#'
#' @param spe a SpatialExperiment.
#' @param cell_type a character, the colData column of cell type labels.
#' @param sample_id a character, the colData column of sample identifiers.
#' @param count a character, the colData column of per-cell total counts.
#' @param area a character, the colData column of per-cell area.
#' @return a samples x cellTypes matrix of centred log size factors.
#' @noRd
.computeLS <- function(spe, cell_type, sample_id, count, area) {
  cd <- as.data.frame(SummarizedExperiment::colData(spe))
  cts <- sort(unique(as.character(cd[[cell_type]])))
  sids <- unique(as.character(cd[[sample_id]]))

  sf <- vapply(cts, function(ct) {
    sub <- cd[cd[[cell_type]] == ct, , drop = FALSE]
    fs <- tapply(sub[[count]], factor(sub[[sample_id]], levels = sids), sum)
    ar <- tapply(sub[[area]], factor(sub[[sample_id]], levels = sids), sum)
    ls <- log2(fs / ar)
    ls - mean(ls, na.rm = TRUE) # centre per cell type
  }, numeric(length(sids)))
  rownames(sf) <- sids
  colnames(sf) <- cts
  return(sf)
}

#' Compute per-sample library-size offsets
#'
#' Derives a per-sample library-size offset (LS) by computing, for each
#' (sample, cell type), a log size factor from total counts over total area,
#' centring these per cell type, and taking the median across cell types for
#' each sample. The offset is written to a per-cell \code{colData} column so it
#' can enter the model design as a nuisance covariate.
#'
#' @param spe a SpatialExperiment.
#' @param cell_type a character, the colData column of cell type labels.
#' @param sample_id a character, the colData column of sample identifiers.
#' @param count a character, the colData column of per-cell total counts
#'   (e.g. total UMI or negative-probe counts).
#' @param area a character, the colData column of per-cell area.
#' @param name a character, the colData column to write the LS offset to
#'   (default "LS").
#' @param ... ignored.
#'
#' @return the input \code{spe} with a per-cell \code{name} column in colData.
#'
#' @examples
#' data(toySpiDE)
#' spe <- toySpiDE
#' spe <- computeSizeFactors(spe, count = "nCount", area = "Area")
#' head(spe$LS)
#'
#' @rdname computeSizeFactors
#' @export
setMethod(
  "computeSizeFactors",
  signature = "ANY",
  definition = function(spe, cell_type = "cell_type", sample_id = "sample_id",
                        count = "nCount", area = "Area", name = "LS", ...) {
    checkSPE(spe, cell_type = cell_type, sample_id = sample_id)
    cd <- SummarizedExperiment::colData(spe)
    for (col in c(count, area)) {
      if (!col %in% colnames(cd)) {
        stop(sprintf("column '%s' not found in colData(spe)", col))
      }
    }

    sf <- .computeLS(spe, cell_type, sample_id, count, area)
    ls_samples <- apply(sf, 1, stats::median, na.rm = TRUE)

    sids <- as.character(cd[[sample_id]])
    SummarizedExperiment::colData(spe)[[name]] <- unname(ls_samples[sids])
    return(spe)
  }
)
