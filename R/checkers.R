# Input validation helpers. These fail early with informative messages so the
# user-facing functions can assume well-formed inputs.

#' @importFrom methods is
#' @importFrom SummarizedExperiment assayNames assay colData
#' @importFrom SingleCellExperiment reducedDimNames
checkSPE <- function(spe, assay = "counts", cell_type = "cell_type", sample_id = "sample_id") {
  if (!is(spe, "SpatialExperiment")) {
    stop("'spe' should be a SpatialExperiment object")
  }
  if (!assay %in% assayNames(spe)) {
    stop(sprintf("assay '%s' not found in 'spe'", assay))
  }
  cd <- colData(spe)
  if (!cell_type %in% colnames(cd)) {
    stop(sprintf("cell type column '%s' not found in colData(spe)", cell_type))
  }
  if (!sample_id %in% colnames(cd)) {
    stop(sprintf("sample id column '%s' not found in colData(spe)", sample_id))
  }
  invisible(TRUE)
}

# Counts must be non-negative (NB GLM assumption). Only a min() reduction is
# forced, which is cheap and DelayedArray-friendly (no full realisation).
checkCounts <- function(Y) {
  mn <- suppressWarnings(min(Y, na.rm = TRUE))
  if (!is.finite(mn) || mn < 0) {
    stop("counts should be non-negative")
  }
  invisible(TRUE)
}

# The condition column must exist and have exactly two levels.
checkCondition <- function(spe, condition) {
  cd <- SummarizedExperiment::colData(spe)
  if (!condition %in% colnames(cd)) {
    stop(sprintf("condition column '%s' not found in colData(spe)", condition))
  }
  vals <- cd[[condition]]
  lvls <- unique(vals[!is.na(vals)])
  if (length(lvls) != 2) {
    stop(sprintf("condition '%s' should have exactly two levels, found %d", condition, length(lvls)))
  }
  invisible(TRUE)
}

# Nuisance covariates must be present in colData.
checkCovariates <- function(spe, covariates) {
  if (length(covariates) == 0) {
    return(invisible(TRUE))
  }
  cd <- SummarizedExperiment::colData(spe)
  missing <- setdiff(covariates, colnames(cd))
  if (length(missing) > 0) {
    stop(sprintf("covariate(s) not found in colData(spe): %s", paste(missing, collapse = ", ")))
  }
  invisible(TRUE)
}

# The requested niche reducedDim must have been built.
checkNiche <- function(spe, sigma, name = "Niche") {
  nms <- SingleCellExperiment::reducedDimNames(spe)
  need <- paste0(name, sigma)
  missing <- setdiff(need, nms)
  if (length(missing) > 0) {
    stop(sprintf(
      "niche reducedDim(s) not found: %s. Run buildNiches() first.",
      paste(missing, collapse = ", ")
    ))
  }
  invisible(TRUE)
}
