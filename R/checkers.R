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

# Patient-level checks for the mixed-effects (random-effects) fit. The condition
# must be a patient-level variable (constant within each sample), otherwise the
# per-sample random intercept is mis-specified, and there must be enough samples
# for the between-sample variance components to be identifiable.
checkSample <- function(spe, condition = NULL, sample_id = "sample_id",
                        covariates = character()) {
  cd <- SummarizedExperiment::colData(spe)
  if (!sample_id %in% colnames(cd)) {
    stop(sprintf("sample id column '%s' not found in colData(spe)", sample_id))
  }
  smp <- as.character(cd[[sample_id]])
  # A condition-free (niche-only) design has no condition to be patient-level,
  # so this check applies only when one was supplied.
  if (!is.null(condition)) {
    cond <- as.character(cd[[condition]])
    n_lvl <- tapply(cond, smp, function(x) length(unique(x[!is.na(x)])))
    if (any(n_lvl > 1)) {
      stop(sprintf(
        "condition '%s' varies within sample(s): %s. The random-effects fit needs a patient-level condition (constant within '%s').",
        condition, paste(names(n_lvl)[n_lvl > 1], collapse = ", "), sample_id
      ))
    }
  }
  # sample-constant covariates are confounded with the per-sample random
  # intercept (which already adjusts for all between-sample nuisance variation)
  const <- covariates[vapply(covariates, function(cv) {
    all(tapply(as.character(cd[[cv]]), smp,
               function(x) length(unique(x[!is.na(x)]))) <= 1)
  }, logical(1))]
  if (length(const) > 0) {
    stop(sprintf(
      "covariate(s) constant within sample: %s. With random='intercept'/'slope' the per-sample random intercept already absorbs all between-sample effects, so drop these sample-level covariates.",
      paste(const, collapse = ", ")
    ))
  }
  n_samples <- length(unique(smp))
  if (n_samples < 3) {
    warning(sprintf(
      "only %d sample(s); random-effect variance components may be unreliable",
      n_samples
    ))
  }
  invisible(TRUE)
}

# fdr must be a single value in (0, 1]; 1 is allowed as the "show everything"
# threshold (see testSpiDE()'s documentation).
checkFdr <- function(fdr) {
  if (!is.numeric(fdr) || length(fdr) != 1 || is.na(fdr) || fdr <= 0 || fdr > 1) {
    stop("'fdr' should be a single numeric value in (0, 1]")
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
