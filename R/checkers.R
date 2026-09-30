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
checkCounts <- function(Y, integer.only = FALSE) {
  mn <- suppressWarnings(min(Y, na.rm = TRUE))
  if (!is.finite(mn) || mn < 0) {
    stop("counts should be non-negative")
  }
  # The negative binomial likelihood is defined on counts. dnbinom() returns
  # -Inf for a non-integer value, so the per-gene Newton fit would reject every
  # step and -- while maximising a constant -Inf -- return the dispersion
  # optimiser's upper bound for EVERY gene, with its own diagnostics reporting
  # success (measured on spiDE 0.99.17's polish: psi 999.96 across the board).
  # Refuse it here, before the fit, rather than after an hour of work.
  if (integer.only) {
    ss <- as.numeric(Y[seq_len(min(nrow(Y), 20L)), , drop = FALSE])
    ss <- ss[is.finite(ss)]
    if (length(ss) && max(abs(ss - round(ss))) > 1e-8) {
      stop("fitSpiDE() requires integer counts, and this assay is not ",
           "integer-valued: the negative binomial likelihood is undefined off ",
           "the integers, so every gene's dispersion would collapse ",
           "to its upper bound.\n  Supply raw counts (assay = \"counts\").",
           call. = FALSE)
    }
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
checkCovariates <- function(spe, covariates, finite.only = FALSE) {
  if (length(covariates) == 0) {
    return(invisible(TRUE))
  }
  cd <- SummarizedExperiment::colData(spe)
  missing <- setdiff(covariates, colnames(cd))
  if (length(missing) > 0) {
    stop(sprintf("covariate(s) not found in colData(spe): %s", paste(missing, collapse = ", ")))
  }
  # model.matrix() drops rows with a missing value, so the design comes back
  # shorter than the patient block built from the full-length sample
  # labels and the two fail to cbind with "number of rows of matrices must
  # match" -- an error that says nothing about which covariate is at fault. A
  # non-finite value is the usual cause and is easy to produce by accident:
  # log() of a zero-valued QC column gives -Inf, and centring that gives NaN.
  #
  # OPT-IN, because it is only the GLM design that breaks this way.
  # a patient-level estimator would deliberately DROP patients with a missing patient-level
  # covariate and reports the dropout, which is a documented behaviour its
  # tests assert; rejecting there would remove a feature.
  if (!finite.only) return(invisible(TRUE))
  cd <- SummarizedExperiment::colData(spe)
  bad <- covariates[vapply(covariates, function(cv) {
    x <- cd[[cv]]
    is.numeric(x) && !all(is.finite(x))
  }, logical(1))]
  if (length(bad) > 0) {
    stop(sprintf(
      "covariate(s) with missing or non-finite values: %s. Every cell needs a finite value (log() of a zero-valued column is the usual cause -- use log1p, or drop the affected cells).",
      paste(bad, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

# Patient-level checks. The condition must be a patient-level variable
# (constant within each sample): spiDE tests it as a between-patient contrast of
# within-patient slopes. Every index type's fit has an intercept per patient,
# which absorbs any covariate constant within a sample, so those are rejected.
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
        "condition '%s' varies within sample(s): %s. spiDE needs a patient-level condition (constant within '%s').",
        condition, paste(names(n_lvl)[n_lvl > 1], collapse = ", "), sample_id
      ))
    }
  }
  # sample-constant covariates are confounded with the per-patient intercepts
  # (which already adjust for all between-patient nuisance variation)
  const <- covariates[vapply(covariates, function(cv) {
    all(tapply(as.character(cd[[cv]]), smp,
               function(x) length(unique(x[!is.na(x)]))) <= 1)
  }, logical(1))]
  if (length(const) > 0) {
    stop(sprintf(
      "covariate(s) constant within sample: %s. The per-patient intercepts already absorb all between-patient effects, so drop these patient-level covariates.",
      paste(const, collapse = ", ")
    ))
  }
  n_samples <- length(unique(smp))
  if (n_samples < 3) {
    warning(sprintf(
      "only %d sample(s); the between-patient tests need several patients per condition",
      n_samples
    ))
  }
  invisible(TRUE)
}

# fdr must be a single value in (0, 1]; 1 is allowed as the "show everything"
# threshold (see testSpiDE()'s documentation).
# A colData column named by an optional argument (NULL = not used).
checkColumn <- function(spe, column, what) {
  if (is.null(column)) return(invisible(TRUE))
  if (!is.character(column) || length(column) != 1L ||
      !column %in% colnames(SummarizedExperiment::colData(spe))) {
    stop(sprintf("%s column '%s' not found in colData(spe)", what, paste(column, collapse = ", ")),
         call. = FALSE)
  }
  invisible(TRUE)
}

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

# The patient intercepts a plot or accessor reads; fits from spiDE <= 0.99.33
# did not keep them.
checkIntercepts <- function(xi, what) {
  if (is.null(xi$intercept)) {
    stop(sprintf("%s needs the patient intercepts that fitSpiDE() keeps from spiDE 0.99.34 on: ", what),
         "refit this object with the installed spiDE", call. = FALSE)
  }
  invisible(TRUE)
}

# A plotting function's input: a current SpiDEResults (or, where the plot needs
# only the fit, a SpiDEFit).
checkResults <- function(x, what, fit.ok = FALSE) {
  if (!(is(x, "SpiDEResults") || (fit.ok && is(x, "SpiDEFit")))) {
    stop(sprintf("%s needs a SpiDEResults from testSpiDE() or spiDE()%s", what,
                 if (fit.ok) ", or a SpiDEFit from fitSpiDE()" else ""), call. = FALSE)
  }
  .assertCurrent(x)
  invisible(TRUE)
}

# Views of each patient's own slopes exist only for the slopes engine.
checkSlopesEngine <- function(fit, what) {
  if (fit@engine != "slopes") {
    stop(sprintf("%s needs the slopes engine's per-patient slopes: fit with engine = \"slopes\"", what),
         call. = FALSE)
  }
  invisible(TRUE)
}

# An index type the fit has, one of its tested niches, genes it tested.
checkTriplet <- function(fit, index, niche = NULL, gene = NULL) {
  if (!is.character(index) || length(index) != 1L || !index %in% names(fit@index)) {
    stop(sprintf("'index' must be one of the fitted index types: %s", paste(names(fit@index), collapse = ", ")),
         call. = FALSE)
  }
  xi <- fit@index[[index]]
  if (!is.null(niche) && (!is.character(niche) || length(niche) != 1L || !niche %in% xi$niches)) {
    stop(sprintf("'niche' must be one of the niches tested in %s: %s", index, paste(xi$niches, collapse = ", ")),
         call. = FALSE)
  }
  if (!is.null(gene)) {
    if (!is.character(gene) || !length(gene)) stop("'gene' must be a character vector", call. = FALSE)
    miss <- setdiff(gene, xi$genes)
    if (length(miss)) {
      stop(sprintf("gene(s) not tested in %s: %s", index, paste(miss, collapse = ", ")), call. = FALSE)
    }
  }
  invisible(TRUE)
}
