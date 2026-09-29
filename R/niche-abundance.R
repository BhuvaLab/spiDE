# The between-patient niche-abundance association, tested at the patient level.
#
# spiDE's niche slopes are within-patient slopes: fitSpiDE() gives every patient
# its own intercept in each index cell type, so the tests compare cells WITHIN a
# patient. What those intercepts absorb is a real association in the data:
# patients whose type-k cells sit, on average, in a denser type-n neighbourhood
# may also have a different mean expression in type k. It is not
# neighbourhood-dependent differential expression (it was found in 0.99.16 when
# it leaked into the mixed model's niche slopes and inflated the triplet FDR on
# the real cohort; research/fdr-ordering/REPORT.md, section 5d). But "do
# patients whose tumour compartment is fibroblast-rich express this gene
# differently IN tumour cells?" is a legitimate question with S experimental
# units, and this is where it is asked -- at the patient level, on pseudobulk
# means, with limma, so the standard error is the between-patient one that the
# question needs. (Up to spiDE 0.99.22 this function was compositionTest().)

#' Per-(sample, index type) pseudobulk profiles and mean niche densities
#'
#' @return a list with, per index type k: \code{Y} (genes x samples of
#'   log2-CPM), \code{niche} (samples x niche types of mean log1p density),
#'   \code{ncells} (samples), all restricted to samples with at least
#'   \code{min.cells} cells of type k.
#' @noRd
.pseudobulkByIndex <- function(Y, NM, ct, smp, index, min.cells, prior.count) {
  cts <- if (is.null(index)) sort(unique(ct)) else intersect(.sanitise(index), sort(unique(ct)))
  out <- list()
  for (k in cts) {
    ik <- which(ct == k)
    tab <- table(smp[ik])
    keep <- names(tab)[tab >= min.cells]
    if (length(keep) < 3L) next
    grp <- factor(smp[ik], levels = keep)
    ok <- !is.na(grp)
    ik <- ik[ok]; grp <- grp[ok]
    # pseudobulk: sum counts over the type-k cells of each sample, then
    # log2-CPM against that sample's own type-k library size
    # a product with a sparse sample indicator, so the counts are never
    # densified whole (only genes x samples is)
    ind <- Matrix::sparseMatrix(i = seq_along(grp), j = as.integer(grp), x = 1,
                                dims = c(length(grp), nlevels(grp)))
    sums <- as.matrix(Y[, ik, drop = FALSE] %*% ind)                 # genes x samples
    colnames(sums) <- levels(grp)
    lib <- colSums(sums)
    lcpm <- log2(t((t(sums) + prior.count) / (lib + 2 * prior.count)) * 1e6)
    nm <- rowsum(log1p(NM[ik, , drop = FALSE]), group = grp, reorder = TRUE) /
      as.numeric(table(grp))
    out[[k]] <- list(Y = lcpm, niche = nm, ncells = as.numeric(table(grp)),
                     samples = keep)
  }
  out
}

#' Test whether patients with more of a niche cell type express genes differently
#'
#' A between-patient question that the within-patient engines cannot see:
#' does a gene's expression in index cell type \eqn{k} differ between patients
#' whose type-\eqn{k} cells sit, on average, among more type-\eqn{n} cells?
#' ([fitSpiDE()] and [testSpiDE()] instead ask whether expression changes with
#' niche density WITHIN each patient.)
#'
#' For each index cell type \eqn{k} and niche cell type \eqn{n}, regresses the
#' per-sample pseudobulk log2-CPM of every gene in the type-\eqn{k} cells on
#' the per-sample mean log1p niche density of type \eqn{n} around those cells,
#' across samples, with a \code{limma} moderated \eqn{t}-test. With a
#' \code{condition} the design is \code{~ niche * condition + covariates} and
#' two terms are reported: \code{"niche"}, the association pooled across
#' conditions, and \code{"condition:niche"}, its difference between conditions
#' -- the patient-level counterpart of the condition-specific niche test that
#' [testSpiDE()] runs within patients. Without a condition only
#' \code{"niche"} is reported.
#'
#' This is the association that [fitSpiDE()]'s per-patient intercepts
#' deliberately absorb. It is a between-patient effect with \eqn{S} experimental units; it is not
#' neighbourhood-dependent differential expression, and the two must not be
#' conflated (see the model vignette).
#' It was \code{compositionTest()} up to spiDE 0.99.22.
#' Samples contributing fewer than
#' \code{min.cells} cells of the index type are dropped for that index type,
#' and an index type with fewer than three remaining samples is skipped.
#'
#' @param spe a SpatialExperiment with niche reducedDims (see [buildNiches()]).
#' @param condition a character, the colData column of the condition (constant
#'   within sample), or \code{NULL}.
#' @param sigma a numeric, the bandwidth (one value).
#' @param index,niche character vectors restricting the index / niche cell
#'   types (NULL = all). An index type is never tested against its own niche.
#' @param covariates a character vector of sample-level colData columns to
#'   adjust for (constant within sample; the first non-missing value per
#'   sample is used).
#' @param assay a character, the counts assay.
#' @param cell_type,sample_id the colData columns of cell type and sample.
#' @param name the niche reducedDim prefix.
#' @param min.cells minimum cells of the index type a sample must contribute.
#' @param prior.count the pseudocount in the log2-CPM.
#' @param verbose report progress.
#' @param ... further arguments passed to the method.
#' @return a data.frame with one row per (gene, index type, niche type, term):
#'   \code{gene}, \code{index}, \code{niche}, \code{term}, \code{estimate}
#'   (log2-CPM per unit log1p density), \code{t}, \code{p},
#'   \code{n_patients}, \code{q} (BH within each (index, niche, term) over
#'   genes) and \code{q.global} (BH over every row).
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = 30)
#' na <- testNicheAbundance(spe, condition = "condition", sigma = 30)
#' head(na[order(na$p), ])
#' @importFrom limma lmFit eBayes
#' @importFrom stats model.matrix p.adjust
#' @rdname testNicheAbundance
#' @export
setMethod(
  "testNicheAbundance",
  signature = "ANY",
  definition = function(spe, condition = NULL, sigma, index = NULL, niche = NULL,
                        covariates = character(), assay = "counts",
                        cell_type = "cell_type", sample_id = "sample_id",
                        name = "Niche", min.cells = 10L, prior.count = 1,
                        verbose = TRUE) {
    checkSPE(spe, assay = assay, cell_type = cell_type, sample_id = sample_id)
    if (!is.null(condition)) checkCondition(spe, condition)
    checkCovariates(spe, covariates)
    if (length(sigma) != 1L) stop("'sigma' must be a single bandwidth")
    checkNiche(spe, sigma, name = name)

    cd <- SummarizedExperiment::colData(spe)
    ct <- .sanitise(as.character(cd[[cell_type]]))
    smp <- as.character(cd[[sample_id]])
    Y <- SummarizedExperiment::assay(spe, assay)
    NM <- as.matrix(SingleCellExperiment::reducedDim(spe, paste0(name, sigma)))
    colnames(NM) <- .sanitise(colnames(NM))
    niches <- if (is.null(niche)) colnames(NM) else intersect(.sanitise(niche), colnames(NM))
    if (!length(niches)) stop("no requested niche cell types found in the niche reducedDim")
    NM <- NM[, niches, drop = FALSE]

    pb <- .pseudobulkByIndex(Y, NM, ct, smp, index, min.cells, prior.count)
    if (!length(pb)) {
      stop("no index cell type has at least three samples with >= min.cells cells")
    }
    # sample-level condition and covariates: first non-missing value per sample
    pats <- unique(smp)
    cond_s <- if (is.null(condition)) NULL else
      factor(.patientValue(as.character(cd[[condition]]), smp, pats))
    cov_s <- lapply(covariates, function(cv) .patientValue(cd[[cv]], smp, pats))
    names(cov_s) <- covariates

    rows <- list()
    for (k in names(pb)) {
      b <- pb[[k]]
      for (n in setdiff(niches, k)) {
        if (verbose) message(sprintf("testNicheAbundance: %s x %s (%d samples)", k, n, length(b$samples)))
        df <- data.frame(niche = as.numeric(b$niche[b$samples, n]))
        if (!is.null(cond_s)) df$condition <- droplevels(cond_s[b$samples])
        for (cv in covariates) df[[cv]] <- cov_s[[cv]][b$samples]
        f <- if (is.null(cond_s)) "~ niche" else "~ niche * condition"
        if (length(covariates)) f <- paste(f, "+", paste(covariates, collapse = " + "))
        X <- tryCatch(stats::model.matrix(stats::as.formula(f), df), error = function(e) NULL)
        # both conditions must be present with at least two samples each for
        # the interaction to be estimable; otherwise report the pooled term only
        if (is.null(X) || nrow(X) < ncol(X) + 2L || nrow(X) < nrow(df)) {
          if (!is.null(cond_s) && nrow(df) >= 4L) {
            f <- "~ niche"; if (length(covariates)) f <- paste(f, "+", paste(covariates, collapse = " + "))
            X <- stats::model.matrix(stats::as.formula(f), df)
            if (nrow(X) < ncol(X) + 2L) next
          } else next
        }
        fit <- limma::lmFit(b$Y, design = X)
        fit <- limma::eBayes(fit, robust = nrow(X) >= 6L)
        terms <- intersect(c("niche", grep("^niche:condition", colnames(X), value = TRUE)), colnames(X))
        for (tm in terms) {
          rows[[length(rows) + 1L]] <- data.frame(
            gene = rownames(b$Y), index = k, niche = n,
            term = if (tm == "niche") "niche" else "condition:niche",
            estimate = fit$coefficients[, tm], t = fit$t[, tm], p = fit$p.value[, tm],
            n_patients = nrow(X), row.names = NULL, stringsAsFactors = FALSE)
        }
      }
    }
    if (!length(rows)) stop("no (index, niche) pair had enough samples to fit")
    out <- do.call(rbind, rows)
    out <- out[is.finite(out$p), , drop = FALSE]
    key <- paste(out$index, out$niche, out$term)
    out$q <- stats::ave(out$p, key, FUN = function(p) stats::p.adjust(p, "BH"))
    out$q.global <- stats::p.adjust(out$p, "BH")
    rownames(out) <- NULL
    out
  }
)

# ---- patient-level helpers (shared with the archived two-stage estimator) -----

#' First non-missing value of a cell-level vector, per patient
#'
#' Indexes the original vector rather than going through tapply(), which
#' unlists a factor into bare integer level codes -- a factor patient
#' covariate would then enter the stage-2 design as a continuous trend in
#' arbitrary codes, silently. Taking the first NON-missing cell also keeps a
#' patient whose first cell happens to be NA while its value is known
#' elsewhere.
#' @noRd
.patientValue <- function(x, pat, pats) {
  idx <- vapply(pats, function(p) {
    i <- which(pat == p & !is.na(x))
    if (length(i)) i[1L] else which(pat == p)[1L]
  }, integer(1))
  stats::setNames(x[idx], pats)
}
