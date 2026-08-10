# Two-stage estimation: per-patient niche slopes, then a patient-level contrast.
#
# WHY THIS IS NOT A `random` MODE
# -------------------------------
# `random = "none"/"intercept"/"slope"` all fit the SAME design matrix and differ
# only in which columns are ridge-penalised. This is a different ESTIMATOR: the
# niche slope is estimated within each patient and the condition contrast is then
# taken over patients. None of `W`, `alpha`, `psi`, `penalty`, `tau2` or the
# Satterthwaite df machinery applies, so slotting it under `random` would hand
# every downstream consumer an object whose fields do not mean what they claim.
#
# WHAT IT FIXES
# -------------
# `condition` is assigned per PATIENT, so patients are the experimental units.
# The one-stage model treats cells as replicates of a patient-level contrast; on
# a patient-label permutation of the YTMA cohort (55 patients, every triplet null
# by construction) that returns ~576 false calls in every replicate, realized FDR
# 1.00 against a nominal 0.05. Two-stage on the same permuted data gives raw
# type-I 0.036 and ZERO false calls, in 1.5 minutes against 128.
#
# WHAT IT DOES NOT FIX
# --------------------
# Multiplicity. The estimator carries real signal (true triplets enriched 5x at
# alpha 0.05 and 65x at alpha 0.001 on a spiked plasmode) but a full-panel
# hypothesis space of ~1.8M triplets buries it: BH needs the best p below
# 0.05/m, and the hierarchical cascade does not help either, because ACAT over a
# gene's ~137 mostly-null triplets dilutes a single true effect to no better than
# Bonferroni. Restrict `index`, `niche` and the gene set to a pre-specified
# hypothesis; ~4,000 tests is the order at which a p ~ 1e-5 effect survives.

#' Per-patient niche slopes (stage 1)
#'
#' For each patient and index cell type, regresses each gene's log-CPM on each
#' niche density across that patient's cells of that type. Vectorised as one
#' matrix product per (patient, index): all genes and niches at once.
#'
#' @param E genes x cells log-CPM matrix.
#' @param nm cells x niches density matrix.
#' @param ct,pat per-cell index cell type and patient labels.
#' @param idx_types index cell types to fit.
#' @param min.cells minimum cells for a patient to contribute a slope.
#' @return a named list of (genes x niches x patients) arrays, one per index.
#' @noRd
.patientSlopes <- function(E, nm, ct, pat, idx_types, min.cells) {
  pats <- sort(unique(pat))
  niches <- colnames(nm)
  out <- setNames(vector("list", length(idx_types)), idx_types)
  for (ix in idx_types) {
    A <- array(NA_real_, c(nrow(E), length(niches), length(pats)),
               dimnames = list(rownames(E), niches, pats))
    for (pp in seq_along(pats)) {
      k <- which(ct == ix & pat == pats[pp])
      if (length(k) < min.cells) next
      Xc <- scale(nm[k, , drop = FALSE], center = TRUE, scale = FALSE)
      ss <- colSums(Xc^2)
      ok <- ss > 1e-8
      if (!any(ok)) next
      Ec <- E[, k, drop = FALSE]
      Ec <- Ec - rowMeans(Ec)
      A[, ok, pp] <- (Ec %*% Xc[, ok, drop = FALSE]) %*%
        diag(1 / ss[ok], sum(ok))
    }
    out[[ix]] <- A
  }
  out
}

#' Patient-level contrast of the per-patient slopes (stage 2)
#'
#' With no covariates this is a Welch two-sample t-test on the patient slopes --
#' valid without assuming equal variances, and the right test when the two
#' groups may differ in how precisely their slopes are estimated.
#'
#' With `patient.covariates` it becomes a linear model on the patient slopes.
#' That is the capability the one-stage path CANNOT offer: `checkSample()`
#' rejects sample-constant covariates under `random != "none"`, because they are
#' collinear with the per-sample random intercept. Here the unit of analysis IS
#' the patient, so Sex, Stage, Histology and treatment arm are ordinary
#' covariates -- which matters on cohorts where composition is confounded with
#' the condition.
#' @noRd
.patientContrast <- function(sl, grp, cov_df = NULL) {
  ok <- !is.na(sl)
  if (sum(ok) < 4L) return(c(NA_real_, NA_real_, NA_real_))
  y <- sl[ok]; g <- grp[ok]
  if (length(unique(g)) < 2L) return(c(NA_real_, NA_real_, NA_real_))
  if (is.null(cov_df)) {
    a <- y[g == levels(g)[1]]; b <- y[g != levels(g)[1]]
    if (length(a) < 2L || length(b) < 2L) return(c(NA_real_, NA_real_, NA_real_))
    va <- stats::var(a); vb <- stats::var(b)
    se <- sqrt(va / length(a) + vb / length(b))
    if (!is.finite(se) || se <= 0) return(c(NA_real_, NA_real_, NA_real_))
    est <- mean(b) - mean(a)
    df <- (va / length(a) + vb / length(b))^2 /
      ((va / length(a))^2 / (length(a) - 1) + (vb / length(b))^2 / (length(b) - 1))
    tt <- est / se
    return(c(est, tt, 2 * stats::pt(-abs(tt), df)))
  }
  d <- data.frame(y = y, g = g, cov_df[ok, , drop = FALSE])
  fit <- try(stats::lm(y ~ ., data = d), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(NA_real_, NA_real_, NA_real_))
  cf <- summary(fit)$coefficients
  r <- grep("^g", rownames(cf))
  if (!length(r)) return(c(NA_real_, NA_real_, NA_real_))
  c(cf[r[1], 1], cf[r[1], 3], cf[r[1], 4])
}

#' Two-stage niche differential expression
#'
#' Estimates each patient's niche slope, then contrasts those slopes between
#' conditions. Because \code{condition} is assigned per patient, patients are
#' the experimental units; the one-stage [fitSpiDE()] path treats cells as
#' replicates of that contrast, which is anti-conservative when between-patient
#' variation in the niche-expression relationship is real.
#'
#' On a patient-label permutation of a 55-patient CosMx cohort, where every
#' triplet is null by construction, [fitSpiDE()] with \code{random = "intercept"}
#' returned ~576 false calls per replicate (realized FDR 1.00 against a nominal
#' 0.05). This returned raw type-I 0.036 and zero false calls, ~100x faster.
#'
#' \strong{It does not solve multiplicity.} On a spiked plasmode the true
#' triplets were enriched 5x at \eqn{\alpha = 0.05} and 65x at
#' \eqn{\alpha = 0.001}, yet nothing survived BH across a full-panel space of
#' ~1.8 million triplets, and the hierarchical cascade did not help (combining a
#' gene's ~137 mostly-null triplets dilutes a single true effect). Restrict
#' \code{index}, \code{niche} and the gene set to a pre-specified hypothesis:
#' roughly 4,000 tests is the order at which a \eqn{p \approx 10^{-5}} effect
#' survives.
#'
#' @param spe a SpatialExperiment with niche reducedDims (see [buildNiches()]).
#' @param condition a character, the colData column of the tested condition,
#'   constant within patient.
#' @param sigma a numeric, the bandwidth (one value; slopes are per bandwidth).
#' @param index,niche character vectors restricting index / niche cell types.
#'   \strong{Use them}: the full space is usually too large to detect anything.
#' @param patient.covariates a character vector of \emph{patient-level} colData
#'   columns to adjust for in stage 2 (e.g. sex, stage, histology). These are
#'   exactly the covariates [fitSpiDE()] rejects under \code{random != "none"}.
#' @param assay,cell_type,sample_id,name column/assay names.
#' @param min.cells minimum cells for a patient to contribute a slope (default
#'   30). Patients below it are dropped for that index type and the count is
#'   reported, rather than contributing an unstable slope.
#' @param fdr the target false discovery rate for the reported table.
#' @param verbose a logical.
#' @return a [SpiDEResults] with the tidy \code{results} table populated.
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = 30)
#' res <- twoStageSpiDE(spe, condition = "condition", sigma = 30,
#'                      min.cells = 10, verbose = FALSE)
#' head(results(res))
#' @rdname twoStageSpiDE
#' @importFrom SummarizedExperiment assay colData
#' @importFrom SingleCellExperiment reducedDim
#' @export
twoStageSpiDE <- function(spe, condition, sigma, index = NULL, niche = NULL,
                          patient.covariates = character(), assay = "counts",
                          cell_type = "cell_type", sample_id = "sample_id",
                          name = "Niche", min.cells = 30L, fdr = 0.05,
                          verbose = TRUE) {
  checkSPE(spe, assay = assay, cell_type = cell_type, sample_id = sample_id)
  checkCondition(spe, condition)
  checkNiche(spe, sigma, name = name)
  cd <- SummarizedExperiment::colData(spe)
  pat <- as.character(cd[[sample_id]])
  ct <- as.character(cd[[cell_type]])
  grp_cell <- as.character(cd[[condition]])
  # the contrast is patient-level; refuse anything else rather than silently
  # averaging a within-patient condition into nonsense
  chk <- tapply(grp_cell, pat, function(z) length(unique(z[!is.na(z)])))
  if (any(chk > 1)) {
    stop("condition '", condition, "' varies within patient; two-stage needs a ",
         "patient-level condition")
  }
  nm <- as.matrix(SingleCellExperiment::reducedDim(spe, paste0(name, sigma)))
  if (!is.null(niche)) {
    keep <- intersect(niche, colnames(nm))
    if (!length(keep)) stop("no requested niche cell types found")
    nm <- nm[, keep, drop = FALSE]
  }
  idx_types <- sort(unique(ct))
  if (!is.null(index)) idx_types <- intersect(index, idx_types)
  if (!length(idx_types)) stop("no requested index cell types found")

  Y <- SummarizedExperiment::assay(spe, assay)
  checkCounts(Y)
  lib <- Matrix::colSums(Y)
  E <- log1p(as.matrix(Y) %*% Matrix::Diagonal(x = mean(pmax(lib, 1)) / pmax(lib, 1)))
  E <- as.matrix(E)
  if (verbose) {
    message(sprintf("stage 1: %d genes x %d index types x %d niches x %d patients",
                    nrow(E), length(idx_types), ncol(nm), length(unique(pat))))
  }
  sl <- .patientSlopes(E, nm, ct, pat, idx_types, min.cells)

  pats <- sort(unique(pat))
  pgrp <- factor(vapply(pats, function(p) grp_cell[match(p, pat)], character(1)))
  cov_df <- NULL
  if (length(patient.covariates)) {
    cov_df <- as.data.frame(lapply(patient.covariates, function(cv) {
      v <- tapply(cd[[cv]], pat, function(z) z[1])
      v[pats]
    }), col.names = patient.covariates, stringsAsFactors = TRUE)
    if (verbose) message("stage 2 adjusting for: ", paste(patient.covariates, collapse = ", "))
  }

  recs <- list()
  for (ix in idx_types) {
    A <- sl[[ix]]
    if (is.null(A)) next
    npat <- sum(apply(!is.na(A[1, , , drop = FALSE]), 3, any))
    if (verbose) message(sprintf("  %-22s %d patients >= %d cells", ix, npat, min.cells))
    for (nn in colnames(nm)) {
      if (identical(nn, ix)) next          # self-interaction, as the GLM design drops
      M <- A[, nn, ]
      st <- t(apply(M, 1, .patientContrast, grp = pgrp, cov_df = cov_df))
      recs[[length(recs) + 1]] <- data.frame(
        gene = rownames(M), ct_index = ix, ct_niche = nn,
        coef = st[, 1], t = st[, 2], p.niche = st[, 3],
        stringsAsFactors = FALSE)
    }
  }
  # An empty result is a legitimate outcome, not an error: if no patient clears
  # min.cells for any index type there is nothing estimable, and the honest
  # answer is an empty table rather than an invented one.
  empty <- data.frame(gene = character(), ct_index = character(),
                      ct_niche = character(), coef = numeric(), t = numeric(),
                      p.niche = numeric(), fdr.niche = numeric(),
                      DirectionNiche = character(), bandwidth.max = numeric(),
                      stringsAsFactors = FALSE)
  out <- if (length(recs)) do.call(rbind, recs) else empty[, seq_len(6)]
  out <- out[is.finite(out$p.niche), , drop = FALSE]
  if (!nrow(out)) {
    if (verbose) message("stage 2: no estimable triplets")
    return(new("SpiDEResults", fits = list(), sigma = sigma,
               condition = condition, mode = "condition", index = idx_types,
               niche = colnames(nm), covariates = patient.covariates,
               coldata = cd, gene.weights = NULL, p.cauchy.pos = NULL,
               p.cauchy.neg = NULL, results = empty, fdr = fdr,
               call = match.call()))
  }
  out$fdr.niche <- stats::p.adjust(out$p.niche, "BH")
  out$DirectionNiche <- ifelse(out$t > 0, "Up", "Down")
  out$bandwidth.max <- sigma
  out <- out[order(out$p.niche), , drop = FALSE]
  rownames(out) <- NULL
  if (verbose) {
    message(sprintf("stage 2: %d triplets tested, %d at FDR %.2g",
                    nrow(out), sum(out$fdr.niche <= fdr), fdr))
  }
  new("SpiDEResults", fits = list(), sigma = sigma, condition = condition,
      mode = "condition", index = idx_types, niche = colnames(nm),
      covariates = patient.covariates, coldata = cd, gene.weights = NULL,
      p.cauchy.pos = NULL, p.cauchy.neg = NULL,
      results = out[out$fdr.niche <= fdr, , drop = FALSE],
      fdr = fdr, call = match.call())
}
