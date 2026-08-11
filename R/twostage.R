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

#' Does an assay look like counts?
#'
#' Counts and Pearson residuals are told apart by their values, not by their
#' name, because passing the wrong assay is silent and expensive: running the
#' variance-stabilisation twice, or not at all, changes power without changing
#' anything visible in the output.
#' @noRd
.looksLikeCounts <- function(Y, n = 2000L) {
  idx <- if (length(Y) > n) sample(length(Y), n) else seq_along(Y)
  v <- as.numeric(Y[idx])
  v <- v[is.finite(v)]
  if (!length(v)) return(FALSE)
  all(v >= 0) && all(abs(v - round(v)) < 1e-8)
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
#' @param stage1 one of "ols" (the default) or "nb".
#'
#'   \strong{"nb" is currently impractical at panel scale and the default is
#'   "ols" on measurement, not preference.} A single \code{fitNB} call on one
#'   (patient, index) subset of the YTMA cohort costs \strong{23 s} even at
#'   \code{maxit.psi = 1} -- and for a 158-cell subset, so the expense is
#'   per-GENE (dispersion estimation and IRLS setup over 13,348 genes), not per
#'   cell. Across 1,320 (patient, index) combinations that is \strong{8.5
#'   hours}, four times the one-stage fit this estimator exists to replace.
#'
#'   "nb" is also not obviously better where it has been checked: on the toy
#'   fixture it ranked the planted effect 17th where "ols" ranked it 2nd, though
#'   that fixture is unfair to it (27-cell subsets, dispersion moderated over 40
#'   genes).
#'
#'   The route to making "nb" viable is FEWER CALLS, not faster ones: one fit per
#'   patient carrying \code{CellType:niche} interactions (55 calls) rather than
#'   one per (patient, index) (1,320), which the cost profile above predicts
#'   would land near the one-stage runtime. Until that is built, "nb" is offered
#'   for small index/niche restrictions where the call count is low.
#' @param winsor,lambda.a,maxit.psi,backend forwarded to
#'   \code{\link[SpaNorm]{fitNB}} by the "nb" path.
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
                          stage1 = c("auto", "nbresid", "ols", "analytic", "nb"),
                          winsor = 4, lambda.a = 0,
                          maxit.psi = 2, backend = "cpu", verbose = TRUE) {
  stage1 <- match.arg(stage1)
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
  # AUTO: counts in -> variance-stabilise them here; residuals in -> trust them.
  # Pearson residuals are the natural input (a SpaNorm normalisation run already
  # produces them via adj.method = "pearson"), and running OLS on residuals that
  # are ALREADY variance-stabilised is the correct thing -- re-standardising
  # would be wrong. Counts are detected rather than declared, because passing the
  # wrong assay silently is the expensive mistake here.
  is_counts <- .looksLikeCounts(Y)
  if (identical(stage1, "auto")) {
    stage1 <- if (is_counts) "nbresid" else "ols"
    if (verbose) {
      message(sprintf("stage1 = \"%s\" (assay '%s' %s)", stage1, assay,
                      if (is_counts) "looks like counts" else
                        "is not counts -- assumed pre-stabilised, e.g. Pearson residuals"))
    }
  }
  if (stage1 %in% c("nbresid", "nb", "analytic") && !is_counts) {
    stop("stage1 = '", stage1, "' needs counts, but assay '", assay,
         "' does not look like counts")
  }
  if (is_counts) checkCounts(Y)
  # log-CPM is only needed by the "ols" path, and densifying 13k x 77k costs
  # 8 GB -- do not pay for it when fitting NB.
  E <- NULL
  if (stage1 == "nbresid") {
    E <- .nbPearsonResiduals(Y, winsor = winsor, lambda.a = lambda.a,
                             maxit.psi = maxit.psi, backend = backend,
                             verbose = verbose)
  } else if (stage1 == "analytic") {
    E <- .pearsonResiduals(Y)
  } else if (stage1 == "ols") {
    lib <- Matrix::colSums(Y)
    E <- as.matrix(log1p(as.matrix(Y) %*%
                           Matrix::Diagonal(x = mean(pmax(lib, 1)) / pmax(lib, 1))))
  }
  if (verbose) {
    message(sprintf("stage 1: %d genes x %d index types x %d niches x %d patients",
                    nrow(E), length(idx_types), ncol(nm), length(unique(pat))))
  }
  sl <- .patientSlopes(Y, E, nm, ct, pat, idx_types, min.cells,
                       stage1 = stage1, winsor = winsor, lambda.a = lambda.a,
                       maxit.psi = maxit.psi, backend = backend,
                       verbose = verbose)

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
