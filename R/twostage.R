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

#' Two-stage niche differential expression
#'
#' Estimates each core's (sample's) niche slope from a one-step response
#' anchored on a SpaNorm fit, pools those slopes per patient, then contrasts
#' the patient slopes between conditions with a variance-moderated limma fit.
#' Because \code{condition} is assigned per patient, patients are the
#' experimental units; the one-stage [fitSpiDE()] path treats cells as
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
#' @param patient a character, the colData column identifying the patient a
#'   sample (core) belongs to. Defaults to \code{NULL}, in which case each
#'   sample is its own patient (no pooling). Multiple samples per patient are
#'   pooled by precision (stage 2) before the condition contrast is taken;
#'   \code{condition} must still be constant within patient.
#' @param patient.covariates a character vector of \emph{patient-level} colData
#'   columns to adjust for in stage 2 (e.g. sex, stage, histology). These are
#'   exactly the covariates [fitSpiDE()] rejects under \code{random != "none"}.
#' @param assay,cell_type,sample_id,name column/assay names.
#' @param min.cells minimum cells for a sample to contribute a slope (default
#'   30). Samples below it are dropped for that index type and the count is
#'   reported in \code{diagnostics$inclusion}, rather than contributing an
#'   unstable slope.
#' @param fdr the target false discovery rate for the reported table.
#' @param stage1 one of "spanorm" (the default), "ols" or "nb", the working
#'   response stage-1 slopes are estimated from.
#'
#'   \strong{"spanorm"} linearises the one-step working response at a stored
#'   [SpaNorm::SpaNorm()] fit (\code{metadata(spe)$SpaNorm}; see \code{stage1}
#'   errors if none is found) and, under \code{epsilon = "addback"}, adds the
#'   fitted biology component back so that only the library-size/batch part of
#'   the fit is removed. This is the intended default: it reuses the
#'   normalisation the rest of the package already relies on and needs a
#'   single joint per-(sample, index) weighted fit, not one \code{fitNB} call
#'   per subset.
#'
#'   \strong{"ols"} regresses log-CPM directly on the niche columns with unit
#'   weights. It needs no stored SpaNorm fit and no dispersion estimate, so it
#'   is the fallback when one is not available (e.g. \code{data(toySpiDE)}).
#'
#'   \strong{"nb" remains the most expensive path.} It fits a fresh
#'   \code{SpaNorm::fitNB} per (sample, index) subset carrying
#'   \code{[1, log-library-size, niche columns]}; at panel scale (thousands of
#'   genes) this cost is per-GENE (dispersion estimation and IRLS setup), not
#'   per cell, so it is offered for small index/niche restrictions where the
#'   subset count is low rather than as a default.
#' @param epsilon one of "addback" (the default) or "residual", only used by
#'   \code{stage1 = "spanorm"}. "addback" linearises at the fitted mean and
#'   adds the biology (non-library/batch) component back to the working
#'   response, removing only the library-size/batch effect; "residual" leaves
#'   the bare working residual, which under-states the niche slope when a
#'   niche column overlaps the biology basis (see
#'   \code{diagnostics$r2}).
#' @param winsor,lambda.a,maxit.psi,backend forwarded to
#'   \code{\link[SpaNorm]{fitNB}} by the "nb" path.
#' @param verbose a logical.
#' @return a [SpiDEResults] with the tidy \code{results} table populated
#'   (unchanged schema: \code{gene}, \code{ct_index}, \code{ct_niche},
#'   \code{coef}, \code{t}, \code{p.niche}, \code{fdr.niche},
#'   \code{DirectionNiche}, \code{bandwidth.max}). \code{@fits} is empty (no
#'   per-bandwidth GLM fit exists for this estimator). Diagnostics are
#'   attached at \code{r@diagnostics}, a list of three tables: \code{r2} (the
#'   niche columns' R2 against the SpaNorm biology basis, per sample x index,
#'   "spanorm" stage1 only), \code{inclusion} (the per-index patient inclusion
#'   table, with a warning when \code{min.cells} dropout is associated with
#'   \code{condition}), and \code{tau2} (the DerSimonian-Laird between-patient
#'   variance per index x niche used to weight the stage-2 contrast).
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = 30)
#' res <- twoStageSpiDE(spe, condition = "condition", sigma = 30,
#'                      min.cells = 10, stage1 = "ols", verbose = FALSE)
#' head(results(res))
#' @rdname twoStageSpiDE
#' @importFrom SummarizedExperiment assay colData
#' @importFrom SingleCellExperiment reducedDim
#' @export
twoStageSpiDE <- function(spe, condition, sigma, index = NULL, niche = NULL,
                          patient = NULL, patient.covariates = character(),
                          assay = "counts", cell_type = "cell_type",
                          sample_id = "sample_id", name = "Niche",
                          min.cells = 30L, fdr = 0.05,
                          stage1 = c("spanorm", "ols", "nb"),
                          epsilon = c("addback", "residual"),
                          winsor = 4, lambda.a = 0, maxit.psi = 2,
                          backend = "cpu", verbose = TRUE) {
  stage1 <- match.arg(stage1); epsilon <- match.arg(epsilon)
  checkSPE(spe, assay = assay, cell_type = cell_type, sample_id = sample_id)
  checkCondition(spe, condition)
  checkNiche(spe, sigma, name = name)
  cd <- SummarizedExperiment::colData(spe)
  smp <- as.character(cd[[sample_id]])
  pat <- if (is.null(patient)) smp else as.character(cd[[patient]])
  ct <- as.character(cd[[cell_type]])
  grp_cell <- as.character(cd[[condition]])
  chk <- tapply(grp_cell, pat, function(z) length(unique(z[!is.na(z)])))
  if (any(chk > 1)) {
    stop("condition '", condition, "' varies within patient; two-stage needs ",
         "a patient-level condition")
  }
  s2p <- stats::setNames(vapply(unique(smp),
                                function(s) pat[match(s, smp)], character(1)),
                         unique(smp))
  nm <- as.matrix(SingleCellExperiment::reducedDim(spe, paste0(name, sigma)))
  if (!is.null(niche)) {
    keep <- intersect(niche, colnames(nm))
    if (!length(keep)) stop("no requested niche cell types found")
    nm <- nm[, keep, drop = FALSE]
  }
  idx_types <- sort(unique(ct))
  if (!is.null(index)) idx_types <- intersect(index, idx_types)
  if (!length(idx_types)) stop("no requested index cell types found")

  Y <- as.matrix(SummarizedExperiment::assay(spe, assay))
  if (stage1 %in% c("spanorm", "nb") && !.looksLikeCounts(Y)) {
    stop("stage1 = '", stage1, "' needs counts, but assay '", assay,
         "' does not look like counts")
  }
  comp <- if (stage1 == "spanorm") .spanormComponents(spe) else NULL
  E <- NULL
  if (stage1 == "ols") {
    lib <- colSums(Y)
    E <- log1p(sweep(Y, 2, mean(pmax(lib, 1)) / pmax(lib, 1), "*"))
  }
  if (verbose) {
    message(sprintf("stage 1 (%s): %d genes x %d index types x %d niches x %d samples",
                    stage1, nrow(Y), length(idx_types), ncol(nm),
                    length(unique(smp))))
  }
  sl <- .sampleSlopes(Y, E, comp, nm, ct, smp, idx_types, min.cells,
                      stage1 = stage1, epsilon = epsilon, winsor = winsor,
                      lambda.a = lambda.a, maxit.psi = maxit.psi,
                      backend = backend, verbose = verbose)

  pats <- unique(unname(s2p))
  pgrp <- factor(stats::setNames(
    vapply(pats, function(p) grp_cell[match(p, pat)], character(1)), pats))
  Xdes <- if (length(patient.covariates)) {
    cv <- as.data.frame(lapply(patient.covariates, function(v) {
      x <- tapply(cd[[v]], pat, function(z) z[1]); x[pats]
    }), col.names = patient.covariates, stringsAsFactors = TRUE)
    stats::model.matrix(~ g + ., data = cbind(data.frame(g = pgrp[pats]), cv))
  } else {
    stats::model.matrix(~ g, data = data.frame(g = pgrp[pats]))
  }
  diag_incl <- .inclusionDiagnostics(sl$ncells, s2p, pgrp, min.cells)

  recs <- list(); tau_tab <- list()
  for (ix in idx_types) {
    pooled <- .poolPatientSlopes(sl$beta[[ix]], sl$var[[ix]], s2p)
    for (nn in colnames(nm)) {
      if (identical(nn, ix)) next            # self-interaction, as the GLM drops
      B <- pooled$beta[, nn, pats, drop = TRUE]
      V <- pooled$var[, nn, pats, drop = TRUE]
      if (is.null(dim(B))) {
        B <- matrix(B, nrow = nrow(Y), dimnames = list(rownames(Y), pats))
        V <- matrix(V, nrow = nrow(Y), dimnames = list(rownames(Y), pats))
      }
      if (all(!is.finite(B))) next
      t2 <- .tau2DL(B, V, Xdes)
      st <- .limmaStage2(B, V, t2, Xdes)
      st <- st[is.finite(st$p), , drop = FALSE]
      if (!nrow(st)) next
      tau_tab[[length(tau_tab) + 1L]] <-
        data.frame(ct_index = ix, ct_niche = nn, tau2 = t2)
      recs[[length(recs) + 1L]] <- data.frame(
        gene = st$gene, ct_index = ix, ct_niche = nn,
        coef = st$coef, t = st$t, p.niche = st$p,
        stringsAsFactors = FALSE)
    }
  }
  empty <- data.frame(gene = character(), ct_index = character(),
                      ct_niche = character(), coef = numeric(), t = numeric(),
                      p.niche = numeric(), fdr.niche = numeric(),
                      DirectionNiche = character(), bandwidth.max = numeric(),
                      stringsAsFactors = FALSE)
  diagnostics <- list(
    r2 = sl$r2, inclusion = diag_incl,
    tau2 = if (length(tau_tab)) do.call(rbind, tau_tab) else
      data.frame(ct_index = character(), ct_niche = character(),
                 tau2 = numeric()))
  if (!length(recs)) {
    if (verbose) message("stage 2: no estimable triplets")
    return(new("SpiDEResults", fits = list(), sigma = sigma,
               condition = condition, mode = "condition", index = idx_types,
               niche = colnames(nm), covariates = patient.covariates,
               coldata = cd, gene.weights = NULL, p.cauchy.pos = NULL,
               p.cauchy.neg = NULL, results = empty, fdr = fdr,
               call = match.call(), diagnostics = diagnostics))
  }
  out <- do.call(rbind, recs)
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
      fdr = fdr, call = match.call(), diagnostics = diagnostics)
}
