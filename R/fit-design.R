# The per-index-type design of spiDE's engines.
#
# Every index cell type k is fitted on its own: one negative binomial GLM per
# gene over the type-k cells of every patient, with
#   [patient intercepts | covariates | niche log-densities | (condition x niche)]
# The patient intercepts are absorbed by Schur complement in SpaNorm::polishNB()
# (unpenalised apart from a 1e-3 ridge that keeps a patient with no counts of a
# gene finite), so every niche coefficient is a WITHIN-patient slope and any
# between-patient association (testNicheAbundance()'s question) is absorbed.
# They also absorb the condition main effect: the condition only enters through
# the condition x niche columns, and only in the sandwich engine's condition
# model. The niche columns are log1p densities of every niche type except the
# index type's own (.isSelfNiche(), which honours mergeNiches() groups).

# The cells of index type k that enter the fit: type-k cells of patients
# contributing at least min.cells of them.
.indexCells <- function(ct, smp, k, min.cells) {
  ik <- which(!is.na(ct) & ct == k)
  tab <- table(smp[ik])
  keep <- names(tab)[tab >= min.cells]
  ik[smp[ik] %in% keep]
}

# Niche columns for index type k: every column that is not k's own niche, with
# positive variance over k's cells. tested = the ones the user asked to test
# (NULL = all of them); the rest stay in the design as adjustments.
.nicheColumns <- function(NM, k, ik, niche = NULL, group_map = NULL) {
  nm <- colnames(NM)
  self <- .isSelfNiche(rep(.sanitise(k), length(nm)), .sanitise(nm), group_map)
  cols <- nm[!self]
  v <- apply(log1p(NM[ik, cols, drop = FALSE]), 2, stats::var)
  cols <- cols[is.finite(v) & v > 0]
  tested <- if (is.null(niche)) cols else intersect(cols, niche)
  list(cols = cols, tested = tested)
}

# The design for one index type. L: cells x niche columns (log1p densities,
# tested columns first). cov: cells x covariates (centred), possibly 0 columns.
# patient: a factor over the index cells. trt: NULL (the pooled model) or 0/1
# per cell (the condition model: condition x niche columns). strata: NULL or a
# patient-level factor per cell, entered as strata x niche nuisance columns (its
# main effect is absorbed by the patient intercepts) -- needed where the
# condition is confounded with a slide or batch, whose niche slopes would
# otherwise load onto the condition contrast.
.indexDesign <- function(L, cov, patient, trt = NULL, strata = NULL, tested = colnames(L)) {
  P <- stats::model.matrix(~ 0 + patient)
  colnames(P) <- paste0("patient:", levels(patient))
  np <- ncol(P)
  C <- if (is.null(cov)) matrix(0, nrow(L), 0) else cov
  if (!is.null(strata)) {
    st <- droplevels(factor(strata))
    if (nlevels(st) > 1L) {
      Sd <- stats::model.matrix(~ st)[, -1, drop = FALSE]
      SN <- do.call(cbind, lapply(seq_len(ncol(L)), function(j) Sd * L[, j]))
      colnames(SN) <- paste0("strata:", rep(levels(st)[-1], ncol(L)), ":", rep(colnames(L), each = ncol(Sd)))
      C <- cbind(C, SN[, colSums(SN != 0) > 0, drop = FALSE])
    }
  }
  if (is.null(trt)) {
    W <- cbind(P, C, L)
    tcols <- tested
  } else {
    TN <- L * trt
    colnames(TN) <- paste0("condition:", colnames(L))
    W <- cbind(P, C, L, TN)
    tcols <- paste0("condition:", tested)
  }
  list(W = W, patient = as.integer(patient), patients = levels(patient), npat = np,
       absorb = c(rep(TRUE, np), rep(FALSE, ncol(W) - np)),
       pen = c(rep(1e-3, np), rep(0, ncol(W) - np)),
       dense = if (ncol(W) > np) (np + 1L):ncol(W) else integer(),
       tested = match(tcols, colnames(W)), niche_cols = match(colnames(L), colnames(W)),
       tested_niche = tested)
}
