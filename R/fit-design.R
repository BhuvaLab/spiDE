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
#
# Without blocks, W is a dense matrix [patient indicators | dense columns] and
# polishNB() absorbs the indicators as 1x1 blocks. With blocks (.depthBlocks():
# the within-patient library-size spline), each patient's intercept and spline
# columns form one absorbed block, held compactly as a SpaNorm::nbBlockDesign()
# ([dense columns | patient 1's q | patient 2's q | ...], each patient's block
# zero-padded to the widest; the 1e-3 ridge keeps the padding at 0), which
# costs one gram of the dense columns plus the block's q per Newton step
# whatever the number of patients. Zs then holds each patient's own, unpadded
# block over its cells, which inference residualises on (.partialBlock()).
.indexDesign <- function(L, cov, patient, trt = NULL, strata = NULL, tested = colnames(L),
                         blocks = NULL) {
  np <- nlevels(patient)
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
    X <- cbind(C, L)
    tcols <- tested
  } else {
    TN <- L * trt
    colnames(TN) <- paste0("condition:", colnames(L))
    X <- cbind(C, L, TN)
    tcols <- paste0("condition:", tested)
  }
  if (!is.null(blocks)) {
    W <- SpaNorm::nbBlockDesign(X, blocks$Zc, patient)
    px <- ncol(X); q <- ncol(blocks$Zc)
    return(list(W = W, patient = as.integer(patient), patients = levels(patient), npat = np,
                absorb = NULL, start = c(rep(FALSE, px), TRUE, rep(FALSE, q - 1L)),
                pen = c(rep(0, px), rep(1e-3, q)), dense = seq_len(px),
                intercept_cols = px + (seq_len(np) - 1L) * q + 1L, Zs = blocks$Zs,
                tested = match(tcols, colnames(W)), niche_cols = match(colnames(L), colnames(W)),
                tested_niche = tested))
  }
  P <- stats::model.matrix(~ 0 + patient)
  colnames(P) <- paste0("patient:", levels(patient))
  W <- cbind(P, X)
  list(W = W, patient = as.integer(patient), patients = levels(patient), npat = np,
       absorb = c(rep(TRUE, np), rep(FALSE, ncol(X))), start = c(rep(TRUE, np), rep(FALSE, ncol(X))),
       pen = c(rep(1e-3, np), rep(0, ncol(X))), dense = if (ncol(X)) np + seq_len(ncol(X)) else integer(),
       intercept_cols = seq_len(np), Zs = NULL,
       tested = match(tcols, colnames(W)), niche_cols = match(colnames(L), colnames(W)),
       tested_niche = tested)
}

# The linear predictor W alpha of one gene, for a dense design or a compact
# SpaNorm::nbBlockDesign().
.linPred <- function(W, alpha) {
  if (inherits(W, "nbBlockDesign")) {
    return(as.numeric(W$X %*% alpha[W$xi]) + as.numeric(W$Zsp %*% alpha[W$zi]))
  }
  as.numeric(W %*% alpha)
}

# The non-absorbed (dense) columns of the design, and each patient's absorbed
# block over its own cells (its intercept alone without depth blocks).
.denseX <- function(des) {
  if (inherits(des$W, "nbBlockDesign")) des$W$X else des$W[, des$dense, drop = FALSE]
}
.patientBlocks <- function(des, rows_of) {
  if (!is.null(des$Zs)) return(des$Zs)
  lapply(seq_along(rows_of), function(s) matrix(1, length(rows_of[[s]]), 1L))
}

# X residualised on a patient's absorbed block Z with weights w (Frisch-Waugh-
# Lovell): X - Z (Z'WZ)^-1 Z'WX. With Z the patient's intercept alone this is
# centring within the patient with weights w. NULL if the block is singular.
.partialBlock <- function(X, w, Z) {
  if (ncol(Z) == 1L && all(Z == 1)) return(sweep(X, 2, colSums(X * w) / sum(w)))
  G <- crossprod(Z * sqrt(w))
  cf <- tryCatch(solve(G, crossprod(Z, X * w)), error = function(e) NULL)
  if (is.null(cf)) return(NULL)
  X - Z %*% cf
}

# The within-patient library-size spline of depth = "spatial_spline": for each
# section of each patient, the columns [l, l * B] over that section's index
# cells, where l is log library size centred within (section, index type) and
# B the section's centred natural-spline tensor basis of position
# (SpaNorm::tpsBasis() on the whole section's coordinates: the basis SpaNorm's
# library-size function h(x, y) is built from), so each gene's depth effect is
# logLS * (a + h(x, y)) within the patient, as in SpaNorm, fitted jointly with
# -- and competing with -- the niche terms. The basis has df x df columns,
# df = 3 for a section with >= 200 index cells, 2 with >= 100, and 0 (l alone)
# below. Returns Zs (per patient: [1 | its sections' columns] over its cells,
# ascending), Zc (every cell's row of its patient's Zs, zero-padded to the
# widest patient) and r2 (patients x niches: the R^2 of each niche column on
# the patient's block, a diagnostic of how much of the niche covariate the
# spline could absorb).
.depthBlocks <- function(ell, xy_all, sec_all, ik, patient, L) {
  if (!all(is.finite(ell))) stop("depth = \"spatial_spline\" needs a finite log library size for every cell", call. = FALSE)
  sec <- sec_all[ik]
  S <- nlevels(patient); pid <- as.integer(patient)
  Zs <- vector("list", S)
  r2 <- matrix(NA_real_, S, ncol(L), dimnames = list(levels(patient), colnames(L)))
  for (s in seq_len(S)) {
    i <- which(pid == s)
    cols <- list()
    for (sg in unique(sec[i])) {
      j <- which(sec[i] == sg)
      l <- ell[i[j]] - mean(ell[i[j]])
      df <- if (length(j) >= 200L) 3L else if (length(j) >= 100L) 2L else 0L
      z <- matrix(0, length(i), 1L + if (df) df^2 else 0L)
      z[j, 1] <- l
      if (df) {
        ref <- which(sec_all == sg)
        B <- SpaNorm::tpsBasis(xy_all[ik[i[j]], 1], xy_all[ik[i[j]], 2], df = c(df, df),
                               ref.x = xy_all[ref, 1], ref.y = xy_all[ref, 2])
        z[j, -1] <- l * B
      }
      cols[[length(cols) + 1L]] <- z
    }
    Zs[[s]] <- cbind(1, do.call(cbind, cols))
    for (k in seq_len(ncol(L))) {
      y <- L[i, k]; tss <- sum((y - mean(y))^2)
      if (tss > 0) r2[s, k] <- 1 - sum(stats::lm.fit(Zs[[s]], y)$residuals^2) / tss
    }
  }
  q <- max(vapply(Zs, ncol, integer(1)))
  Zc <- matrix(0, length(ell), q, dimnames = list(NULL, c("intercept", paste0("depth", seq_len(q - 1L)))))
  for (s in seq_len(S)) Zc[which(pid == s), seq_len(ncol(Zs[[s]]))] <- Zs[[s]]
  list(Zc = Zc, Zs = Zs, r2 = r2)
}
