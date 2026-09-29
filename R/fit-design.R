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
# otherwise load onto the condition contrast. blocks: NULL or list(Z, group),
# further per-patient absorbed columns (.depthBlocks(): the within-patient
# library-size spline), Z zero outside its patient's cells, group the patient
# (1..S) of each column. Each patient's intercept and its block columns are one
# absorbed block (SpaNorm::polishNB()'s grouped absorption); block_cols lists
# them per patient, and inference residualises on them (.partialBlock()).
.indexDesign <- function(L, cov, patient, trt = NULL, strata = NULL, tested = colnames(L),
                         blocks = NULL) {
  P <- stats::model.matrix(~ 0 + patient)
  colnames(P) <- paste0("patient:", levels(patient))
  np <- ncol(P)
  if (!is.null(blocks)) {
    P <- cbind(P, blocks$Z)
    grp <- c(seq_len(np), blocks$group)
  } else {
    grp <- seq_len(np)
  }
  nb <- ncol(P)
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
  absorb <- if (nb > np) c(grp, rep(NA_integer_, ncol(W) - nb)) else c(rep(TRUE, np), rep(FALSE, ncol(W) - np))
  list(W = W, patient = as.integer(patient), patients = levels(patient), npat = np,
       absorb = absorb, start = c(rep(TRUE, np), rep(FALSE, ncol(W) - np)),
       block_cols = split(seq_len(nb), factor(grp, levels = seq_len(np))),
       pen = c(rep(1e-3, nb), rep(0, ncol(W) - nb)),
       dense = if (ncol(W) > nb) (nb + 1L):ncol(W) else integer(),
       tested = match(tcols, colnames(W)), niche_cols = match(colnames(L), colnames(W)),
       tested_niche = tested)
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
# B the section's centred natural-spline tensor basis of position (the basis
# SpaNorm's library-size function h(x, y) is built from; .sectionBasis()), so
# each gene's depth effect is logLS * (a + h(x, y)) within the patient, as in
# SpaNorm, fitted jointly with -- and competing with -- the niche terms. The
# basis has df x df columns, df = 3 for a section with >= 200 index cells,
# 2 with >= 100, and 0 (l alone) below. Returns list(Z, group, r2) with r2 the
# patients x niches R^2 of each niche column on the patient's block (a
# diagnostic: how much of the niche covariate the spline could absorb).
.depthBlocks <- function(ell, xy_all, sec_all, ik, patient, L) {
  if (!all(is.finite(ell))) stop("depth = \"spatial_spline\" needs a finite log library size for every cell", call. = FALSE)
  sec <- sec_all[ik]
  S <- nlevels(patient); pid <- as.integer(patient)
  cols <- list(); grp <- integer(); nm <- character()
  r2 <- matrix(NA_real_, S, ncol(L), dimnames = list(levels(patient), colnames(L)))
  for (s in seq_len(S)) {
    Zs <- list()
    for (sg in unique(sec[pid == s])) {
      i <- which(pid == s & sec == sg)
      n <- length(i)
      l <- ell[i] - mean(ell[i])
      df <- if (n >= 200L) 3L else if (n >= 100L) 2L else 0L
      z <- matrix(0, length(ell), 1L + if (df) df^2 else 0L)
      z[i, 1] <- l
      if (df) {
        ref <- which(sec_all == sg)
        B <- .sectionBasis(xy_all[ik[i], 1], xy_all[ik[i], 2], df, xy_all[ref, 1], xy_all[ref, 2])
        z[i, -1] <- l * B
      }
      colnames(z) <- paste0("depth:", levels(patient)[s], ":", sg, ":", c("l", if (df) paste0("lB", seq_len(df^2))))
      Zs[[length(Zs) + 1L]] <- z
    }
    Zs <- do.call(cbind, Zs)
    cols[[s]] <- Zs; grp <- c(grp, rep(s, ncol(Zs)))
    i <- which(pid == s)
    Zi <- cbind(1, Zs[i, , drop = FALSE])
    for (j in seq_len(ncol(L))) {
      y <- L[i, j]; tss <- sum((y - mean(y))^2)
      if (tss > 0) r2[s, j] <- 1 - sum(stats::lm.fit(Zi, y)$residuals^2) / tss
    }
  }
  list(Z = do.call(cbind, cols), group = grp, r2 = r2)
}

# The centred natural-spline tensor basis of position over a section: ns(x, df)
# x ns(y, df), knots at the quantiles of the whole section's coordinates
# (ref.x, ref.y; all cell types) and columns centred over the section. Equal to
# SpaNorm::tpsBasis(x, y, df = c(df, df), ref.x, ref.y) (SpaNorm >= 1.7.15;
# tests/testthat/test-depth.R checks it when available).
.sectionBasis <- function(x, y, df, ref.x, ref.y) {
  bx <- splines::ns(ref.x, df = df); by <- splines::ns(ref.y, df = df)
  ev <- function(b, v) stats::predict(b, v)
  tens <- function(X, Y) do.call(cbind, lapply(seq_len(df), function(a) X[, a] * Y))
  Bref <- tens(bx, by)
  B <- tens(ev(bx, x), ev(by, y))
  sweep(B, 2, colMeans(Bref))
}
