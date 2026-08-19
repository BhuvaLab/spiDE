# Stage-1 machinery for the two-stage estimator: SpaNorm fit extraction, the
# one-step working response, the joint weighted slope fit, and the basis-R2
# diagnostic. See design/specs/2026-08-11-twostage-fixes-design.md.

#' Extract and validate the stored SpaNorm fit
#'
#' The fit supplies the technical (ls + batch) columns whose linear predictor
#' becomes the stage-1 OFFSET, plus the biology split kept for the basis-R2
#' diagnostic. Kept as small pieces -- full genes x cells matrices are only
#' materialised per (sample, index) subset by .stage1Offset().
#'
#' Why an offset and not the old "addback" working response: SpaNorm models
#' biology only to ANCHOR the normalisation -- without a biology term the LS
#' spline would over-correct, because LS confounds biology. The fitted biology
#' component is a very smooth catch-all, not successfully modelled biology,
#' and must not be read back into a response. spiDE therefore models all
#' biology from scratch and assumes ONLY that the LS (and within-sample batch,
#' e.g. field-of-view) effects are appropriately modelled -- which is exactly
#' what fixing their linear predictor at coefficient 1 expresses.
#' @param spe a SpatialExperiment previously normalised with SpaNorm().
#' @return list(alpha, gmean, psi, W, bio, off) with `bio`/`off` logicals over
#'   columns of `W` marking the biology block and the ls+batch (offset) block.
#' @noRd
.spanormComponents <- function(spe) {
  fit <- S4Vectors::metadata(spe)$SpaNorm
  if (is.null(fit)) {
    stop("no SpaNorm fit found in metadata(spe)$SpaNorm; run ",
         "SpaNorm::SpaNorm(spe) first (stage1 = \"spanorm\" needs it), or ",
         "use stage1 = \"ols\"")
  }
  if (fit$ncells != ncol(spe) || fit$ngenes != nrow(spe)) {
    stop("the stored SpaNorm fit (", fit$ngenes, " x ", fit$ncells,
         ") does not match spe (", nrow(spe), " x ", ncol(spe), "); ",
         "re-run SpaNorm::SpaNorm() on this object")
  }
  wt <- as.character(fit$wtype)
  off <- wt %in% c("ls", "batch")
  if (!any(off)) {
    stop("the stored SpaNorm fit has no 'ls' or 'batch' columns to form the ",
         "stage-1 offset from (wtype: ", paste(unique(wt), collapse = ", "),
         ")")
  }
  list(alpha = fit$alpha, gmean = fit$gmean, psi = fit$psi, W = fit$W,
       bio = wt == "biology", off = off)
}

#' The fixed technical offset for a cell subset
#'
#' The linear predictor of the ls + batch columns only:
#' O = gmean + alpha[, off] %*% t(W[cells, off]), genes x cells. Passed to
#' SpaNorm::fitNB() as a fixed offset (coefficient 1 by construction), so it
#' cannot absorb depth-correlated biology the way a fitted covariate can.
#' gmean is included for scale; the stage-1 intercept absorbs any per-gene
#' constant either way.
#' @noRd
.stage1Offset <- function(comp, cells) {
  Wc <- comp$W[cells, comp$off, drop = FALSE]
  comp$gmean + tcrossprod(comp$alpha[, comp$off, drop = FALSE], Wc)
}

#' Joint WLS of the working response on all niche columns, per gene
#'
#' Design [1, X] (the intercept absorbs each gene's level, so no per-gene
#' weighted centring is needed). Per-gene weights preclude one shared Gram
#' matrix, so cross-products are batched as matrix products over the
#' upper-triangle column pairs, then each gene's small (k+1) system is
#' solved in a loop. Degenerate columns (weighted variance ~ 0) are dropped
#' up front and return NA. Variance is the usual WLS sandwich-free form
#' phi_g * (X'W_gX)^{-1} with phi_g the weighted RSS over n - rank.
#' @param eps,w genes x cells response and weights (Task 3).
#' @param X cells x k centred niche matrix.
#' @return list(beta, var): genes x k, NA for dropped columns.
#' @noRd
.jointSlopes <- function(eps, w, X) {
  G <- nrow(eps); n <- ncol(eps); k <- ncol(X)
  out_b <- out_v <- matrix(NA_real_, G, k,
                           dimnames = list(rownames(eps), colnames(X)))
  keep <- apply(X, 2, function(x) stats::var(x) > 1e-10)
  if (!any(keep)) return(list(beta = out_b, var = out_v))
  D <- cbind(`(Intercept)` = 1, X[, keep, drop = FALSE])
  p <- ncol(D)
  # per-gene Gram entries: (X'W_gX)_{ij} = sum_c w_gc D_ci D_cj
  ut <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  P <- D[, ut[, 1], drop = FALSE] * D[, ut[, 2], drop = FALSE]  # n x npairs
  Gm <- w %*% P                                                  # G x npairs
  Rhs <- (w * eps) %*% D                                         # G x p
  wrss_tot <- rowSums(w * eps^2)
  for (g in seq_len(G)) {
    A <- matrix(0, p, p)
    A[cbind(ut[, 1], ut[, 2])] <- Gm[g, ]
    A[cbind(ut[, 2], ut[, 1])] <- Gm[g, ]
    Ai <- tryCatch(solve(A), error = function(e) NULL)
    if (is.null(Ai)) next
    bg <- Ai %*% Rhs[g, ]
    df <- max(n - p, 1)
    phi <- max((wrss_tot[g] - sum(bg * Rhs[g, ])) / df, 1e-12)
    out_b[g, keep] <- bg[-1]
    out_v[g, keep] <- phi * diag(Ai)[-1]
  }
  list(beta = out_b, var = out_v)
}

#' R2 of each niche column on the biology basis, within a subset
#'
#' Measures the attenuation the "residual" response would suffer and the
#' smooth-trend overlap the "addback" response is exposed to. Reported per
#' (sample, index) subset in the diagnostics. A constant/degenerate column
#' (tss ~ 0) reports NA rather than flooring tss and reading off r2 ~ 1 --
#' that floor previously manufactured a perfect-fit false positive for a
#' column .jointSlopes() already drops as NA, so the two diagnostics
#' contradicted each other.
#' @noRd
.nicheBasisR2 <- function(X, B) {
  f <- stats::lm.fit(cbind(1, B), X)
  res <- as.matrix(f$residuals)
  tss <- colSums(sweep(X, 2, colMeans(X))^2)
  degenerate <- tss < 1e-12
  r2 <- rep(NA_real_, length(tss))
  r2[!degenerate] <- pmin(pmax(
    1 - colSums(res^2)[!degenerate] / tss[!degenerate], 0), 1)
  stats::setNames(r2, colnames(X))
}

#' Pooled per-sample dispersion for the spanorm stage-1 path
#'
#' One dispersion estimation per SAMPLE, shared by every (sample, index)
#' subset fit, instead of one per subset. Profiled on the real cohort
#' (1,280-cell subset, 800 genes): estimateDisp is ~half of each subset
#' fitNB, so pooling cuts the dominant cost ~an order of magnitude in call
#' count while estimating each psi from ALL the sample's cells rather than
#' one type's 30-2,000.
#'
#' The pooling design carries CELL-TYPE MEANS (plus the ls/batch offset) but
#' deliberately NOT the niche columns: unmodelled niche variation inflates
#' psi slightly, which errs conservative -- the direction wanted while the
#' offset arm measures anti-conservative. The design is nested in the subset
#' fitting design, which is the condition a supplied psi needs (see
#' SpaNorm::fitNB's psi docs).
#' @return named list of per-gene psi vectors, keyed by sample.
#' @noRd
.pooledPsi <- function(Y, comp, ct, smp, winsor, lambda.a, maxit.psi,
                       backend, verbose = FALSE) {
  out <- list()
  for (ss in sort(unique(smp))) {
    cells <- which(smp == ss)
    cts <- factor(ct[cells])
    D <- if (nlevels(cts) > 1L) {
      stats::model.matrix(~ 0 + cts)
    } else {
      matrix(1, length(cells), 1L, dimnames = list(NULL, "(Intercept)"))
    }
    O <- .stage1Offset(comp, cells)
    f <- try(SpaNorm::fitNB(Y[, cells, drop = FALSE], D, offset = O,
                            lambda.a = lambda.a, winsor = winsor,
                            maxit.psi = maxit.psi, backend = backend,
                            verbose = FALSE), silent = TRUE)
    if (inherits(f, "try-error")) {
      if (verbose) message("  pooled psi failed for ", ss,
                           "; its subsets fall back to per-subset estimation")
      next
    }
    out[[ss]] <- f$psi
  }
  out
}

#' Per-(sample, index) joint niche slopes for every stage-1 path
#'
#' Replaces .patientSlopes(). All paths run the same joint fit and return a
#' variance for every slope (stage 2 needs it): "spanorm" uses the one-step
#' response/weights, "ols" uses log-CPM with unit weights, "nb" fits the
#' per-subset NB GLM on [1, loglib, niches] and reads Fisher variances.
#' @return list(beta, var, r2, ncells) -- see Interfaces in the plan.
#' @noRd
.sampleSlopes <- function(Y, E, comp, nm, ct, smp, idx_types, min.cells,
                          stage1 = c("spanorm", "ols", "nb"),
                          winsor = 4, lambda.a = 0, maxit.psi = 2,
                          pool.psi = TRUE, backend = "cpu", verbose = FALSE) {
  stage1 <- match.arg(stage1)
  # pooled dispersion needs SpaNorm::fitNB(psi=); fall back silently on an
  # older SpaNorm so the path stays runnable (per-subset estimation is the
  # pre-pooling behaviour, not an error)
  psi_pool <- NULL
  if (stage1 == "spanorm" && isTRUE(pool.psi) &&
      "psi" %in% names(formals(SpaNorm::fitNB))) {
    psi_pool <- .pooledPsi(Y, comp, ct, smp, winsor, lambda.a, maxit.psi,
                           backend, verbose)
    if (verbose) message(sprintf("  pooled psi for %d of %d samples",
                                 length(psi_pool), length(unique(smp))))
  }
  smps <- sort(unique(smp)); niches <- colnames(nm)
  gn <- rownames(Y)
  # Matrix::colSums(), not the bare generic: base::colSums() has no method
  # for a sparse Matrix unless the Matrix package is attached (not just
  # loaded), and Matrix::colSums() also handles plain matrices and
  # DelayedArray correctly, so one call covers every Y this function accepts
  # without densifying it. Only the "nb" path needs it (a fitted covariate
  # COLUMN of its GLM design below -- a free coefficient, not an offset), so
  # it's a full pass over Y that the other two paths would otherwise pay for
  # and never use.
  loglib <- if (stage1 == "nb") log(pmax(Matrix::colSums(Y), 1)) else NULL
  beta <- var <- stats::setNames(vector("list", length(idx_types)), idx_types)
  r2 <- list(); ncl <- list()
  for (ix in idx_types) {
    A <- V <- array(NA_real_, c(nrow(Y), length(niches), length(smps)),
                    dimnames = list(gn, niches, smps))
    for (ss in smps) {
      cells <- which(ct == ix & smp == ss)
      ncl[[length(ncl) + 1L]] <- data.frame(sample = ss, index = ix,
                                            n = length(cells))
      if (length(cells) < min.cells) next
      X <- sweep(nm[cells, , drop = FALSE], 2,
                 colMeans(nm[cells, , drop = FALSE]))
      # Drop the index cell type's own niche column, as the one-stage GLM
      # design drops symmetric self-interactions: A's local density of A is
      # not a meaningful niche predictor for A cells. `A`/`V` above keep the
      # full niches dimension for reporting -- the self column is simply
      # never written, so it stays at its initial NA.
      if (ix %in% colnames(X)) {
        X <- X[, setdiff(colnames(X), ix), drop = FALSE]
      }
      if (stage1 == "spanorm") {
        # NB GLM on RAW counts with the SpaNorm ls+batch linear predictor as a
        # FIXED offset: design [1, niches], log mu = a0 + X b + O. Structurally
        # the "nb" path below, with the free loglib covariate replaced by the
        # offset -- the free coefficient is exactly the leak being closed.
        keep <- apply(X, 2, function(x) stats::var(x) > 1e-10)
        Wn <- cbind(`(Intercept)` = 1, X[, keep, drop = FALSE])
        if (qr(Wn)$rank < ncol(Wn)) next
        O <- .stage1Offset(comp, cells)
        # ONE fit per subset. A two-pass variant (refit at the returned psi)
        # was built and REVERTED on measurement: 1.70x cost on the unpooled
        # path (747s vs 439s on the toy) for no likelihood gain -- median
        # -0.41 across 6 subsets, better in only 2 of 6.
        #
        # The estimating and supplied-psi paths DO disagree here (up to 64
        # loglik across 6 toy subsets), but not because either fails to
        # converge: scored at a COMMON psi the single-pass fit wins 4 of 6, so
        # a refit does not recover anything. On these subsets the gap closes at
        # winsor = Inf (median +0.04, all within 1.2), consistent with the two
        # paths optimising different WINSORISED objectives -- `winsor` caps
        # against the current fitted mu, which supplying psi changes. That is
        # NOT the whole story though: on a well-conditioned synthetic design
        # the gap survives winsor = Inf unchanged, so conditioning (these
        # subsets are 12 niche columns over 23-40 cells) is likely also
        # involved. See research/notes/fitnb-offset-psi-disagreement.R, which
        # prints both regimes including the counter-example.
        #
        # PRACTICAL POINT: pool.psi is therefore an ESTIMATOR choice, not just
        # a speedup -- validate it on type-I/power, never on likelihood.
        args <- list(Y[, cells, drop = FALSE], Wn, offset = O,
                     lambda.a = lambda.a, winsor = winsor,
                     maxit.psi = maxit.psi, backend = backend,
                     verbose = FALSE)
        if (!is.null(psi_pool[[ss]])) args$psi <- psi_pool[[ss]]
        f <- try(do.call(SpaNorm::fitNB, args), silent = TRUE)
        if (inherits(f, "try-error")) next
        mu <- SpaNorm::calculateMu(f$gmean, f$alpha, Wn, offset = O)
        js <- list(beta = matrix(NA_real_, nrow(Y), ncol(X),
                                 dimnames = list(gn, colnames(X))),
                   var = matrix(NA_real_, nrow(Y), ncol(X),
                                dimnames = list(gn, colnames(X))))
        # Per-gene Fisher information. The weights wnb are the GENE's own, so
        # a thin or collinear subset can make `info` singular for SOME genes
        # and not others -- measured on the real cohort at sigma = 10, where
        # the full 13,348-gene run died on rcond 2.7e-20 while a 200-gene
        # subset of the SAME design completed. An unguarded inversion there
        # kills an 80-minute run over one gene, so failures leave NA (stage 2
        # already requires finite, positive variances and drops them).
        wnb <- mu / (1 + f$psi * mu)
        nsing <- 0L
        for (g in seq_len(nrow(Y))) {
          info <- crossprod(Wn * wnb[g, ], Wn)
          vc <- tryCatch(diag(SpaNorm::invert_mat(info)),
                         error = function(e) NULL)
          if (is.null(vc)) { nsing <- nsing + 1L; next }
          js$beta[g, keep] <- f$alpha[g, -1]
          js$var[g, keep] <- vc[-1]
        }
        if (nsing > 0L && verbose)
          message(sprintf("    %s/%s: %d of %d genes had a singular information matrix",
                          ix, ss, nsing, nrow(Y)))
        # Basis-R2 diagnostics for BOTH SpaNorm blocks. The biology overlap is
        # EXPECTED to be high (a real niche effect is smooth spatial variation,
        # so it lives partly in that basis) and is reported for context, not as
        # a defect. The LS overlap is the one that matters under the offset
        # construction: an over-flexible LS field could absorb the same spatial
        # variation the niche covariate carries -- the LS-side analogue of the
        # simulator finding that depth structure can swallow a planted effect.
        Bbio <- comp$W[cells, comp$bio, drop = FALSE]
        Bls <- comp$W[cells, comp$off, drop = FALSE]
        r2[[length(r2) + 1L]] <- rbind(
          data.frame(sample = ss, index = ix, niche = colnames(X),
                     basis = "biology", r2 = as.numeric(.nicheBasisR2(X, Bbio))),
          data.frame(sample = ss, index = ix, niche = colnames(X),
                     basis = "ls", r2 = as.numeric(.nicheBasisR2(X, Bls))))
      } else if (stage1 == "ols") {
        Ec <- E[, cells, drop = FALSE]
        js <- .jointSlopes(Ec, matrix(1, nrow(Ec), ncol(Ec)), X)
      } else {                                        # "nb" reference path
        keep <- apply(X, 2, function(x) stats::var(x) > 1e-10)
        Wn <- cbind(`(Intercept)` = 1,
                    loglib = loglib[cells] - mean(loglib[cells]),
                    X[, keep, drop = FALSE])
        if (qr(Wn)$rank < ncol(Wn)) next
        f <- try(SpaNorm::fitNB(Y[, cells, drop = FALSE], Wn,
                                lambda.a = lambda.a, winsor = winsor,
                                maxit.psi = maxit.psi, backend = backend,
                                verbose = FALSE), silent = TRUE)
        if (inherits(f, "try-error")) next
        mu <- SpaNorm::calculateMu(f$gmean, f$alpha, Wn)
        js <- list(beta = matrix(NA_real_, nrow(Y), ncol(X),
                                 dimnames = list(gn, colnames(X))),
                   var = matrix(NA_real_, nrow(Y), ncol(X),
                                dimnames = list(gn, colnames(X))))
        # Per-gene Fisher information. The weights wnb are the GENE's own, so
        # a thin or collinear subset can make `info` singular for SOME genes
        # and not others -- measured on the real cohort at sigma = 10, where
        # the full 13,348-gene run died on rcond 2.7e-20 while a 200-gene
        # subset of the SAME design completed. An unguarded inversion there
        # kills an 80-minute run over one gene, so failures leave NA (stage 2
        # already requires finite, positive variances and drops them).
        wnb <- mu / (1 + f$psi * mu)
        nsing <- 0L
        for (g in seq_len(nrow(Y))) {
          info <- crossprod(Wn * wnb[g, ], Wn)
          vc <- tryCatch(diag(SpaNorm::invert_mat(info)),
                         error = function(e) NULL)
          if (is.null(vc)) { nsing <- nsing + 1L; next }
          js$beta[g, keep] <- f$alpha[g, -(1:2)]
          js$var[g, keep] <- vc[-(1:2)]
        }
        if (nsing > 0L && verbose)
          message(sprintf("    %s/%s: %d of %d genes had a singular information matrix",
                          ix, ss, nsing, nrow(Y)))
      }
      A[, colnames(js$beta), ss] <- js$beta
      V[, colnames(js$var), ss] <- js$var
    }
    beta[[ix]] <- A; var[[ix]] <- V
    if (verbose) {
      message(sprintf("  %-22s %d samples >= %d cells", ix,
                      sum(apply(!is.na(A[1, , , drop = FALSE]), 3, any)),
                      min.cells))
    }
  }
  list(beta = beta, var = var,
       r2 = if (length(r2)) do.call(rbind, r2) else
         data.frame(sample = character(), index = character(),
                    niche = character(), basis = character(), r2 = numeric()),
       ncells = do.call(rbind, ncl))
}
