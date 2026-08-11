# Stage-1 machinery for the two-stage estimator: SpaNorm fit extraction, the
# one-step working response, the joint weighted slope fit, and the basis-R2
# diagnostic. See docs/superpowers/specs/2026-08-11-twostage-fixes-design.md.

#' Extract and validate the stored SpaNorm fit
#'
#' The fit supplies the biology/LS design split (`wtype`), the coefficients
#' and the dispersion that stage 1 linearises around. Kept as small pieces --
#' full genes x cells matrices are only materialised per (sample, index)
#' subset by .stage1Epsilon().
#' @param spe a SpatialExperiment previously normalised with SpaNorm().
#' @return list(alpha, gmean, psi, W, bio) with `bio` a logical over columns
#'   of `W` marking the biology block.
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
  list(alpha = fit$alpha, gmean = fit$gmean, psi = fit$psi, W = fit$W,
       bio = wt == "biology")
}

#' One-step working response and weights for a cell subset
#'
#' Linearises log counts at the FULL fitted mean (the correct expansion
#' point: handles zeros, no baseline leakage) and, under "addback", returns
#' the biology component so that, net, only the LS (and batch) effect is
#' removed. Materialises genes x subset matrices only.
#' @param Y counts matrix (genes x all cells; dense, sparse, or DelayedArray
#'   -- only the `cells` subset is ever densified).
#' @param comp the .spanormComponents() list.
#' @param cells integer/logical index of the subset's cells.
#' @param epsilon "addback" (default; eta_bio + z) or "residual" (z).
#' @return list(eps, w): the response and NB working weights, genes x cells.
#' @noRd
.stage1Epsilon <- function(Y, comp, cells, epsilon = c("addback", "residual")) {
  epsilon <- match.arg(epsilon)
  Wc <- comp$W[cells, , drop = FALSE]
  eta <- comp$gmean + tcrossprod(comp$alpha, Wc)
  mu <- exp(eta)
  # Densify only this (sample, index) subset -- Y itself may be sparse/
  # DelayedArray (twoStageSpiDE() no longer densifies the whole matrix).
  Yc <- as.matrix(Y[, cells, drop = FALSE])
  z <- (Yc - mu) / mu
  w <- mu / (1 + comp$psi * mu)
  eps <- if (epsilon == "addback") {
    comp$gmean + tcrossprod(comp$alpha[, comp$bio, drop = FALSE],
                            Wc[, comp$bio, drop = FALSE]) + z
  } else {
    z
  }
  list(eps = eps, w = w)
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
#' (sample, index) subset in the diagnostics.
#' @noRd
.nicheBasisR2 <- function(X, B) {
  f <- stats::lm.fit(cbind(1, B), X)
  res <- as.matrix(f$residuals)
  tss <- colSums(sweep(X, 2, colMeans(X))^2)
  r2 <- 1 - colSums(res^2) / pmax(tss, 1e-12)
  stats::setNames(pmin(pmax(r2, 0), 1), colnames(X))
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
                          epsilon = c("addback", "residual"),
                          winsor = 4, lambda.a = 0, maxit.psi = 2,
                          backend = "cpu", verbose = FALSE) {
  stage1 <- match.arg(stage1); epsilon <- match.arg(epsilon)
  smps <- sort(unique(smp)); niches <- colnames(nm)
  gn <- rownames(Y)
  # Matrix::colSums(), not the bare generic: base::colSums() has no method
  # for a sparse Matrix unless the Matrix package is attached (not just
  # loaded), and Matrix::colSums() also handles plain matrices and
  # DelayedArray correctly, so one call covers every Y this function accepts
  # without densifying it.
  loglib <- log(pmax(Matrix::colSums(Y), 1))
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
      if (stage1 == "spanorm") {
        se <- .stage1Epsilon(Y, comp, cells, epsilon)
        js <- .jointSlopes(se$eps, se$w, X)
        Bbio <- comp$W[cells, comp$bio, drop = FALSE]
        r2[[length(r2) + 1L]] <- data.frame(
          sample = ss, index = ix, niche = niches,
          r2 = as.numeric(.nicheBasisR2(X, Bbio)))
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
        js <- list(beta = matrix(NA_real_, nrow(Y), length(niches),
                                 dimnames = list(gn, niches)),
                   var = matrix(NA_real_, nrow(Y), length(niches),
                                dimnames = list(gn, niches)))
        wnb <- mu / (1 + f$psi * mu)
        for (g in seq_len(nrow(Y))) {
          info <- crossprod(Wn * wnb[g, ], Wn)
          vc <- diag(SpaNorm::invert_mat(info))
          js$beta[g, keep] <- f$alpha[g, -(1:2)]
          js$var[g, keep] <- vc[-(1:2)]
        }
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
                    niche = character(), r2 = numeric()),
       ncells = do.call(rbind, ncl))
}
