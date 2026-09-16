# The per-gene Newton, restructured so a block of genes shares every read of
# the design.
#
# Nothing here is a new estimator. .polishGene() computes the per-gene penalised
# NB optimum and .newtonSolver() already absorbs the nested indicator block; this
# file computes the same thing for B genes at once, and .polishGene() stays in
# the tree as the reference implementation and the test oracle.
#
# What batching buys, measured at the cohort shape (n = 77,454, px = 398;
# FINDINGS.md 2026-09-15). Every gene in a block shares W, so the two
# design-sized matrix-vector products per Newton step -- the linear predictor
# and the score -- become one GEMM each:
#
#   B genes        separate matvecs   one GEMM   per gene
#   1                    0.07 s         0.08 s     79 ms
#   64                   7.38 s         0.11 s    1.7 ms
#   256                 29.33 s         0.22 s    0.8 ms
#
# The gram is NOT batched here. Each gene has its own working weights, so it has
# its own information matrix, and on the CPU that is irreducible; it stays a
# per-gene call into the solver. After the factorisation memoisation it is paid
# on roughly one step in three, and it is what the device path exists to attack.
#
# The dispersion search is also left per gene. At 0.16 s per optimize() against
# ~18 s of gram per gene it is a rounding error, and keeping stats::optimize()
# keeps this a pure restructuring: no deliberate numerical divergence to argue
# about while the control flow is being rebuilt.

# Live gene x cell matrices inside one batched Newton iteration: the counts
# slice, mu, the residual/weight matrix, the candidate coefficients' mean, the
# accepted mean, plus headroom for R's copy-on-modify. Sized by inspection, so
# it errs high.
SPIDE_POLISH_GENE_CELL_MATS <- 6

#' Genes per batched Newton, from a memory budget
#'
#' A gene block is sized to bound densification of the counts (2,000 genes); the
#' batched Newton's working set is gene x cell and must be bounded separately.
#' At the cohort's 77,454 cells a 2,000-gene batch would allocate over a
#' terabyte, so the block size cannot be the batch size.
#'
#' @param ncells cells in the design.
#' @param budget bytes available to ONE worker for the batched working set --
#'   per worker, not a total to be divided among them. \code{.covBatchSize()}'s
#'   budget is documented as a total and then claimed independently by every
#'   forked worker, which at 64 workers is a 128 GB claim in a stage already
#'   OOM-killed once at 503 GB. This one says what it means.
#' @return genes per batch, at least 1.
#' @noRd
.polishBatchSize <- function(ncells,
                             budget = getOption("spiDE.polish.mem.budget", 1e9)) {
  per_gene <- 8 * as.numeric(ncells) * SPIDE_POLISH_GENE_CELL_MATS
  max(1L, as.integer(floor(budget / per_gene)))
}

#' Converge a block of genes to their own penalised NB optima
#'
#' @param Yb counts for the block (genes x cells), dense.
#' @param W the design (cells x columns).
#' @param A0 fitNB's coefficients for the block (genes x columns).
#' @param psi0 fitNB's dispersions (length nrow(Yb), or recycled).
#' @param pen the per-column ridge penalty.
#' @param solver a \code{.newtonSolver()}.
#' @param maxit,tol,ct_cols,psi.range,psi.method,warm as in \code{.polishGene()}.
#'
#' There is no compaction knob. Every batched quantity is built from the active
#' rows, so the work already scales with the active set and there is no
#' masked-but-computed waste to repack. The copy-versus-mask question arises on
#' the device path, where the allocation is the tensor, and belongs there.
#' @return a list of per-gene results in the shape \code{.polishFit()} expects:
#'   \code{alpha} (genes x columns) and the vectors \code{psi}, \code{loglik},
#'   \code{iterations}, \code{restarted}, \code{capped}, \code{singular},
#'   \code{psi_bound}, \code{polished}.
#' @importFrom stats dnbinom optimize
#' @noRd
.polishBatch <- function(Yb, W, A0, psi0, pen, solver, maxit = 50L, tol = 1e-8,
                         ct_cols = NULL, psi.range = c(1e-3, 1e3),
                         psi.method = c("profile", "moderated"), warm = FALSE,
                         shared.factor = FALSE, nested = NULL) {
  psi.method <- match.arg(psi.method)
  B <- nrow(Yb)
  p <- ncol(W)
  psi0 <- rep_len(as.numeric(psi0), B)
  tW <- t(W)
  has_factor <- is.function(solver$factor)
  # Factorisation accounting, for Phase 2e. The per-gene policy refreshes one
  # gene's information matrix when THAT gene is stale; a single shared tensor
  # factorisation cannot, and must refresh the whole active stack whenever any
  # gene in it is stale. `sync` is what that would have cost, counted as the
  # engine runs, so the design question is answered by measurement rather than
  # by argument. One integer per refresh point; it changes no result.
  n_fac <- 0L
  n_fac_sync <- 0L
  # One factorisation for the active set instead of a list of per-gene ones:
  # the state a device can hold. It is refreshed when ANY active gene is stale,
  # which is a different path from per-gene staleness -- measured at 11% more
  # factorisations at a 128-gene batch (FINDINGS, 2026-09-16) and gated on the
  # objective, not on equality.
  solverB <- if (shared.factor) .newtonSolverBatch(W, pen, nested) else NULL

  # --- batched kernels -------------------------------------------------------
  # psi is length nrow(M): a matrix is column-major, so a per-gene vector
  # recycles down each column and reaches element (i, j) as psi[i]. That is the
  # whole reason these read as if psi were scalar.
  mu_of <- function(A) pmax(exp(A %*% tW), .MU_FLOOR)
  ll_of <- function(Y, M, ps, A) {
    rowSums(stats::dnbinom(Y, size = 1 / ps, mu = M, log = TRUE)) -
      0.5 * as.numeric((A^2) %*% pen)
  }

  # --- the damped Newton, over a set of genes --------------------------------
  # Mirrors .polishGene()'s newton() exactly, including the staleness policy,
  # the 1e-9 acceptance slack, the 1e-6 step floor and the rebuild-retry that
  # consumes an iteration. Every one of those is per gene and is carried as a
  # vector indexed by position within `rows`.
  newton <- function(rows, A, ps, maxit) {
    m <- length(rows)
    Y <- Yb[rows, , drop = FALSE]
    Ai <- A[rows, , drop = FALSE]
    pi_ <- ps[rows]
    Mu <- mu_of(Ai)
    ll <- ll_of(Y, Mu, pi_, Ai)
    it <- integer(m); conv <- logical(m); sing <- logical(m)
    stale <- integer(m); fac <- vector("list", m)
    st <- NULL; st_rows <- integer(0)
    act <- seq_len(m)

    while (length(act)) {
      it[act] <- it[act] + 1L
      Ya <- Y[act, , drop = FALSE]; Ma <- Mu[act, , drop = FALSE]
      Aa <- Ai[act, , drop = FALSE]; pa <- pi_[act]
      R <- (Ya - Ma) / (1 + pa * Ma)
      S <- R %*% W - sweep(Aa, 2L, pen, `*`)

      if (shared.factor) {
        # the stack is aligned to `act`; genes only ever LEAVE the active set,
        # so shrinkage is a subset of the stack rather than a rebuild
        if (!is.null(st) && !identical(st_rows, act)) {
          st <- .subsetState(st, match(act, st_rows))
          st_rows <- act
        }
        if (is.null(st) || any(stale[act] >= 3L)) {
          Wt <- Ma / (1 + pa * Ma)
          st <- solverB$factor(Wt)
          st_rows <- act
          stale[act] <- 0L
          n_fac <<- n_fac + length(act)
          n_fac_sync <<- n_fac_sync + length(act)
        }
        D <- solverB$solve(st, S)
        # same verdict as the per-gene path's `!all(is.finite(d))`: NA from a
        # singular slice, but Inf and NaN too
        bad <- !apply(D, 1L, function(r) all(is.finite(r)))
        D[bad, ] <- 0
      } else {
      refresh <- vapply(act, function(k) is.null(fac[[k]]), logical(1)) | stale[act] >= 3L
      if (any(refresh)) {
        kk <- act[refresh]
        Wt <- Ma[refresh, , drop = FALSE] / (1 + pa[refresh] * Ma[refresh, , drop = FALSE])
        for (j in seq_along(kk)) {
          w <- Wt[j, ]
          fac[[kk[j]]] <- if (has_factor) solver$factor(w) else w
        }
        stale[kk] <- 0L
        n_fac <<- n_fac + length(kk)
        n_fac_sync <<- n_fac_sync + length(act)
      }

      # the step is per gene: its own information, its own right-hand side
      D <- matrix(0, length(act), p)
      bad <- logical(length(act))
      for (j in seq_along(act)) {
        d <- solver$solve(fac[[act[j]]], S[j, ])
        if (is.null(d) || !all(is.finite(d))) { bad[j] <- TRUE; next }
        D[j, ] <- d
      }
      }
      if (any(bad)) {
        sing[act[bad]] <- TRUE
        act <- act[!bad]
        if (!length(act)) break
        D <- D[!bad, , drop = FALSE]; S <- S[!bad, , drop = FALSE]
        if (shared.factor) {
          st <- .subsetState(st, which(!bad))
          st_rows <- act
        }
      }

      # --- the line search, one trial round for every pending gene ----------
      na <- length(act)
      step <- rep(1, na); halv <- integer(na); ok <- logical(na)
      A1 <- Ai[act, , drop = FALSE]; M1 <- Mu[act, , drop = FALSE]
      L1 <- ll[act]
      pend <- seq_len(na)
      while (length(pend)) {
        cand <- Ai[act[pend], , drop = FALSE] + step[pend] * D[pend, , drop = FALSE]
        mu_c <- mu_of(cand)
        ll_c <- ll_of(Yb[rows[act[pend]], , drop = FALSE], mu_c, pi_[act[pend]], cand)
        ref <- ll[act[pend]]
        acc <- is.finite(ll_c) & ll_c >= ref - 1e-9 * abs(ref)
        if (any(acc)) {
          take <- pend[acc]
          A1[take, ] <- cand[acc, , drop = FALSE]
          M1[take, ] <- mu_c[acc, , drop = FALSE]
          L1[take] <- ll_c[acc]
          ok[take] <- TRUE
        }
        fail <- pend[!acc]
        step[fail] <- step[fail] / 2
        halv[fail] <- halv[fail] + 1L
        pend <- fail[step[fail] > 1e-6]
      }

      # --- what each gene does next ------------------------------------------
      # a stale matrix can give a bad direction: rebuild once and retry, which
      # costs an iteration, exactly as the per-gene loop does via `next`
      retry <- !ok & stale[act] > 0L
      if (any(retry)) {
        kk <- act[retry]
        if (shared.factor) {
          # one stack: a retry for any gene rebuilds it for all of them
          Wt <- Mu[act, , drop = FALSE] / (1 + pi_[act] * Mu[act, , drop = FALSE])
          st <- solverB$factor(Wt)
          st_rows <- act
          stale[act] <- 0L
        } else {
          Wt <- Mu[kk, , drop = FALSE] / (1 + pi_[kk] * Mu[kk, , drop = FALSE])
          for (j in seq_along(kk)) {
            w <- Wt[j, ]
            fac[[kk[j]]] <- if (has_factor) solver$factor(w) else w
          }
          stale[kk] <- 0L
        }
        n_fac <<- n_fac + length(kk)
        n_fac_sync <<- n_fac_sync + length(act)
      }
      stop_now <- !ok & !retry                      # line search exhausted

      leave_converged <- integer(0)
      if (any(ok)) {
        kk <- act[ok]
        gain <- L1[ok] - ll[kk]
        Ai[kk, ] <- A1[ok, , drop = FALSE]
        Mu[kk, ] <- M1[ok, , drop = FALSE]
        ll[kk] <- L1[ok]
        stale[kk] <- ifelse(halv[ok] > 2L, 3L, stale[kk] + 1L)
        # the convergence test reads the NEW log-likelihood, as the per-gene
        # loop does (`ll <- ll1` before `if (gain < tol * abs(ll))`)
        done <- gain < tol * abs(ll[kk])
        conv[kk[done]] <- TRUE
        leave_converged <- kk[done]
      }

      act <- setdiff(act, c(act[stop_now], leave_converged))
      act <- act[it[act] < maxit]
      if (shared.factor && length(act) && !identical(st_rows, act)) {
        st <- .subsetState(st, match(act, st_rows))
        st_rows <- act
      }
    }
    list(A = Ai, Mu = Mu, ll = ll, it = it, converged = conv, singular = sing)
  }

  # --- per-gene helpers, vectorised where they are shared --------------------
  sane_start <- function(rows) {
    A <- matrix(0, length(rows), p)
    ct <- if (is.null(ct_cols)) integer(0) else which(ct_cols)
    if (length(ct)) {
      for (j in ct) {
        cells <- W[, j] != 0
        A[, j] <- if (any(cells)) log(rowMeans(Yb[rows, cells, drop = FALSE]) + 1e-3) else 0
      }
    } else {
      A[, 1] <- log(rowMeans(Yb[rows, , drop = FALSE]) + 1e-3)
    }
    A
  }
  degenerate <- function(rows, A) {
    out <- logical(length(rows))
    fin <- apply(A, 1L, function(a) all(is.finite(a)))
    out[!fin] <- TRUE
    if (any(fin)) {
      kk <- which(fin)
      Eta <- A[kk, , drop = FALSE] %*% tW
      Yk <- Yb[rows[kk], , drop = FALSE]
      pos <- Yk > 0
      mins <- vapply(seq_along(kk), function(j) {
        if (!any(pos[j, ])) return(Inf)
        min(Eta[j, pos[j, ]])
      }, numeric(1))
      out[kk] <- mins < -10
    }
    out
  }
  psi_ml <- function(rows, Mu) {
    lo <- log(psi.range[1]); hi <- log(psi.range[2])
    est <- numeric(length(rows)); bound <- logical(length(rows))
    for (j in seq_along(rows)) {
      y <- Yb[rows[j], ]; mu <- Mu[j, ]
      o <- stats::optimize(function(lp) {
        -sum(stats::dnbinom(y, size = 1 / exp(lp), mu = mu, log = TRUE))
      }, c(lo, hi))
      bound[j] <- (o$minimum - lo) < 1e-3 * (hi - lo) ||
        (hi - o$minimum) < 1e-3 * (hi - lo)
      est[j] <- exp(o$minimum)
    }
    list(psi = est, at_bound = bound)
  }

  # --- the flow, mirroring .polishGene() ------------------------------------
  alpha <- A0; psi <- psi0
  loglik <- rep(NA_real_, B); iters <- integer(B)
  restarted <- capped <- singular <- psi_bound <- polished <- logical(B)

  if (warm) {
    fin <- apply(A0, 1L, function(a) all(is.finite(a)))
    if (any(fin)) {
      rows <- which(fin)
      f <- newton(rows, A0, psi0, maxit)
      good <- !f$singular & apply(f$A, 1L, function(a) all(is.finite(a))) & is.finite(f$ll)
      kk <- rows[good]
      alpha[kk, ] <- f$A[good, , drop = FALSE]
      loglik[kk] <- f$ll[good]
      iters[kk] <- f$it[good]
      capped[kk] <- !f$converged[good]
      polished[kk] <- TRUE
      singular[rows[!good]] <- f$singular[!good]
    }
    return(list(alpha = alpha, psi = psi, loglik = loglik, iterations = iters,
                restarted = restarted, capped = capped, singular = singular,
                psi_bound = psi_bound, polished = polished))
  }

  A <- A0
  deg <- degenerate(seq_len(B), A)
  if (any(deg)) {
    A[deg, ] <- sane_start(which(deg))
    restarted[deg] <- TRUE
  }
  f <- newton(seq_len(B), A, psi0, maxit)
  fin <- apply(f$A, 1L, function(a) all(is.finite(a)))
  mx <- apply(f$Mu, 1L, max)
  need <- !restarted & (f$singular | !fin | mx > 1e10)
  if (any(need)) {
    rows <- which(need)
    A[rows, ] <- sane_start(rows)
    restarted[rows] <- TRUE
    f2 <- newton(rows, A, psi0, maxit)
    # newton() indexes its result by POSITION within `rows`, not by gene id
    f$A[rows, ] <- f2$A
    f$Mu[rows, ] <- f2$Mu
    f$ll[rows] <- f2$ll; f$it[rows] <- f2$it
    f$converged[rows] <- f2$converged; f$singular[rows] <- f2$singular
  }
  fin <- apply(f$A, 1L, function(a) all(is.finite(a)))
  fell <- f$singular | !fin | !is.finite(f$ll)
  singular[fell] <- f$singular[fell]
  keep <- which(!fell)
  if (!length(keep)) {
    return(list(alpha = alpha, psi = psi, loglik = loglik, iterations = iters,
                restarted = restarted, capped = capped, singular = singular,
                psi_bound = psi_bound, polished = polished))
  }

  A <- f$A; Mu <- f$Mu; ll <- f$ll; it_total <- f$it; conv <- f$converged
  ps <- psi0
  if (psi.method == "profile") {
    live <- keep
    for (k in 1:2) {
      if (!length(live)) break
      pm <- psi_ml(live, Mu[live, , drop = FALSE])
      hit <- pm$at_bound
      psi_bound[live[hit]] <- TRUE
      live <- live[!hit]
      if (!length(live)) break
      ps[live] <- pm$psi[!hit]
      f2 <- newton(live, A, ps, 20L)
      A[live, ] <- f2$A
      Mu[live, ] <- f2$Mu
      ll[live] <- f2$ll
      it_total[live] <- it_total[live] + f2$it
      conv[live] <- conv[live] & f2$converged
    }
  }
  alpha[keep, ] <- A[keep, , drop = FALSE]
  psi[keep] <- ps[keep]
  loglik[keep] <- ll[keep]
  iters[keep] <- it_total[keep]
  capped[keep] <- !conv[keep]
  polished[keep] <- TRUE
  structure(list(alpha = alpha, psi = psi, loglik = loglik, iterations = iters,
                 restarted = restarted, capped = capped, singular = singular,
                 psi_bound = psi_bound, polished = polished),
            factorisations = c(pergene = n_fac, sync = n_fac_sync))
}
