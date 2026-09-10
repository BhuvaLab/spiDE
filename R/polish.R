# Per-gene convergence of the penalised NB fit.
#
# SpaNorm::fitNB fits every gene in one IRLS loop: it shares a single
# gene-averaged cell weight vector across genes, decides step-halving and
# convergence on the AGGREGATE log-likelihood, and clamps coefficient columns
# across genes. That is what makes a 13,000-gene fit affordable, and for the
# great majority of genes it is indistinguishable from the per-gene optimum. For
# a bright, cell-type-restricted gene -- whose own working weights look nothing
# like the average -- the aggregate criterion is met long before that gene's own
# score is zero: measured on the YTMA cohort, every gene in the top 5% by
# expression sat 1-4 production standard errors from its own penalised-NB
# optimum, with log-likelihood gaps of 1e4-1e6 and an edgeR dispersion ~1.6x too
# large.
#
# This stage removes the fit from the inference question. It runs AFTER fitNB
# (so the cross-gene dispersion moderation still happens on the whole gene set,
# which is why genes must not be blocked at fit time), and because a polished
# gene depends only on its own counts it is blockable and parallel -- the same
# split .blockedInference() uses.

#' Penalised NB log-likelihood of one gene
#'
#' @param y counts (length ncells).
#' @param mu the fitted mean (length ncells).
#' @param psi the NB dispersion (scalar).
#' @param a the coefficients (length ncol(W)).
#' @param pen the per-column ridge penalty.
#' @return a numeric scalar.
#' @importFrom stats dnbinom
#' @noRd
# The fitted mean is floored here and at inference, at the same value. A
# linear predictor below about -745 underflows exp() to exactly 0, and the
# negative binomial quantities built from it then divide by zero: the Pearson
# working dispersion (y - mu)^2 / (mu + psi mu^2) becomes 0/0 = NaN, which
# propagates to the standard error, the statistic and every gene-set test
# downstream. exp(-30) is far below any mean the model can meaningfully
# estimate, so flooring there is a numerical guard, not a statistical clamp --
# and applying the SAME floor in both stages keeps the polish and the inference
# on one mean function, which is the point of not winsorising here.
.MU_FLOOR <- exp(-30)

.nbPenLoglik <- function(y, mu, psi, a, pen) {
  sum(stats::dnbinom(y, size = 1 / psi, mu = mu, log = TRUE)) -
    0.5 * sum(pen * a^2)
}

#' Newton solver for the penalised NB information, with indicator columns absorbed
#'
#' The nested (sample x cell type) random intercepts are 0/1 indicators that
#' partition the cells, so their blocks of the information matrix are cheap:
#' with per-cell weights \code{w} and the design split into dense \code{X} and
#' indicators \code{Z}, \code{A = X' diag(w) X + diag(pen_x)},
#' \code{B = X' diag(w) Z} is one \code{rowsum} pass, and
#' \code{C = Z' diag(w) Z + diag(pen_z)} is DIAGONAL. The step then comes from
#' the Schur complement \code{S = A - B C^-1 B'}, whose cost is one
#' \code{ncol(X)}-column gram regardless of how many groups there are -- the
#' difference between a 345-column and a 1,000-column gram per iteration on a
#' real cohort. \code{S^-1} is also the X-block of the full penalised
#' covariance, which is what the tested columns need.
#'
#' @param W the full design (cells x columns).
#' @param pen the per-column ridge penalty.
#' @param nested a logical over the columns of \code{W} marking the indicator
#'   block; all-FALSE (or \code{NULL}) selects the plain dense path.
#' @return a list with \code{solve(w, s)} (the Newton step over all columns) and
#'   \code{xcov(w)} (the X-block of the penalised covariance). Either returns
#'   \code{NULL} on a singular system.
#' @noRd
.newtonSolver <- function(W, pen, nested = NULL) {
  if (is.null(nested)) nested <- rep(FALSE, ncol(W))
  if (!any(nested)) {
    return(list(
      solve = function(w, s) {
        info <- crossprod(W * sqrt(w))
        diag(info) <- diag(info) + pen
        tryCatch(solve(info, s), error = function(e) NULL)
      },
      xcov = function(w) {
        info <- crossprod(W * sqrt(w))
        diag(info) <- diag(info) + pen
        tryCatch(solve(info), error = function(e) NULL)
      }
    ))
  }

  xi <- which(!nested)
  zi <- which(nested)
  X <- W[, xi, drop = FALSE]
  pen_x <- pen[xi]
  pen_z <- pen[zi]
  # The indicator each cell belongs to. The columns are 0/1 and partition the
  # cells, so this dot product recovers the group index -- via round(), not
  # as.integer(): a floating-point product of 7 can come back as 6.9999999,
  # which as.integer() would truncate to the WRONG group, silently.
  Zblk <- W[, zi, drop = FALSE]
  # the absorption is exact only if C = Z' diag(w) Z is diagonal, i.e. only if
  # every cell belongs to exactly one group. Check that before relying on it.
  rs <- rowSums(Zblk)
  if (anyNA(rs) || max(abs(rs - 1)) > 1e-8) {
    stop("the nested random-effect columns are not 0/1 indicators partitioning ",
         "the cells; .newtonSolver() cannot absorb them", call. = FALSE)
  }
  gidx <- round(as.numeric(Zblk %*% seq_along(zi)))
  gf <- factor(gidx, levels = seq_along(zi))

  parts <- function(w) {
    # ONE weighted copy of X, reused for both the gram and the group sums.
    # Writing it as crossprod(X * sqrt(w)) plus rowsum(cbind(w, X * w), ...)
    # allocates two full n x ncol(X) temporaries, and at realistic sizes the
    # memory traffic -- not the flop count -- is what this stage costs.
    Xw <- X * w
    A <- crossprod(X, Xw)
    diag(A) <- diag(A) + pen_x
    cvec <- as.numeric(rowsum(w, group = gf, reorder = TRUE)) + pen_z
    B <- t(rowsum(Xw, group = gf, reorder = TRUE))   # ncol(X) x G
    list(S = A - B %*% (t(B) / cvec), B = B, cvec = cvec)
  }

  list(
    solve = function(w, s) {
      p <- parts(w)
      rhs <- s[xi] - as.numeric(p$B %*% (s[zi] / p$cvec))
      dx <- tryCatch(solve(p$S, rhs), error = function(e) NULL)
      if (is.null(dx)) return(NULL)
      dz <- (s[zi] - as.numeric(crossprod(p$B, dx))) / p$cvec
      out <- numeric(ncol(W))
      out[xi] <- dx
      out[zi] <- dz
      out
    },
    xcov = function(w) {
      tryCatch(solve(parts(w)$S), error = function(e) NULL)
    }
  )
}

#' Converge one gene to its own penalised NB optimum
#'
#' Damped Newton on the penalised log-likelihood at fixed \code{psi}, then
#' profile ML for \code{psi} at the converged mean, then a re-polish -- twice.
#' The information matrix is reused for up to three consecutive steps (only the
#' score is recomputed), since it changes slowly and each rebuild is essentially
#' the whole cost of an iteration.
#'
#' Starting from fitNB's coefficients is right for almost every gene, but for a
#' degenerate fit (a fitted mean below exp(-10) somewhere) Newton from that
#' point DIVERGES -- measured: fitted log-means reaching +57 to +358 and
#' predicted one-step gains of 1e6-1e8 against actual gains of 1e2-1e5. Those
#' genes restart from a sane point instead: the per-cell-type log mean, every
#' other coefficient zero, which converges in 5-29 iterations.
#'
#' @param y counts for this gene (length ncells).
#' @param W the design.
#' @param a0,psi0 fitNB's coefficients and dispersion for this gene.
#' @param pen the per-column ridge penalty.
#' @param solver a \code{.newtonSolver()} for this \code{W} and \code{pen}.
#' @param maxit,tol iteration cap and relative log-likelihood tolerance.
#' @param ct_cols a logical over the columns of \code{W} marking the cell-type
#'   intercepts (from the design's covtype tags).
#' @param psi.range the search interval for the profile-ML dispersion.
#' @param warm logical; \code{a0}/\code{psi0} are a CONVERGED fit at a nearby
#'   penalty (the re-polish after a variance-component step). A warm polish is
#'   a few damped Newton steps at the held dispersion: no profile-psi search
#'   (its optimum moves at second order in the penalty change) and no
#'   log-mean restart check -- a converged fit legitimately has fitted
#'   log-means below -10 where a gene is absent from a cell type, and treating
#'   that like fitNB's degenerate output threw ~100 of 769 genes back to the
#'   sane start on every re-polish pass at bandwidth 10 on the cohort (3-5 in
#'   the cold pass), each redoing a full cold polish.
#' @return a list with \code{alpha}, \code{psi}, \code{loglik},
#'   \code{iterations}, \code{restarted}, \code{capped}, \code{singular},
#'   \code{psi_bound} and \code{polished}.
#' @importFrom stats optimize
#' @noRd
.polishGene <- function(y, W, a0, psi0, pen, solver, maxit = 50L, tol = 1e-8,
                        ct_cols = NULL, psi.range = c(1e-3, 1e3),
                        psi.method = c("profile", "moderated"), warm = FALSE) {
  psi.method <- match.arg(psi.method)
  restarted <- FALSE
  singular <- FALSE

  # the sane start: cell-type (or, absent a cell-type block, overall) log means.
  # `ct_cols` comes from the design's own covtype tags. It used to be recovered
  # by a regex on colnames(W), which is a second, weaker parser of a convention
  # .tagCovtype() already owns: a user covariate literally named "CellTypeScore"
  # matched it and was assigned a log mean as though it were an indicator, and a
  # cell-type label containing ":" did not match at all.
  sane_start <- function() {
    a <- numeric(ncol(W))
    ct <- if (is.null(ct_cols)) integer(0) else which(ct_cols)
    if (length(ct)) {
      for (j in ct) {
        cells <- W[, j] != 0
        a[j] <- if (any(cells)) log(mean(y[cells]) + 1e-3) else 0
      }
    } else {
      a[1] <- log(mean(y) + 1e-3)
    }
    a
  }

  newton <- function(a, psi, maxit) {
    mu <- pmax(as.numeric(exp(W %*% a)), .MU_FLOOR)
    ll <- .nbPenLoglik(y, mu, psi, a, pen)
    it <- 0L
    converged <- FALSE
    stale <- 0L
    w <- NULL
    while (it < maxit) {
      it <- it + 1L
      s <- as.numeric(crossprod(W, (y - mu) / (1 + psi * mu))) - pen * a
      if (is.null(w) || stale >= 3L) {
        w <- mu / (1 + psi * mu)
        stale <- 0L
      }
      d <- solver$solve(w, s)
      # solve() only ERRORS below rcond ~1e-7; between that and well-conditioned
      # it returns a finite but numerically meaningless answer, which the line
      # search can accept because a badly scaled step in roughly the right
      # direction still raises the objective. A NaN/Inf right-hand side does not
      # error either. Treat both as singular rather than letting a wrong number
      # through as a converged coefficient.
      if (is.null(d) || !all(is.finite(d))) {
        singular <<- TRUE
        break
      }
      step <- 1
      ok <- FALSE
      halvings <- 0L
      while (step > 1e-6) {
        a1 <- a + step * d
        mu1 <- pmax(as.numeric(exp(W %*% a1)), .MU_FLOOR)
        ll1 <- .nbPenLoglik(y, mu1, psi, a1, pen)
        if (is.finite(ll1) && ll1 >= ll - 1e-9 * abs(ll)) {
          ok <- TRUE
          break
        }
        step <- step / 2
        halvings <- halvings + 1L
      }
      if (!ok) {
        # a stale information matrix can give a bad direction; rebuild it once
        # before giving up
        if (stale > 0L) {
          w <- mu / (1 + psi * mu)
          stale <- 0L
          next
        }
        break
      }
      # a hard line search or a stale matrix both call for a rebuild next step
      stale <- if (halvings > 2L) 3L else stale + 1L
      gain <- ll1 - ll
      a <- a1
      mu <- mu1
      ll <- ll1
      if (gain < tol * abs(ll)) {
        converged <- TRUE
        break
      }
    }
    list(a = a, mu = mu, ll = ll, it = it, converged = converged)
  }

  # Profile ML for psi at a fixed mean. The optimiser needs a bounded interval,
  # so an under-dispersed or near-empty gene lands ON a bound -- measured: a
  # Poisson gene returns 0.0046 and an all-zero gene 976, neither of which is an
  # estimate. Storing a bound as though it were one is worse than not polishing
  # the dispersion at all, because on the fixed-effects path psi scales the
  # standard error directly. Report it instead, and let the caller keep fitNB's
  # moderated value.
  psi_ml <- function(mu) {
    lo <- log(psi.range[1]); hi <- log(psi.range[2])
    o <- stats::optimize(function(lp) {
      -sum(stats::dnbinom(y, size = 1 / exp(lp), mu = mu, log = TRUE))
    }, c(lo, hi))
    edge <- (o$minimum - lo) < 1e-3 * (hi - lo) ||
      (hi - o$minimum) < 1e-3 * (hi - lo)
    list(psi = exp(o$minimum), at_bound = edge)
  }

  fallback <- function(why) {
    # Hand back fitNB's own estimate, NOT the sane start. Returning the sane
    # start would replace a usable fit with cell-type log means and exact zeros
    # on every tested coefficient, which inference then reports as t = 0 with a
    # finite SE -- a confident null for a gene that was never converged.
    list(alpha = a0, psi = psi0, loglik = NA_real_, iterations = 0L,
         restarted = FALSE, capped = FALSE, singular = identical(why, "singular"),
         psi_bound = FALSE, polished = FALSE)
  }
  a <- a0
  psi <- psi0
  if (warm) {
    if (!all(is.finite(a))) return(fallback("nonfinite"))
    f <- newton(a, psi, maxit)
    if (singular || !all(is.finite(f$a)) || !is.finite(f$ll)) {
      return(fallback(if (singular) "singular" else "nonfinite"))
    }
    return(list(alpha = f$a, psi = psi, loglik = f$ll, iterations = f$it,
                restarted = FALSE, capped = !f$converged, singular = FALSE,
                psi_bound = FALSE, polished = TRUE))
  }
  if (!all(is.finite(a)) || min(as.numeric(W %*% a)) < -10) {
    a <- sane_start()
    restarted <- TRUE
  }
  f <- newton(a, psi, maxit)
  if (!restarted && (singular || !all(is.finite(f$a)) || max(f$mu) > 1e10)) {
    singular <- FALSE
    restarted <- TRUE
    f <- newton(sane_start(), psi, maxit)
  }
  # both starts failed: keep fitNB's fit rather than an unconverged guess
  if (singular || !all(is.finite(f$a)) || !is.finite(f$ll)) {
    return(fallback(if (singular) "singular" else "nonfinite"))
  }
  converged <- f$converged
  it_total <- f$it
  psi_bound <- FALSE
  # "moderated" keeps fitNB's cross-gene moderated dispersion and converges
  # only the mean under it; "profile" re-estimates psi by profile ML at the
  # converged mean and re-polishes, twice
  if (psi.method == "profile") for (k in 1:2) {
    pm <- psi_ml(f$mu)
    if (pm$at_bound) {
      # the dispersion is not identified for this gene; keep the moderated one
      psi_bound <- TRUE
      break
    }
    psi <- pm$psi
    f2 <- newton(f$a, psi, 20L)
    it_total <- it_total + f2$it
    converged <- converged && f2$converged
    f <- f2
  }
  list(alpha = f$a, psi = psi, loglik = f$ll, iterations = it_total,
       restarted = restarted, capped = !converged, singular = singular,
       psi_bound = psi_bound, polished = TRUE)
}

#' Run this worker's BLAS single-threaded
#'
#' Called inside a \code{bplapply()} worker of a multi-worker BPPARAM. A
#' forked worker inherits the parent's OpenBLAS thread count, so N workers on
#' N cores run N x threads BLAS threads: measured at the cohort's design shape
#' (77,454 cells, 345 dense + 660 nested columns), 4 workers x 4 threads on 4
#' cores take 2.3-2.6 s per Newton step against 0.24-0.28 s at one thread each,
#' while one worker gains only 1.5x from 4 threads -- the per-gene gram is
#' memory-bound. The 0.99.19 cohort runs (4 workers x 8 threads on 8 cores)
#' spent 245-395 min on the cold polish pass this way. A no-op without
#' RhpcBLASctl (Suggests), and never called in the parent process.
#' @return the previous thread count, invisibly (\code{NA} without RhpcBLASctl).
#' @noRd
.workerBLAS <- function() {
  if (!requireNamespace("RhpcBLASctl", quietly = TRUE)) return(invisible(NA_integer_))
  prev <- RhpcBLASctl::blas_get_num_procs()
  if (isTRUE(prev > 1L)) RhpcBLASctl::blas_set_num_threads(1L)
  invisible(prev)
}

#' Converge every gene's fit, blocked over genes
#'
#' @param Y counts (genes x cells).
#' @param W the design.
#' @param alpha,psi fitNB's estimates.
#' @param pen the per-column ridge penalty (scalar or length ncol(W)).
#' @param re_group the per-column random-effect group, used only to locate the
#'   indicator block; \code{NULL} for a fixed-effects fit.
#' @param covtype the design's per-column covariate tags, used to locate the
#'   cell-type intercepts for the restart. \code{NULL} falls back to the
#'   first column.
#' @param maxit,tol forwarded to \code{.polishGene()}.
#' @param block.size,BPPARAM gene blocking and dispatch, as in [testSpiDE()].
#'   With more than one worker, each worker's BLAS is set single-threaded
#'   (\code{.workerBLAS()}) when RhpcBLASctl is installed: forked workers
#'   inherit the parent's OpenBLAS thread count, and \code{workers x threads}
#'   on \code{workers} cores was measured at 9x per Newton step at the
#'   cohort's design shape (4 x 4 on 4 cores: 2.3-2.6 s against 0.24-0.28 s).
#' @param verbose logical.
#' @param warm logical, forwarded to \code{.polishGene()}: a re-polish of a
#'   converged fit at a nearby penalty.
#' @return a list with \code{alpha}, \code{psi}, \code{loglik} and a per-gene
#'   \code{polish} data.frame.
#' @importFrom BiocParallel bplapply SerialParam bpnworkers
#' @noRd
.polishFit <- function(Y, W, alpha, psi, pen, re_group = NULL, covtype = NULL,
                       maxit = 50L, tol = 1e-8, block.size = NULL,
                       BPPARAM = BiocParallel::SerialParam(), verbose = FALSE,
                       psi.method = "profile", warm = FALSE) {
  ng <- nrow(alpha)
  if (!length(pen) %in% c(1L, ncol(W))) {
    stop("'lambda.a' must be a single value or one per design column (",
         ncol(W), " here, ", length(pen), " supplied). Note re.celltype = TRUE ",
         "adds one column per non-empty (sample, cell type), so a vector sized ",
         "for an earlier design is now too short.", call. = FALSE)
  }
  if (length(pen) == 1L) pen <- rep(pen, ncol(W))
  if (length(psi) == 1L) psi <- rep(psi, ng)
  # The negative binomial likelihood is defined on counts. On a non-integer
  # assay dnbinom() returns -Inf for every cell, so the line search rejects
  # every step, alpha is left exactly at fitNB's value and the dispersion
  # optimiser -- maximising a constant -Inf -- returns its upper bound for EVERY
  # gene. Measured: psi 999.96 across the board, with `capped` and `singular`
  # both reporting success. That is silent, total corruption, and the documented
  # real-cohort object (counts <- 2^logcounts - 1) is exactly such an assay, so
  # refuse it here rather than let it through.
  chk <- as.numeric(Y[seq_len(min(nrow(Y), 20L)), , drop = FALSE])
  chk <- chk[is.finite(chk)]
  if (length(chk) && max(abs(chk - round(chk))) > 1e-8) {
    stop("the polish stage needs integer counts: the negative binomial ",
         "likelihood is undefined otherwise, and every gene's dispersion would ",
         "silently collapse to its upper bound.\n  The assay passed is not ",
         "integer-valued (e.g. a back-transform such as 2^logcounts - 1).\n  ",
         "Use the raw counts, or skip polishSpiDE().", call. = FALSE)
  }
  ct_cols <- if (is.null(covtype)) NULL else as.character(covtype) == "CellType"
  nested <- if (is.null(re_group)) {
    rep(FALSE, ncol(W))
  } else {
    !is.na(re_group) & re_group == "SampleCellTypeInt"
  }
  solver <- .newtonSolver(W, pen, nested)

  # .chunkGenes(ng, NULL) is ONE block, which would hand every gene to a single
  # worker however many BPPARAM has -- a silent loss of the parallelism the
  # caller asked for. Absent an explicit block.size, split at least one block
  # per worker (this stage is exact per gene, so blocking never changes the
  # answer -- test-polish.R asserts that).
  # .chunkGenes(ng, NULL) is ONE block. That would (a) hand every gene to a
  # single worker however many BPPARAM has, and (b) -- worse -- densify the
  # WHOLE counts matrix at line `as.matrix(Y[gi, ])`: 13,348 x 77,454 doubles is
  # 8.3 GB, inside fitSpiDE(), under its own default SerialParam(). The
  # architecture's invariant is that the counts matrix is never eagerly
  # densified, so cap the block regardless of worker count, the way
  # .blockLoglik() does with its 2000-gene default.
  nw <- max(1L, BiocParallel::bpnworkers(BPPARAM))
  if (is.null(block.size)) {
    block.size <- max(1L, min(2000L, ceiling(ng / nw)))
  }
  blocks <- .chunkGenes(ng, block.size)
  single_blas <- .singleBLAS(BPPARAM)
  if (verbose) {
    message(sprintf("  %s %d genes per gene (%d block%s%s)",
                    if (warm) "re-polishing" else "converging", ng,
                    length(blocks), if (length(blocks) == 1L) "" else "s",
                    if (single_blas) ", one BLAS thread per worker" else ""))
  }
  # A whole-transcriptome polish is hours of work, so report progress rather
  # than going silent after the opening message. Blocks are timed as they
  # finish; under a parallel BPPARAM they complete out of order, so the count
  # is of blocks retired, not a position in the gene list.
  t0 <- Sys.time()
  nb <- length(blocks)
  step <- max(1L, nb %/% 20L)
  res <- .bplapplySingleBLAS(seq_along(blocks), function(b) {
    gi <- blocks[[b]]
    # densify the whole block once: Y may be sparse or a DelayedArray, where a
    # per-gene read costs a round trip each time (the invariant is that the
    # WHOLE matrix is never densified, not that a block is never densified --
    # .blockedInference() does exactly the same).
    Yb <- as.matrix(Y[gi, , drop = FALSE])
    out <- lapply(seq_along(gi), function(i) {
      g <- gi[[i]]
      .polishGene(as.numeric(Yb[i, ]), W, alpha[g, ], psi[[g]], pen, solver,
                  maxit = maxit, tol = tol, ct_cols = ct_cols,
                  psi.method = psi.method, warm = warm)
    })
    if (verbose && (b %% step == 0L || b == nb)) {
      message(sprintf("    block %d/%d (%.1f min elapsed)", b, nb,
                      as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    }
    out
  }, BPPARAM = BPPARAM)
  res <- unlist(res, recursive = FALSE)
  if (verbose) {
    message(sprintf(
      "  polished %d/%d genes (%d restarted, %d not converged, %d singular, %d dispersion at a bound)",
      sum(vapply(res, `[[`, logical(1), "polished")), ng,
      sum(vapply(res, `[[`, logical(1), "restarted")),
      sum(vapply(res, `[[`, logical(1), "capped")),
      sum(vapply(res, `[[`, logical(1), "singular")),
      sum(vapply(res, `[[`, logical(1), "psi_bound"))))
  }

  out_alpha <- alpha
  out_psi <- psi
  for (g in seq_len(ng)) {
    out_alpha[g, ] <- res[[g]]$alpha
    out_psi[[g]] <- res[[g]]$psi
  }
  polish <- data.frame(
    iterations = vapply(res, `[[`, integer(1), "iterations"),
    psi_fitnb = psi,
    restarted = vapply(res, `[[`, logical(1), "restarted"),
    capped = vapply(res, `[[`, logical(1), "capped"),
    singular = vapply(res, `[[`, logical(1), "singular"),
    # the dispersion optimum sat on its search bound, so fitNB's moderated psi
    # was kept for this gene rather than a boundary value stored as an estimate
    psi_bound = vapply(res, `[[`, logical(1), "psi_bound"),
    # FALSE when the gene fell back to fitNB's fit entirely
    polished = vapply(res, `[[`, logical(1), "polished"),
    row.names = rownames(alpha)
  )
  list(alpha = out_alpha, psi = out_psi,
       loglik = vapply(res, `[[`, numeric(1), "loglik"), polish = polish)
}

#' The per-column ridge the polish must respect
#'
#' A mixed fit carries its penalty vector; a fixed-effects fit carries none
#' (\code{@penalty} is what selects the mixed inference path, so it cannot be
#' set on a fixed fit) and the ridge is whatever \code{lambda.a} the fit was
#' made with.
#' @noRd
.polishPenalty <- function(penalty, lambda.a, p) {
  if (!is.null(penalty)) return(penalty)
  if (length(lambda.a) == 1L) rep(lambda.a, p) else lambda.a
}

#' The variance-component fixed-point iteration, Steffensen-accelerated
#'
#' Schall's update is a fixed-point map whose convergence is linear, at a
#' rate that is the fraction of "missing information" in the random-effect
#' block: fast for a well-identified per-sample intercept (within 1% of its
#' limit in two steps on the cohort), slow for the (sample x cell type) block
#' (ratio 0.5-0.7 per step on the cohort, where a cap of three left it 6-12%
#' above its extrapolated limit in thirteen of fifteen runs; on the clustered
#' fixture, whose nested component is truly zero, the map creeps toward the
#' floor sublinearly). Aitken's delta-squared extrapolation after every two
#' plain steps (Steffensen's method) reaches a geometric limit exactly and a
#' near-geometric one in a few steps; it is applied per component, and only
#' when the last two steps contract monotonically (ratio in (0, 0.95)), so an
#' oscillating or stalled sequence is left to the plain map. The extrapolated
#' value is clamped to \code{range}.
#'
#' @param tau2 named list, the starting components.
#' @param step a function of the current components returning the next Schall
#'   values (same shape), or \code{NULL} to stop and keep the current ones.
#' @param maxit,tol iteration cap and tolerance on the largest change in
#'   \code{log(tau2)} across components.
#' @param accelerate logical; \code{FALSE} is the plain map.
#' @param range the clamp for an extrapolated value.
#' @param verbose report each step.
#' @return a list with \code{tau2}, \code{iterations} (steps taken) and
#'   \code{converged}.
#' @noRd
.tau2Iterate <- function(tau2, step, maxit = 10L, tol = 1e-2, accelerate = TRUE,
                         range = c(1e-8, 1e4), verbose = FALSE) {
  x0 <- NULL          # the point two plain steps ago (the Steffensen window)
  x1 <- NULL
  converged <- FALSE
  it <- 0L
  for (k in seq_len(maxit)) {
    tau2_new <- step(tau2)
    if (is.null(tau2_new)) break
    if (!all(is.finite(unlist(tau2_new)))) {
      stop("a variance component is non-finite after a Schall step (",
           paste(sprintf("%s=%s", names(tau2_new),
                         format(unlist(tau2_new), digits = 3)), collapse = ", "),
           "). tau2 is shared by every gene of the bandwidth, so the stage ",
           "stops here rather than carry a bad component into every gene's ",
           "penalty and reference df; a non-finite coefficient in the polished ",
           "fit (see @polish) is the usual cause.", call. = FALSE)
    }
    it <- k
    note <- ""
    if (accelerate && !is.null(x0)) {
      # x0 -> x1 -> tau2_new are two plain steps from the window's start
      for (g in names(tau2_new)) {
        d1 <- x1[[g]] - x0[[g]]
        d2 <- tau2_new[[g]] - x1[[g]]
        r <- if (d1 != 0) d2 / d1 else NA_real_
        if (is.finite(r) && r > 0 && r < 0.95) {
          ext <- tau2_new[[g]] - d2^2 / (d2 - d1)
          tau2_new[[g]] <- min(max(ext, range[1]), range[2])
          note <- " (extrapolated)"
        }
      }
      x0 <- NULL
      x1 <- NULL
    } else if (accelerate) {
      if (is.null(x0)) {
        x0 <- tau2
        x1 <- tau2_new
      }
    }
    delta <- max(abs(log(unlist(tau2_new)) - log(unlist(tau2))))
    if (verbose) {
      message(sprintf("  tau2 from the converged fit: %s%s",
                      paste(sprintf("%s=%.3g", names(tau2_new), unlist(tau2_new)),
                            collapse = ", "), note))
    }
    tau2 <- tau2_new
    if (delta < tol) {
      converged <- TRUE
      break
    }
  }
  list(tau2 = tau2, iterations = it, converged = converged)
}

#' Polish one SpiDEFit in place: converged coefficients, inference invalidated
#' @noRd
.polishSpiDEFit <- function(f, Y, lambda.a = 0, maxit = 50L, tol = 1e-8,
                            block.size = NULL,
                            BPPARAM = BiocParallel::SerialParam(),
                            verbose = TRUE, psi.method = "profile",
                            tau2 = TRUE, tau2.maxit = 10L, tau2.tol = 1e-2,
                            tau2.accelerate = TRUE,
                            tau2.range = c(1e-8, 1e4)) {
  f <- updateObject(f)
  Yf <- Y[rownames(f@alpha), , drop = FALSE]
  pen <- .polishPenalty(f@penalty, lambda.a, ncol(f@W))
  run_polish <- function(alpha0, psi0, pen_now, warm = FALSE) {
    .polishFit(Yf, f@W, alpha0, psi0, pen_now, f@re_group,
               covtype = as.character(f@covtype),
               maxit = maxit, tol = tol, block.size = block.size,
               BPPARAM = BPPARAM, verbose = verbose, psi.method = psi.method,
               warm = warm)
  }
  pol <- run_polish(f@alpha, f@psi, pen)
  alpha <- pol$alpha
  dimnames(alpha) <- dimnames(f@alpha)
  psi <- as.numeric(pol$psi)
  # the diagnostics are the cold pass's (its iterations, restarts and fitNB's
  # psi); the re-polish passes below add their Newton iterations to one column
  polish <- pol$polish
  repolish_it <- integer(nrow(polish))
  repolish_capped <- logical(nrow(polish))
  repolish_singular <- logical(nrow(polish))

  # The variance components from the CONVERGED fit. The fit's Schall loop
  # reads the shared fit's own coefficients and dispersion, which can sit far
  # from every gene's optimum (research fdr-ordering/FINDINGS.md, 2026-09-08:
  # a between-sample variance of 10 against a planted 0.49 on the clustered
  # fixture, from sample intercepts three times too wide). One Schall step on
  # the polished coefficients with the gene-averaged weights at the polished
  # mean and dispersion, then a re-polish at the new penalty, iterated to
  # tolerance; the Satterthwaite df follows below.
  mixed <- !is.null(f@re_group) && !is.null(f@tau2) && length(f@tau2)
  if (tau2 && mixed) {
    # the components the current coefficients are converged at: the fit's,
    # after the cold pass. A step only re-polishes when they change.
    pen_tau2 <- f@tau2
    repolish_at <- function(tau2_now) {
      for (g in names(tau2_now)) pen[which(f@re_group == g)] <<- 1 / tau2_now[[g]]
      pol <- run_polish(alpha, psi, pen, warm = TRUE)
      alpha <<- pol$alpha
      dimnames(alpha) <<- dimnames(f@alpha)
      repolish_it <<- repolish_it + pol$polish$iterations
      # a warm pass that hit its cap or a singular system leaves that gene at
      # its previous converged fit; the flags must reach @polish, not only a
      # verbose message
      repolish_capped <<- repolish_capped | pol$polish$capped
      repolish_singular <<- repolish_singular | pol$polish$singular
      pen_tau2 <<- tau2_now
    }
    schall <- function(tau2_now) {
      if (any(unlist(tau2_now) != unlist(pen_tau2))) repolish_at(tau2_now)
      wbar <- .repWeights(Yf, alpha, f@W, psi, winsor = Inf)
      A <- crossprod(f@W * sqrt(wbar))
      minv <- tryCatch(SpaNorm::invert_mat(A + diag(pen)), error = function(e) NULL)
      if (is.null(minv)) {
        warning("the penalised information at the converged fit is singular; ",
                "the variance components are left where the loop reached",
                call. = FALSE)
        return(NULL)
      }
      .schallStep(alpha, minv, f@re_group, tau2_now, tau2.range)
    }
    loop <- .tau2Iterate(f@tau2, schall, maxit = tau2.maxit, tol = tau2.tol,
                         accelerate = tau2.accelerate, range = tau2.range,
                         verbose = verbose)
    tau2_now <- loop$tau2
    # the coefficients at the REPORTED penalty: the loop ends on a step (or an
    # extrapolation) it has not re-polished at, so one more warm pass -- a
    # Newton step or two per gene -- makes alpha, penalty and df consistent
    if (any(unlist(tau2_now) != unlist(pen_tau2))) repolish_at(tau2_now)
    f@tau2 <- tau2_now
    f@penalty <- pen
    # the reference df reads the components and the penalty
    if (!is.null(names(f@df))) {
      wbar <- .repWeights(Yf, alpha, f@W, psi, winsor = Inf)
      A <- crossprod(f@W * sqrt(wbar))
      minv <- tryCatch(SpaNorm::invert_mat(A + diag(pen)), error = function(e) NULL)
      tested <- match(names(f@df), colnames(f@W))
      df_new <- if (is.null(minv)) NULL else
        .satterthwaiteDF(A, minv, pen, f@re_group, tau2_now, tested, ncol(Yf), names(f@df))
      if (is.null(df_new)) {
        warning("the Satterthwaite reference df could not be refreshed at the ",
                "reported penalty (singular penalised information); it is left ",
                "as the fit computed it, at the fit's variance components",
                call. = FALSE)
      } else {
        f@df <- df_new
      }
    }
  }
  polish$repolish.iterations <- repolish_it
  polish$repolish.capped <- repolish_capped
  polish$repolish.singular <- repolish_singular
  rownames(polish) <- rownames(f@alpha)
  f@alpha <- alpha
  f@psi <- psi
  f@polish <- polish
  f@loglik <- as.numeric(.blockLoglik(Yf, alpha, f@W, f@psi, winsor = Inf))
  # everything inference derived from the old coefficients is now stale
  f@t_stat <- NULL
  f@se <- NULL
  f@p.combined.pos <- NULL
  f@p.combined.neg <- NULL
  f@se_patient <- numeric(0)
  f@rho <- numeric(0)
  f@two.sided <- FALSE
  f
}

#' The polish stage: converge each gene, set its dispersion, re-estimate the
#' variance components
#'
#' The stage between [fitSpiDE()] and [testSpiDE()] (the pipeline is fit ->
#' polish -> test -> gsea; [spiDE()] runs it by default). It converges every
#' gene to its own penalised negative-binomial optimum by damped Newton,
#' re-estimates its dispersion at the converged mean (\code{psi}), and for a
#' mixed fit re-estimates the variance components from the converged
#' coefficients and refreshes the Satterthwaite reference df (\code{tau2}).
#' Run it, or skip it, according to the data: it needs integer counts.
#'
#' Why it exists: \code{fitNB} fits all genes in one IRLS loop with a shared
#' cell-weight vector and an aggregate convergence criterion, and for bright,
#' cell-type-restricted genes it stops one to four standard errors short of
#' the gene's own optimum, with a dispersion estimated off the optimum that is
#' too large. The fit's Schall loop reads those same unconverged coefficients
#' and dispersion, so its variance components inherit the error: on the
#' clustered test fixture it reports a between-sample variance of 10 against
#' a planted 0.49, which this stage brings back to the planted value.
#'
#' Every inference slot derived from the old coefficients (\code{t_stat},
#' \code{se}, the combined p-values, the results table and the cross-bandwidth
#' weights) is cleared; call [testSpiDE()] again afterwards.
#'
#' @param object a SpiDEResults from [fitSpiDE()] (one GLM fit per bandwidth).
#'   An object without per-gene GLM fits is refused.
#' @param spe the SpatialExperiment the fit was made from (for the counts).
#' @param assay a character, the counts assay.
#' @param lambda.a the ridge the fit was made with, for a fixed-effects fit
#'   (\code{random = "none"}), which does not record it; ignored for a mixed
#'   fit, which carries its own penalty vector.
#' @param maxit,tol iteration cap and relative log-likelihood tolerance per
#'   gene (the \code{converge.maxit} / \code{converge.tol} of [fitSpiDE()]).
#' @param block.size genes per block; \code{NULL} splits one block per
#'   \code{BPPARAM} worker.
#' @param BPPARAM a BiocParallelParam; the stage is blocked over genes. With
#'   more than one worker, each worker runs its BLAS single-threaded when
#'   RhpcBLASctl is installed: forked workers inherit the parent's OpenBLAS
#'   thread count, and oversubscribing the cores that way was measured at 9x
#'   per Newton step. The parent process is left as it was.
#' @param verbose report progress.
#' @param ... further arguments passed to the method.
#' @param psi how the dispersion is set at the converged mean:
#'   \code{"profile"} (the default) re-estimates each gene's dispersion by
#'   profile maximum likelihood at the converged mean; \code{"moderated"}
#'   keeps \code{fitNB}'s cross-gene moderated value and converges only the
#'   coefficients under it. The moderated value is whatever the shared fit
#'   left, which can be far from the gene's own (fifteen times, on the
#'   clustered fixture), and the variance-component step below needs a
#'   dispersion consistent with the converged mean.
#' @param tau2 logical; for a mixed fit, re-estimate the variance components
#'   from the converged fit (a Schall step on the polished coefficients, then
#'   a re-polish at the new penalty, iterated to a fixed point), and refresh
#'   the Satterthwaite reference df. The fit's own loop reads the shared fit's
#'   unconverged coefficients, which over-estimates the between-sample
#'   variance badly where that fit is off its optimum. Only the first pass is
#'   a cold polish; each re-polish is a few damped Newton steps per gene at
#'   the held dispersion from the previous converged fit (they agree with a
#'   cold re-polish to 1e-3 on the coefficients), and the stage ends with one
#'   more such pass at the reported penalty so coefficients, penalty and df
#'   are consistent.
#' @param tau2.maxit,tau2.tol iteration cap and tolerance on the largest
#'   change in \code{log(tau2)} across components for that re-estimate. The
#'   tolerance is where the estimates stop mattering downstream: on the
#'   clustered fixture a 5\% change in a component moves individual
#'   t-statistics by at most 0.05 (median below 0.001), a tenfold error in a
#'   near-zero nested component by up to 0.13.
#' @param tau2.accelerate logical; Steffensen-accelerate the fixed-point
#'   iteration (Aitken's extrapolation after every two plain steps, per
#'   component, only while the steps contract monotonically). Schall's map
#'   converges linearly and slowly for the (sample x cell type) block: capped
#'   at three plain steps it stopped 6-12\% above its limit on the cohort.
#'   \code{FALSE} is the plain map.
#' @return the object with converged \code{alpha} and \code{psi}, per-gene
#'   diagnostics in \code{@polish}, and inference cleared.
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = 20)
#' fit0 <- fitSpiDE(spe, condition = "condition", sigma = 20, random = "none",
#'                  verbose = FALSE)
#' fit1 <- polishSpiDE(fit0, spe, verbose = FALSE)
#' head(fits(fit1)[[1]]@polish)
#' @seealso [fitSpiDE()] for the stage before, [testSpiDE()] for the stage
#'   after.
#' @rdname polishSpiDE
#' @export
setMethod(
  "polishSpiDE",
  signature = "SpiDEResults",
  definition = function(object, spe, assay = "counts",
                        psi = c("profile", "moderated"), tau2 = TRUE,
                        tau2.maxit = 10L, tau2.tol = 1e-2,
                        tau2.accelerate = TRUE, lambda.a = 0,
                        maxit = 50L, tol = 1e-8, block.size = NULL,
                        BPPARAM = BiocParallel::SerialParam(), verbose = TRUE) {
    object <- updateObject(object)
    psi <- match.arg(psi)
    if (!length(object@fits)) {
      stop("nothing to polish: the object carries no per-gene GLM fit", call. = FALSE)
    }
    checkSPE(spe, assay = assay)
    Y <- SummarizedExperiment::assay(spe, assay)
    checkCounts(Y, integer.only = TRUE)
    missing_genes <- setdiff(rownames(object@fits[[1]]@alpha), rownames(Y))
    if (length(missing_genes)) {
      stop("the fit's genes are not all in the counts: e.g. ",
           paste(utils::head(missing_genes, 3), collapse = ", "), call. = FALSE)
    }
    object@fits <- lapply(seq_along(object@fits), function(i) {
      if (verbose) message(sprintf("Polishing bandwidth sigma = %s",
                                   object@sigma[i]))
      .polishSpiDEFit(object@fits[[i]], Y, lambda.a = lambda.a,
                      psi.method = psi, tau2 = tau2, tau2.maxit = tau2.maxit,
                      tau2.tol = tau2.tol, tau2.accelerate = tau2.accelerate,
                      maxit = maxit, tol = tol, block.size = block.size,
                      BPPARAM = BPPARAM, verbose = verbose)
    })
    names(object@fits) <- names(updateObject(object)@fits)
    # cross-bandwidth combination and the results table are stale too
    object@gene.weights <- NULL
    object@p.cauchy.pos <- NULL
    object@p.cauchy.neg <- NULL
    object@results <- data.frame()
    object@fdr <- NA_real_
    object
  }
)
