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
    A <- crossprod(X * sqrt(w))
    diag(A) <- diag(A) + pen_x
    agg <- rowsum(cbind(w, X * w), group = gf, reorder = TRUE)
    cvec <- agg[, 1] + pen_z              # the diagonal of C, length G
    B <- t(agg[, -1, drop = FALSE])       # ncol(X) x G
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
#' @return a list with \code{alpha}, \code{psi}, \code{loglik},
#'   \code{iterations}, \code{restarted}, \code{capped}, \code{singular},
#'   \code{psi_bound} and \code{polished}.
#' @importFrom stats optimize
#' @noRd
.polishGene <- function(y, W, a0, psi0, pen, solver, maxit = 50L, tol = 1e-8,
                        ct_cols = NULL, psi.range = c(1e-3, 1e3)) {
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
    mu <- as.numeric(exp(W %*% a))
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
        mu1 <- as.numeric(exp(W %*% a1))
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
  for (k in 1:2) {
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
#' @param verbose logical.
#' @return a list with \code{alpha}, \code{psi}, \code{loglik} and a per-gene
#'   \code{polish} data.frame.
#' @importFrom BiocParallel bplapply SerialParam bpnworkers
#' @noRd
.polishFit <- function(Y, W, alpha, psi, pen, re_group = NULL, covtype = NULL,
                       maxit = 50L, tol = 1e-8, block.size = NULL,
                       BPPARAM = BiocParallel::SerialParam(), verbose = FALSE) {
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
    stop("converge = TRUE needs integer counts: the negative binomial ",
         "likelihood is undefined otherwise, and every gene's dispersion would ",
         "silently collapse to its upper bound.\n  The assay passed is not ",
         "integer-valued (e.g. a back-transform such as 2^logcounts - 1).\n  ",
         "Use the raw counts, or pass converge = FALSE.", call. = FALSE)
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
  if (verbose) {
    message(sprintf("  converging %d genes per gene (%d block%s)", ng,
                    length(blocks), if (length(blocks) == 1L) "" else "s"))
  }
  # A whole-transcriptome polish is hours of work, so report progress rather
  # than going silent after the opening message. Blocks are timed as they
  # finish; under a parallel BPPARAM they complete out of order, so the count
  # is of blocks retired, not a position in the gene list.
  t0 <- Sys.time()
  nb <- length(blocks)
  step <- max(1L, nb %/% 20L)
  res <- BiocParallel::bplapply(seq_along(blocks), function(b) {
    gi <- blocks[[b]]
    # densify the whole block once: Y may be sparse or a DelayedArray, where a
    # per-gene read costs a round trip each time (the invariant is that the
    # WHOLE matrix is never densified, not that a block is never densified --
    # .blockedInference() does exactly the same).
    Yb <- as.matrix(Y[gi, , drop = FALSE])
    out <- lapply(seq_along(gi), function(i) {
      g <- gi[[i]]
      .polishGene(as.numeric(Yb[i, ]), W, alpha[g, ], psi[[g]], pen, solver,
                  maxit = maxit, tol = tol, ct_cols = ct_cols)
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
