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
#
# The engine itself -- the per-gene and batched damped Newton, the Schur
# absorption of the random block, the profile dispersion -- lives in SpaNorm
# (SpaNorm::polishNB(), SpaNorm::nbProfilePsi(), SpaNorm::nbNewtonSolver()).
# This file keeps what is spiDE's: which columns a SpiDEFit absorbs and starts
# from, the penalty it must respect, and the variance-component loop around the
# polish.

#' Which columns the Newton solver should absorb, and how they block
#'
#' Returns what \code{SpaNorm::nbNewtonSolver()}'s and
#' \code{SpaNorm::polishNB()}'s \code{absorb} argument accepts:
#' \code{NULL} for a fixed-effects fit, a LOGICAL selecting the nested
#' indicators when there are no random slopes, and a per-column SAMPLE grouping
#' when there are.
#'
#' Why not always the sample grouping. Without slopes the nested indicators are
#' 0/1 and partition the cells, so their \code{C} is diagonal and the scalar
#' path inverts it by reciprocal -- that is the production path for
#' \code{random = "intercept"} and it is cheaper than a block Cholesky. Adding
#' the per-sample intercepts to the absorbed set would buy little (there are S
#' of them) and would move a path that is already measured.
#'
#' With slopes the nested indicators alone are not the absorbable set: a
#' sample's slope columns are not orthogonal to its intercepts, so the whole
#' random block has to go in, grouped by sample.
#'
#' This is the one place that decides, so a new random-effect group is handled
#' by \code{re_sample} carrying it rather than by another literal match on
#' \code{re_group} at each call site.
#' @noRd
.absorbSpec <- function(fit) {
  rg <- fit@re_group
  if (is.null(rg)) return(NULL)
  has_slope <- any(!is.na(rg) & rg == "SampleSlope")
  rs <- if (methods::.hasSlot(fit, "re_sample")) fit@re_sample else NULL
  if (!has_slope || is.null(rs) || length(rs) != length(rg)) {
    # no slopes, or a fit saved before re_sample existed: the nested indicators
    return(!is.na(rg) & rg == "SampleCellTypeInt")
  }
  rs
}

#' The absorption the shared-factor batched solver gets
#'
#' \code{SpaNorm::polishNB()}'s \code{absorb.batch}. The shared-factor batched
#' solver (the device path) absorbs 1x1 blocks only, so it gets the nested
#' (sample x cell type) indicators even when \code{.absorbSpec()} hands the
#' per-gene solver a slope fit's whole random block grouped by sample. Both are
#' exact, so this costs the device path a wider dense block and moves no
#' result. Without slopes it is the same logical \code{.absorbSpec()} returns;
#' for a fixed-effects fit it is \code{NULL}.
#' @noRd
.absorbBatchSpec <- function(fit) {
  if (is.null(fit@re_group)) NULL else
    !is.na(fit@re_group) & fit@re_group == "SampleCellTypeInt"
}

#' The columns the polish's sane start fills with cell-type log means
#'
#' \code{SpaNorm::polishNB()}'s \code{start.cols}: a gene whose starting point
#' is degenerate restarts from the log mean over each cell type's cells on its
#' cell-type intercept, every other coefficient zero. The one place the start
#' columns are read off the design's tags.
#' @noRd
.testedStartCols <- function(fit) {
  as.character(fit@covtype) == "CellType"
}

#' Refuse a ridge penalty of the wrong length, naming re.celltype
#'
#' \code{SpaNorm::polishNB()} checks the length too, but its message cannot
#' name a spiDE argument. The usual way to get here is a \code{lambda.a}
#' vector sized for a design before \code{re.celltype = TRUE} added its
#' columns, so spiDE checks first and says so.
#' @noRd
.checkPolishPenalty <- function(pen, p) {
  if (!length(pen) %in% c(1L, p)) {
    stop("'lambda.a' must be a single value or one per design column (",
         p, " here, ", length(pen), " supplied). Note re.celltype = TRUE ",
         "adds one column per non-empty (sample, cell type), so a vector sized ",
         "for an earlier design is now too short.", call. = FALSE)
  }
  invisible(TRUE)
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
#' @param maxit,tol iteration cap and tolerance on the estimated distance of
#'   \code{log(tau2)} from its fixed point, the largest over components.
#' @param accelerate logical; \code{FALSE} is the plain map.
#' @param range the clamp for an extrapolated value.
#' @param verbose report each step.
#' @return a list with \code{tau2}, \code{iterations} (steps taken),
#'   \code{converged} and \code{trace} (one row per step: the components,
#'   whether the step was extrapolated, the log-step and the estimated gap).
#' @noRd
.tau2Iterate <- function(tau2, step, maxit = 10L, tol = 1e-2, accelerate = TRUE,
                         range = c(1e-8, 1e4), verbose = FALSE) {
  x0 <- NULL          # the point two plain steps ago (the Steffensen window)
  x1 <- NULL
  d_prev <- NULL      # the previous plain log-step per component, for the gap estimate
  converged <- FALSE
  it <- 0L
  trace <- list()
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
    extrapolated <- FALSE
    plain <- tau2_new
    if (accelerate && !is.null(x0)) {
      # x0 -> x1 -> tau2_new are two plain steps from the window's start
      for (g in names(tau2_new)) {
        d1 <- x1[[g]] - x0[[g]]
        d2 <- tau2_new[[g]] - x1[[g]]
        r <- if (d1 != 0) d2 / d1 else NA_real_
        if (is.finite(r) && r > 0 && r < 0.95) {
          ext <- tau2_new[[g]] - d2^2 / (d2 - d1)
          tau2_new[[g]] <- min(max(ext, range[1]), range[2])
          extrapolated <- TRUE
        }
      }
      x0 <- NULL
      x1 <- NULL
    } else if (accelerate && is.null(x0)) {
      x0 <- tau2
      x1 <- tau2_new
    }
    # The stopping rule is the estimated DISTANCE to the fixed point, not the
    # size of the last step. For a linearly contracting map with ratio r the
    # gap left after a step d is about d r / (1 - r), which a step-size rule
    # ignores: at r = 0.9 it would declare convergence with ten times the
    # tolerance still to go. The ratio comes from the last two plain steps;
    # without one (the first step, or the step after an extrapolation) the
    # step itself stands in for the gap, and a ratio at or above the Aitken
    # ceiling means the gap cannot be bounded, so the loop goes on.
    d_cur <- log(unlist(tau2_new)) - log(unlist(tau2))
    gap <- abs(d_cur)
    if (!extrapolated && !is.null(d_prev)) {
      r <- d_cur / d_prev
      ok <- is.finite(r) & r > 0 & r < 0.95
      gap[ok] <- abs(d_cur[ok]) * r[ok] / (1 - r[ok])
      gap[is.finite(r) & r >= 0.95] <- Inf
    }
    d_prev <- if (extrapolated) NULL else d_cur
    delta <- max(abs(d_cur))
    trace[[k]] <- data.frame(iteration = k, as.list(unlist(tau2_new)),
                             extrapolated = extrapolated, delta = delta,
                             gap = max(gap), check.names = FALSE)
    if (verbose) {
      note <- if (extrapolated) sprintf(" (extrapolated; the plain step gave %s)",
                                        paste(sprintf("%s=%.3g", names(plain), unlist(plain)), collapse = ", ")) else ""
      message(sprintf("  tau2 from the converged fit: %s%s",
                      paste(sprintf("%s=%.3g", names(tau2_new), unlist(tau2_new)),
                            collapse = ", "), note))
    }
    tau2 <- tau2_new
    if (max(gap) < tol) {
      converged <- TRUE
      break
    }
  }
  list(tau2 = tau2, iterations = it, converged = converged,
       trace = if (length(trace)) do.call(rbind, trace) else NULL)
}

#' Polish one SpiDEFit in place: converged coefficients, inference invalidated
#' @noRd
.polishSpiDEFit <- function(f, Y, lambda.a = 0, maxit = 50L, tol = 1e-8,
                            block.size = NULL,
                            BPPARAM = BiocParallel::SerialParam(),
                            verbose = TRUE, psi.method = "profile",
                            tau2 = TRUE, tau2.maxit = 10L, tau2.tol = 1e-2,
                            tau2.accelerate = TRUE,
                            tau2.range = c(1e-8, 1e4),
                            engine = c("batch", "gene"), batch.size = NULL,
                            backend = c("cpu", "auto", "gpu"),
                            gpu.mem.budget = NULL) {
  engine <- match.arg(engine)
  backend <- match.arg(backend)
  psi.method <- match.arg(psi.method, c("profile", "moderated"))
  f <- updateObject(f)
  Yf <- Y[rownames(f@alpha), , drop = FALSE]
  pen <- .polishPenalty(f@penalty, lambda.a, ncol(f@W))
  # "moderated" keeps fitNB's dispersion, which SpaNorm calls holding it fixed
  psi_rule <- if (psi.method == "moderated") "fixed" else "profile"
  # gpu.mem.budget is not forwarded: it sizes the inference stage's device
  # batches, and SpaNorm::polishNB() sizes its own
  run_polish <- function(alpha0, psi0, pen_now, warm = FALSE) {
    .checkPolishPenalty(pen_now, ncol(f@W))
    SpaNorm::polishNB(
      Yf, f@W, alpha0, psi0, lambda.a = pen_now,
      absorb = .absorbSpec(f), absorb.batch = .absorbBatchSpec(f),
      start.cols = .testedStartCols(f), psi.method = psi_rule,
      warm = warm, maxit = maxit, tol = tol, engine = engine,
      batch.size = batch.size, block.size = block.size, backend = backend,
      BPPARAM = BPPARAM, verbose = verbose)
  }
  pol <- run_polish(f@alpha, f@psi, pen)
  alpha <- pol$alpha
  psi <- as.numeric(pol$psi)
  # the diagnostics are the cold pass's (its iterations, restarts and fitNB's
  # psi); the re-polish passes below add their Newton iterations to one column
  polish <- pol$polish
  repolish_it <- integer(nrow(polish))
  repolish_capped <- logical(nrow(polish))
  repolish_singular <- logical(nrow(polish))
  # a gene the cold pass could not polish keeps fitNB's fit through the loop:
  # a warm Newton from fitNB's degenerate point has none of the cold path's
  # guards, and @polish says "polished = FALSE" for it
  cold_ok <- pol$polish$polished
  loop <- NULL

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
      pen <<- .penaltyFromTau2(pen, f@re_group, tau2_now)
      pol <- run_polish(alpha, psi, pen, warm = TRUE)
      new_alpha <- pol$alpha
      new_alpha[!cold_ok, ] <- alpha[!cold_ok, ]
      alpha <<- new_alpha
      repolish_it <<- repolish_it + ifelse(cold_ok, pol$polish$iterations, 0L)
      # a warm pass that hit its cap, met a singular system or fell back on a
      # non-finite result leaves that gene at its previous converged fit; the
      # flags must reach @polish, not only a verbose message
      repolish_capped <<- repolish_capped | pol$polish$capped
      repolish_singular <<- repolish_singular | pol$polish$singular | !pol$polish$polished
      pen_tau2 <<- tau2_now
    }
    schall <- function(tau2_now) {
      if (any(unlist(tau2_now) != unlist(pen_tau2))) repolish_at(tau2_now)
      info <- .penalisedInfo(Yf, alpha, f@W, psi, pen, winsor = Inf)
      if (is.null(info$minv)) {
        warning("the penalised information at the converged fit is singular; ",
                "the variance components are left where the loop reached",
                call. = FALSE)
        return(NULL)
      }
      .schallStep(alpha, info$minv, f@re_group, tau2_now, tau2.range)
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
      info <- .penalisedInfo(Yf, alpha, f@W, psi, pen, winsor = Inf)
      tested <- match(names(f@df), colnames(f@W))
      df_new <- if (is.null(info$minv)) NULL else
        .satterthwaiteDF(info$A, info$minv, pen, f@re_group, tau2_now, tested, ncol(Yf), names(f@df))
      if (is.null(df_new)) {
        # fail toward validity, as the fit does: the between-sample scalar
        f@df <- .betweenDF(f@re_group, .fitMode(f), ncol(Yf))
        warning("the Satterthwaite reference df could not be refreshed at the ",
                "reported penalty (singular penalised information); degraded to ",
                "the conservative between-sample df ", format(f@df),
                call. = FALSE)
      } else {
        # the same bound the fit applies: a between-patient contrast cannot
        # out-run its patients, and this site would otherwise overwrite it
        f@df <- .boundPatientDF(df_new, f@W, f@re_group, tested,
                                .betweenDF(f@re_group, .fitMode(f), ncol(Yf)),
                                .patientsPerTested(f@W, f@re_group, f@coefmap,
                                                   tested),
                                grepl("Response", as.character(f@covtype)[tested]))
      }
    }
    # the warm passes held each gene's dispersion; one profile pass at the
    # final coefficients puts @psi at the reported mean (on the clustered
    # fixture the held value sat 1-5% below its optimum after a loop that
    # moved the nested component 35-fold, worth 0.014 in t at most)
    if (psi.method == "profile" && loop$iterations > 0L) {
      psi <- SpaNorm::nbProfilePsi(Yf, f@W, alpha, psi, block.size = block.size,
                                   BPPARAM = BPPARAM)
    }
  }
  polish$repolish.iterations <- repolish_it
  polish$repolish.capped <- repolish_capped
  polish$repolish.singular <- repolish_singular
  if (!is.null(loop)) {
    attr(polish, "tau2") <- list(iterations = loop$iterations,
                                 converged = loop$converged, trace = loop$trace)
  }
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
#' @param engine \code{"batch"} (the default) converges a block of genes
#'   together, so the design is read once per batch rather than once per gene;
#'   \code{"gene"} is the original per-gene loop, kept as the reference
#'   implementation. The two agree to ~5e-13 on the coefficients with identical
#'   convergence flags and iteration counts, but \code{"batch"} is not
#'   invariant to the batch boundary at machine precision -- BLAS blocks a
#'   many-row product differently from a one-row one.
#' @param backend \code{"cpu"} (the default), or \code{"gpu"}/\code{"auto"}
#'   to run the batched Newton on an accelerator when one is present. The
#'   device path requires \code{engine = "batch"} -- the per-gene engine keeps
#'   one factorisation object per gene, which is what cannot go to a device --
#'   and refuses a single-precision device outright. It is an accelerator,
#'   never a requirement: without one, \code{"gpu"} gives the CPU answer.
#' @param gpu.mem.budget bytes for the device, or NULL to auto-detect.
#' @param batch.size genes per batched Newton. The default comes from a memory
#'   budget (\code{options(spiDE.polish.mem.budget = )}, bytes per worker),
#'   because the batched working set is gene x cell: at 77,454 cells a
#'   2,000-gene block would allocate over a terabyte, so the gene block size
#'   cannot be the batch size.
#' @param BPPARAM a BiocParallelParam; the stage is blocked over genes. With
#'   more than one worker, each worker runs its BLAS single-threaded when
#'   RhpcBLASctl is installed (forked workers inherit the parent's thread
#'   count, and oversubscribing the cores that way costs an order of
#'   magnitude per gene); the parent process is left as it was.
#' @param verbose report progress.
#' @param ... further arguments passed to the method.
#' @param psi how the dispersion is set at the converged mean:
#'   \code{"profile"} (the default) re-estimates each gene's dispersion by
#'   profile maximum likelihood at the converged mean; \code{"moderated"}
#'   keeps \code{fitNB}'s cross-gene moderated value and converges only the
#'   coefficients under it. The moderated value is whatever the shared fit
#'   left, which can be far from the gene's own on a small fixture. The two
#'   rules were measured from one fit polished both ways, on the synthetic
#'   benchmark and on a real cohort, and are equivalent in calibration and
#'   power; \code{"moderated"} is the cheaper one, since it skips the
#'   dispersion search and the two re-polishes it triggers. The numbers are in
#'   \code{vignette("spiDE-calibration")}. Under \code{"profile"} the stage
#'   ends with one profile pass at the reported coefficients, so \code{@psi}
#'   is the optimum at the mean the stage reports.
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
#'   diagnostics in \code{@polish} (the cold pass's \code{iterations},
#'   \code{restarted}, \code{capped}, \code{singular}, \code{psi_bound},
#'   \code{polished} and \code{psi_fitnb}, plus \code{repolish.iterations},
#'   \code{repolish.capped} and \code{repolish.singular} accumulated over the
#'   warm passes, and for a mixed fit the attribute \code{"tau2"}: the loop's
#'   \code{iterations}, \code{converged} and per-step \code{trace}), and
#'   inference cleared.
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
                        BPPARAM = BiocParallel::SerialParam(), verbose = TRUE,
                        engine = c("batch", "gene"), batch.size = NULL,
                        backend = c("cpu", "auto", "gpu"),
                        gpu.mem.budget = NULL) {
    object <- updateObject(object)
    psi <- match.arg(psi)
    engine <- match.arg(engine)
    backend <- match.arg(backend)
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
                      BPPARAM = BPPARAM, verbose = verbose,
                      engine = engine, batch.size = batch.size,
                      backend = backend, gpu.mem.budget = gpu.mem.budget)
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
