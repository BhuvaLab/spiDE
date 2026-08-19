# Two-stage estimation: per-core niche slopes, pooled by precision within each
# patient, then a patient-level contrast via weighted moderated limma.
#
# WHY THIS IS NOT A `random` MODE
# -------------------------------
# `random = "none"/"intercept"/"slope"` all fit the SAME design matrix and differ
# only in which columns are ridge-penalised. This is a different ESTIMATOR: the
# niche slope is estimated per core (sample), the per-core slopes are pooled by
# inverse-variance precision within each patient, and the condition contrast is
# then taken over the pooled patient slopes with limma. None of `W`, `alpha`,
# `psi`, `penalty`, `tau2` or the Satterthwaite df machinery applies, so
# slotting it under `random` would hand every downstream consumer an object
# whose fields do not mean what they claim.
#
# WHAT IT FIXES
# -------------
# `condition` is assigned per PATIENT, so patients are the experimental units.
# The one-stage model treats cells as replicates of a patient-level contrast; on
# a patient-label permutation of the YTMA cohort (55 patients, every triplet null
# by construction) that returns ~576 false calls in every replicate, realized FDR
# 1.00 against a nominal 0.05. NOTE the two-stage numbers below are HISTORICAL,
# measured on this function's predecessor (the nbresid/Welch estimator) -- see
# the roxygen: that predecessor gave raw type-I 0.036 and ZERO false calls on
# the same permuted data, in 1.5 minutes against 128. Re-measurement on the
# current estimator is pending.
#
# WHAT IT DOES NOT FIX
# --------------------
# Multiplicity. The estimator carries real signal (true triplets enriched 5x at
# alpha 0.05 and 65x at alpha 0.001 on a spiked plasmode -- also predecessor
# measurements, pending re-measurement) but a full-panel
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
#' anything visible in the output. Samples a compact genes x cells BLOCK
#' (never a full-length linear index, which forces a full realisation for a
#' sparse/DelayedArray assay) and densifies only that block, so this check
#' never materialises the whole matrix.
#' @noRd
.looksLikeCounts <- function(Y, n = 2000L) {
  ng <- nrow(Y); nc <- ncol(Y)
  if (!ng || !nc) return(FALSE)
  if (as.double(ng) * as.double(nc) <= n) {
    v <- as.numeric(as.matrix(Y))
  } else {
    nr <- max(1L, min(ng, as.integer(ceiling(sqrt(n)))))
    ncl <- max(1L, min(nc, as.integer(ceiling(n / nr))))
    # evenly spaced, not sampled: sample.int() would advance the global RNG
    # stream (a seeded pipeline's downstream draws would then depend on
    # whether this gate ran at all), and a borderline assay could pass the
    # check on one run and fail it on the next
    ri <- unique(as.integer(round(seq(1L, ng, length.out = nr))))
    ci <- unique(as.integer(round(seq(1L, nc, length.out = ncl))))
    v <- as.numeric(as.matrix(Y[ri, ci, drop = FALSE]))
  }
  v <- v[is.finite(v)]
  if (!length(v)) return(FALSE)
  all(v >= 0) && all(abs(v - round(v)) < 1e-8)
}

#' First non-missing value of a cell-level vector, per patient
#'
#' Indexes the original vector rather than going through tapply(), which
#' unlists a factor into bare integer level codes -- a factor patient
#' covariate would then enter the stage-2 design as a continuous trend in
#' arbitrary codes, silently. Taking the first NON-missing cell also keeps a
#' patient whose first cell happens to be NA while its value is known
#' elsewhere.
#' @noRd
.patientValue <- function(x, pat, pats) {
  idx <- vapply(pats, function(p) {
    i <- which(pat == p & !is.na(x))
    if (length(i)) i[1L] else which(pat == p)[1L]
  }, integer(1))
  stats::setNames(x[idx], pats)
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
#' 0.05). \strong{These next two numbers are historical}, measured on this
#' function's predecessor (the nbresid/Welch two-stage estimator, since
#' replaced by the SpaNorm-anchored joint estimator documented here), not on
#' the current implementation; re-measurement on the current estimator is
#' pending. That predecessor returned raw type-I 0.036 and zero false calls,
#' ~100x faster.
#'
#' \strong{It does not solve multiplicity.} Also measured on that same
#' predecessor estimator, pending re-measurement here: on a spiked plasmode the
#' true triplets were enriched 5x at \eqn{\alpha = 0.05} and 65x at
#' \eqn{\alpha = 0.001}, yet nothing survived BH across a full-panel space of
#' ~1.8 million triplets, and the hierarchical cascade did not help (combining a
#' gene's ~137 mostly-null triplets dilutes a single true effect). Restrict
#' \code{index}, \code{niche} and the gene set to a pre-specified hypothesis:
#' roughly 4,000 tests is the order at which a \eqn{p \approx 10^{-5}} effect
#' survives.
#'
#' @seealso The model vignette (`vignette("spiDE-model")`) documents the
#'   estimator's two stages with full equations; the *Two-stage estimation*
#'   benchmark report on the spiDE-research site
#'   (<https://bhuvalab.github.io/spiDE-research/>) reports its measured operating
#'   characteristics against the published simulation study, paired on the
#'   same simulated datasets.
#' @param spe a SpatialExperiment with niche reducedDims (see [buildNiches()]).
#' @param condition a character, the colData column of the tested condition,
#'   constant within patient.
#' @param sigma a numeric, the bandwidth (one value; slopes are per bandwidth).
#' @param index,niche character vectors restricting index / niche cell types.
#'   \strong{Use them}: the full space is usually too large to detect anything.
#'   Stage-1 slopes are estimated \emph{jointly} across niche columns (see
#'   \code{stage1}), so restricting \code{niche} changes the estimated slope
#'   for every remaining niche, not merely which rows of the results table are
#'   reported.
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
#' @param stage1 one of "spanorm" (the default), "ols" or "nb", the model
#'   stage-1 slopes are estimated from.
#'
#'   \strong{"spanorm"} fits an NB GLM on the raw counts per (sample, index)
#'   subset, with design \code{[1, niche columns]} and the stored
#'   [SpaNorm::SpaNorm()] fit's library-size and batch linear predictor as a
#'   \emph{fixed offset} (read from \code{metadata(spe)$SpaNorm}; a clear
#'   error names \code{SpaNorm::SpaNorm()} when none is found). The rationale:
#'   SpaNorm models biology only to \emph{anchor} its normalisation -- the
#'   fitted biology term is a smooth catch-all, not modelled biology -- so
#'   spiDE models all biology from scratch and assumes only that the
#'   library-size (and within-sample batch, e.g. field-of-view) effects are
#'   right. Fixing their predictor at coefficient 1 is that assumption made
#'   literal: unlike a fitted library-size covariate, an offset cannot absorb
#'   depth-correlated biology. (Before spiDE 0.99.16 this option instead used
#'   a one-step "addback" working response linearised at the fitted mean;
#'   results from that construction are not comparable and should be
#'   re-computed.)
#'
#'   \strong{"ols"} regresses log-CPM directly on the niche columns with unit
#'   weights. It needs no stored SpaNorm fit and no dispersion estimate, so it
#'   is the fallback when one is not available (e.g. \code{data(toySpiDE)}).
#'
#'   \strong{"nb"} fits a fresh \code{SpaNorm::fitNB} per (sample, index)
#'   subset carrying \code{[1, log-library-size, niche columns]}. Note the
#'   library-size term is a \emph{fitted covariate} here, so it can absorb a
#'   depth-correlated effect; prefer "spanorm" when a SpaNorm fit is
#'   available. Both NB paths price per-GENE (dispersion estimation and IRLS
#'   setup), not per cell.
#' @param epsilon \strong{Deprecated and ignored.} The former
#'   "addback"/"residual" choice applied to the pre-0.99.16 working-response
#'   construction, which no longer exists; supplying the argument raises a
#'   warning.
#' @param pool.psi logical (default TRUE): for \code{stage1 = "spanorm"},
#'   estimate one per-gene dispersion per SAMPLE (design: cell-type means plus
#'   the ls/batch offset) and supply it to every (sample, index) subset fit,
#'   instead of re-estimating per subset. Cuts the dominant stage-1 cost
#'   (dispersion estimation, measured at ~half of each subset fit on the real
#'   cohort) and estimates psi from all of a sample's cells rather than one
#'   type's few. The pooling design omits the niche columns, which errs
#'   slightly conservative. Requires SpaNorm with \code{fitNB(psi=)}; on an
#'   older SpaNorm the option is silently inert and per-subset estimation is
#'   used.
#' @param winsor,lambda.a,maxit.psi,backend forwarded to
#'   \code{\link[SpaNorm]{fitNB}} by the "nb" path.
#' @param verbose a logical.
#' @return a [SpiDEResults] with the tidy \code{results} table populated
#'   (unchanged schema: \code{gene}, \code{ct_index}, \code{ct_niche},
#'   \code{coef}, \code{t}, \code{p.niche}, \code{fdr.niche},
#'   \code{DirectionNiche}, \code{bandwidth.max}). \code{@fits} is empty (no
#'   per-bandwidth GLM fit exists for this estimator). Diagnostics are
#'   attached at \code{r@diagnostics}, a list of three tables: \code{r2} (the
#'   niche columns' R2 against the SpaNorm biology AND ls bases -- a `basis`
#'   column distinguishes them -- per sample x index, "spanorm" stage1 only;
#'   high biology overlap is expected, high LS overlap warns that the LS field
#'   may absorb the niche signal), \code{inclusion} (the per-index patient inclusion
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
                          epsilon = NULL,
                          winsor = 4, lambda.a = 0, maxit.psi = 2,
                          pool.psi = TRUE,
                          backend = "cpu", verbose = TRUE) {
  stage1 <- match.arg(stage1)
  if (!is.null(epsilon)) {
    warning("'epsilon' is deprecated and ignored: stage1 = \"spanorm\" now ",
            "fits an NB GLM with the SpaNorm ls+batch predictor as a fixed ",
            "offset, so there is no working response to add back to")
  }
  if (length(sigma) != 1L || !is.finite(sigma)) {
    stop("twoStageSpiDE() takes a single bandwidth (stage-1 slopes are per ",
         "bandwidth); got sigma of length ", length(sigma),
         ". Fit one bandwidth at a time -- spiDE()/fitSpiDE() are the ",
         "multi-bandwidth path.")
  }
  checkFdr(fdr)
  if (length(min.cells) != 1L || is.na(min.cells) || min.cells < 1) {
    stop("min.cells must be a single number >= 1")
  }
  checkSPE(spe, assay = assay, cell_type = cell_type, sample_id = sample_id)
  checkCondition(spe, condition)
  checkNiche(spe, sigma, name = name)
  checkCovariates(spe, patient.covariates)
  cd <- SummarizedExperiment::colData(spe)
  if (!is.null(patient) && !patient %in% colnames(cd)) {
    stop(sprintf("patient column '%s' not found in colData(spe)", patient))
  }
  smp <- as.character(cd[[sample_id]])
  pat <- if (is.null(patient)) smp else as.character(cd[[patient]])
  ct <- as.character(cd[[cell_type]])
  grp_cell <- as.character(cd[[condition]])
  chk <- tapply(grp_cell, pat, function(z) length(unique(z[!is.na(z)])))
  if (any(chk > 1)) {
    stop("condition '", condition, "' varies within patient; two-stage needs ",
         "a patient-level condition")
  }
  if (length(patient.covariates)) {
    # A patient-level covariate that actually varies within patient would
    # otherwise reach tapply(cd[[v]], pat, function(z) z[1]) below, which
    # silently keeps only the first cell's value per patient -- an arbitrary,
    # wrong covariate value rather than an error. Mirror the condition check
    # above.
    varies <- vapply(patient.covariates, function(v) {
      n_lvl <- tapply(cd[[v]], pat, function(z) length(unique(z[!is.na(z)])))
      any(n_lvl > 1)
    }, logical(1))
    if (any(varies)) {
      stop("patient.covariate(s) varies within patient: ",
           paste(patient.covariates[varies], collapse = ", "),
           "; two-stage needs patient-level covariates")
    }
  }
  # the precision pooling below attributes each core wholly to one patient;
  # a sample spanning patients (a colData merge error) would silently hand
  # the whole core to whichever patient its first cell carries
  n_pat_per_smp <- tapply(pat, smp, function(z) length(unique(z)))
  if (any(n_pat_per_smp > 1)) {
    stop("sample(s) span more than one patient: ",
         paste(names(n_pat_per_smp)[n_pat_per_smp > 1], collapse = ", "),
         "; each sample (core) must belong to exactly one patient")
  }
  s2p <- stats::setNames(vapply(unique(smp),
                                function(s) pat[match(s, smp)], character(1)),
                         unique(smp))
  nm <- as.matrix(SingleCellExperiment::reducedDim(spe, paste0(name, sigma)))
  if (!is.null(niche)) {
    keep <- intersect(niche, colnames(nm))
    if (!length(keep)) stop("no requested niche cell types found")
    miss <- setdiff(niche, keep)
    if (length(miss)) {
      # stage-1 slopes are JOINT across niche columns, so silently dropping
      # a misspelled niche changes every remaining slope, not just the
      # reported rows -- surface it
      warning("requested niche cell type(s) not found, dropped: ",
              paste(miss, collapse = ", "), call. = FALSE)
    }
    nm <- nm[, keep, drop = FALSE]
  }
  idx_types <- sort(unique(ct))
  if (!is.null(index)) {
    found <- intersect(index, idx_types)
    miss <- setdiff(index, found)
    if (length(miss)) {
      warning("requested index cell type(s) not found, dropped: ",
              paste(miss, collapse = ", "), call. = FALSE)
    }
    idx_types <- found
  }
  if (!length(idx_types)) stop("no requested index cell types found")

  # Y is left as assay() returns it (dense, sparse, or DelayedArray) -- at
  # panel scale a full as.matrix() here costs ~8 GB (13k genes x 77k cells).
  # Only "ols" needs a fully dense working response (E); "spanorm" and "nb"
  # densify per (sample, index) subset only, inside fitNB().
  Y <- SummarizedExperiment::assay(spe, assay)
  if (stage1 %in% c("spanorm", "nb")) {
    if (!.looksLikeCounts(Y)) {
      stop("stage1 = '", stage1, "' needs counts, but assay '", assay,
           "' does not look like counts")
    }
    # .looksLikeCounts() only samples a compact block for speed and can miss
    # a negative value elsewhere; checkCounts() scans the full matrix via a
    # single Matrix-/DelayedArray-aware min() reduction (no densification --
    # see the checkCounts() source) to catch it anywhere.
    checkCounts(Y)
  }
  comp <- if (stage1 == "spanorm") .spanormComponents(spe) else NULL
  E <- NULL
  if (stage1 == "ols") {
    Yd <- as.matrix(Y)                  # ols is the one path that needs it
    lib <- colSums(Yd)
    E <- log1p(sweep(Yd, 2, mean(pmax(lib, 1)) / pmax(lib, 1), "*"))
  }
  if (verbose) {
    message(sprintf("stage 1 (%s): %d genes x %d index types x %d niches x %d samples",
                    stage1, nrow(Y), length(idx_types), ncol(nm),
                    length(unique(smp))))
  }
  sl <- .sampleSlopes(Y, E, comp, nm, ct, smp, idx_types, min.cells,
                      stage1 = stage1, winsor = winsor,
                      lambda.a = lambda.a, maxit.psi = maxit.psi,
                      pool.psi = pool.psi, backend = backend,
                      verbose = verbose)

  pats <- unique(unname(s2p))
  # condition per patient, first non-missing cell. checkCondition() allows NA
  # condition values; an NA patient would otherwise reach model.matrix(),
  # which silently drops that row and misaligns Xdes against the slope
  # matrices (the same failure mode handled for patient.covariates below).
  cond_pat <- .patientValue(grp_cell, pat, pats)
  if (anyNA(cond_pat)) {
    dropped <- pats[is.na(cond_pat)]
    if (verbose) {
      message(sprintf("stage 2: dropping %d patient(s) with missing '%s': %s",
                      length(dropped), condition,
                      paste(dropped, collapse = ", ")))
    }
    pats <- pats[!is.na(cond_pat)]
    cond_pat <- cond_pat[pats]
    lvl <- unique(grp_cell[!is.na(grp_cell)])
    if (any(table(factor(cond_pat, levels = lvl)) < 2)) {
      stop("fewer than 2 patients remain in a condition after dropping ",
           "patient(s) with missing '", condition, "'")
    }
  }
  cv_list <- NULL
  if (length(patient.covariates)) {
    # NA in a patient-level covariate would otherwise reach model.matrix(),
    # which silently drops that row -- misaligning Xdes against the (still
    # full-length) B/V slope matrices and crashing limma downstream with an
    # opaque "(subscript) logical subscript too long". Drop such patients
    # from pats up front instead: every later use of `pats` (pgrp, Xdes, and
    # the pooled$beta/var[, , pats] slices in the loop below) then stays in
    # sync automatically.
    cv_list <- lapply(patient.covariates, function(v)
      .patientValue(cd[[v]], pat, pats))
    bad <- Reduce(`|`, lapply(cv_list, is.na))
    if (any(bad)) {
      dropped <- pats[bad]
      if (verbose) {
        message(sprintf(
          "stage 2: dropping %d patient(s) with missing patient.covariates (%s): %s",
          length(dropped), paste(patient.covariates, collapse = ", "),
          paste(dropped, collapse = ", ")))
      }
      pats <- pats[!bad]
      cv_list <- lapply(cv_list, `[`, !bad)
      lvl <- unique(grp_cell[!is.na(grp_cell)])
      remaining <- table(factor(cond_pat[pats], levels = lvl))
      if (any(remaining < 2)) {
        stop("fewer than 2 patients remain in a condition after dropping ",
             "patient(s) with missing patient.covariates (",
             paste(patient.covariates, collapse = ", "),
             "); cannot fit the stage-2 contrast")
      }
    }
  }
  pgrp <- factor(stats::setNames(unname(cond_pat[pats]), pats))
  Xdes <- if (length(patient.covariates)) {
    cv <- as.data.frame(stats::setNames(cv_list, patient.covariates),
                        stringsAsFactors = TRUE)
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
  out <- if (length(recs)) {
    o <- do.call(rbind, recs)
    o$fdr.niche <- stats::p.adjust(o$p.niche, "BH")
    o$DirectionNiche <- ifelse(o$t > 0, "Up", "Down")
    o$bandwidth.max <- sigma
    o <- o[order(o$p.niche), , drop = FALSE]
    rownames(o) <- NULL
    o
  } else {
    empty
  }
  if (verbose) {
    if (nrow(out)) {
      message(sprintf("stage 2: %d triplets tested, %d at FDR %.2g",
                      nrow(out), sum(out$fdr.niche <= fdr), fdr))
    } else {
      message("stage 2: no estimable triplets")
    }
  }
  new("SpiDEResults", fits = list(), sigma = sigma, condition = condition,
      mode = "condition", index = idx_types, niche = colnames(nm),
      covariates = patient.covariates, coldata = cd, gene.weights = NULL,
      p.cauchy.pos = NULL, p.cauchy.neg = NULL,
      results = out[out$fdr.niche <= fdr, , drop = FALSE],
      fdr = fdr, call = match.call(), diagnostics = diagnostics)
}
