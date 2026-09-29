#!/usr/bin/env Rscript
# Calibration + validity diagnostic for spiDE (>= 0.99.30) test results.
#
#   Rscript calibration_check.R [--expr fit.rds] real.rds [...] \
#           [--block grid01.rds ...] [--perm perm001.rds ...]
#
# Inputs: .rds files each holding a SpiDEResults (testSpiDE() / spiDE()) or a
# data.frame from results(res, test = "both") (or a list with it under
# $results / $res). Files before any flag are REAL runs, each scored on its
# own. Files after --block (niche field shifted toroidally per section) or
# --perm (condition labels permuted across patients within slide) are NULL
# grids of one configuration, scored together; --null is a generic null set
# (no pass/fail limits). --expr names a SpiDEFit or SpiDEResults whose
# per-index mean expression defines the bands for data.frame inputs (the mean
# expression is the same for every block and perm grid of one fit).
#
# Prints, per set:
#   (a) global spread of z per test: sd, RMS, frac |z| > 1.96
#       z = sign(t) * qnorm(p / 2, lower.tail = FALSE), a standard normal under
#       the null whatever the test's df
#   (b) per-band RMS of z: RMS per (gene, index) over niches (and over grids),
#       median within each quintile of mean expression within index type
#   (c) per-index breakdown of (a)
#   (d) null sets only: fraction of grids with >= 1 BH call at 0.05, per test /
#       family (pooled; condition over every triplet; condition over the
#       filtered family), plus real-vs-null call counts when both are given
#   (e) real runs only: condition dropout -- per index type, patients per
#       condition carrying a usable slope, and a Fisher test of inclusion vs
#       condition
#
# Limits (research/simplify, pre-registered; vignettes/spiDE-calibration.Rmd):
#   bands     median RMS of null z in [0.90, 1.10] in every band: condition
#             test on both nulls, pooled test on the block null
#   block     fraction of grids with >= 1 BH call at 0.05 <= 0.10 (both tests)
#   perm      the same fraction <= 0.07, condition test over every triplet and
#             over the filtered family (the pooled test never sees the labels)
#
# Reads slots with attr() and never loads spiDE: a stray is()/inherits() on
# the S4 object makes R load whichever spiDE the default library holds (an old
# mixed-model build on this cluster).
suppressPackageStartupMessages(library(stats))
BAND <- c(0.90, 1.10); TAIL <- c(block = 0.10, perm = 0.07); ALPHA <- 0.05

# ---- arguments ----------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  stop("usage: calibration_check.R [--expr fit.rds] real.rds ... [--block g.rds ...] [--perm p.rds ...]")
}
sets <- list(real = character(), block = character(), perm = character(), null = character())
mode <- "real"; expr_file <- NULL; i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a %in% c("--real", "--block", "--perm", "--null")) {
    mode <- sub("^--", "", a)
  } else if (a == "--expr") {
    i <- i + 1L; expr_file <- args[i]
  } else {
    if (!file.exists(a)) stop("no such file: ", a)
    sets[[mode]] <- c(sets[[mode]], a)
  }
  i <- i + 1L
}

# ---- reading ------------------------------------------------------------------
.a <- function(o, n) attr(o, n, exact = TRUE)
.cls <- function(o) if (isS4(o)) .a(o, "class")[1] else class(o)[1]
.legacy <- function(o) {
  an <- names(attributes(o))
  any(c("covtype", "coefmap", "fits", "results.celltype") %in% an) || !any(c("engine", "table") %in% an)  # legacy slots
}
.legacyStop <- function(path) {
  stop(path, " is a spiDE <= 0.99.22 mixed-model object. Read it with ",
       "spiDEmixed::readSpiDE() (research/mixed) and score it with the archived ",
       "research/mixed/claude/skills/calibration-check/scripts/calibration_check.R.", call. = FALSE)
}

loadOne <- function(path) {
  o <- readRDS(path)
  if (!isS4(o) && is.list(o) && !is.data.frame(o)) {
    hit <- intersect(c("results", "res"), names(o))
    if (!length(hit)) stop(path, ": a list without $results or $res")
    o <- o[[hit[1]]]
  }
  if (isS4(o)) {
    cl <- .cls(o)
    if (.legacy(o)) .legacyStop(path)
    if (cl == "SpiDEFit") stop(path, ": a SpiDEFit carries no tests; run testSpiDE() on it (or pass it with --expr)")
    if (cl != "SpiDEResults") stop(path, ": unexpected S4 class ", cl)
    fit <- .a(o, "fit")
    if (.legacy(fit)) .legacyStop(path)
    return(list(path = path, tab = .a(o, "table"), fit = fit, engine = .a(fit, "engine"),
                condition = .a(o, "condition"), contrast = .a(o, "contrast"),
                procedure = .a(o, "procedure"), fdr = .a(o, "fdr")))
  }
  if (!is.data.frame(o)) stop(path, ": neither a SpiDEResults nor a data.frame")
  need <- c("gene", "index", "niche", "test", "t", "p")
  miss <- setdiff(need, colnames(o))
  if (length(miss)) {
    stop(path, ": lacks column(s) ", paste(miss, collapse = ", "),
         if (any(c("ct_index", "p.niche") %in% colnames(o))) " (a spiDE <= 0.99.22 table: use the archived script)")
  }
  list(path = path, tab = o, fit = NULL, engine = NA_character_, condition = character(),
       contrast = character(), procedure = NA_character_, fdr = ALPHA)
}

bandsOf <- function(fit) {
  if (is.null(fit)) return(NULL)
  idx <- .a(fit, "index")
  do.call(rbind, lapply(names(idx), function(k) {
    me <- idx[[k]]$mean_expr
    if (is.null(me) || !length(me)) return(NULL)
    q <- rank(me, ties.method = "average") / length(me)
    data.frame(gene = idx[[k]]$genes, index = k, mean_expr = unname(me),
               band = cut(q, c(0, .2, .4, .6, .8, 1), labels = paste0("Q", 1:5)),
               stringsAsFactors = FALSE)
  }))
}

addZ <- function(tb) {
  tb$z <- ifelse(is.finite(tb$p) & is.finite(tb$t),
                 sign(tb$t) * qnorm(pmax(tb$p, 1e-300) / 2, lower.tail = FALSE), NA_real_)
  tb
}

# ---- summaries ----------------------------------------------------------------
spread <- function(z) {
  z <- z[is.finite(z)]
  data.frame(n = length(z), sd_z = if (length(z) > 1) sd(z) else NA_real_,
             rms_z = sqrt(mean(z^2)), frac_1.96 = mean(abs(z) > 1.96))
}
byGroup <- function(tb, keys) {
  g <- split(tb, tb[keys], drop = TRUE)
  do.call(rbind, lapply(g, function(d) cbind(d[1, keys, drop = FALSE], spread(d$z))))
}
flag <- function(v, lim) ifelse(is.na(v), "", ifelse(v >= lim[1] & v <= lim[2], "", "*"))

# the families of .bhFamilies(), recomputed from p so any procedure/engine scores alike
families <- function(tb, filt = ALPHA) {
  pl <- tb$test == "pooled" & is.finite(tb$p)
  cl <- tb$test == "condition" & is.finite(tb$p)
  qp <- rep(NA_real_, nrow(tb)); qp[pl] <- p.adjust(tb$p[pl], "BH")
  key <- paste(tb$gene, tb$index, tb$niche, sep = "\r")
  pass <- key[pl][qp[pl] < filt]
  qa <- rep(NA_real_, nrow(tb)); qa[cl] <- p.adjust(tb$p[cl], "BH")
  cf <- cl & key %in% pass
  qf <- rep(NA_real_, nrow(tb)); qf[cf] <- p.adjust(tb$p[cf], "BH")
  out <- c(pooled = sum(qp <= ALPHA, na.rm = TRUE),
           condition_all = if (any(cl)) sum(qa <= ALPHA, na.rm = TRUE) else NA,
           condition_filtered = if (any(cl)) sum(qf <= ALPHA, na.rm = TRUE) else NA)
  # cross-check the recomputed filtered family against the stored in_family
  # (only when the run used the filtered procedure: then not every condition
  # row with a finite p is in the family)
  if ("in_family" %in% names(tb) && any(cl)) {
    inf <- tb$test == "condition" & tb$in_family %in% TRUE
    if (!all(inf[cl])) attr(out, "mismatch") <- sum(inf[cl] != cf[cl])
  }
  out
}
checkFamily <- function(fc, label) {
  m <- attr(fc, "mismatch")
  if (!is.null(m) && m > 0) {
    cat(sprintf("   WARNING %s: %d condition row(s) differ between the recomputed filtered family and\n", label, m))
    cat("   the stored in_family (a filter level other than the fdr assumed, or a table not from testSpiDE())\n")
  }
}

bandTable <- function(runs, bands, scored_tests) {
  if (is.null(bands)) {
    cat("   skipped: no per-index mean expression (give a SpiDEResults, or --expr fit.rds)\n")
    return(invisible(NULL))
  }
  # RMS per (gene, index, test) within each grid over niches, then over grids
  per <- do.call(rbind, lapply(seq_along(runs), function(r) {
    tb <- runs[[r]]$tab; tb <- tb[is.finite(tb$z), ]
    if (!nrow(tb)) return(NULL)
    ag <- aggregate(z ~ gene + index + test, tb, function(z) mean(z^2))
    ag$grid <- r; ag
  }))
  ms <- aggregate(z ~ gene + index + test, per, mean)
  ms$rms <- sqrt(ms$z)
  ms <- merge(ms, bands, by = c("gene", "index"))
  if (!nrow(ms)) { cat("   skipped: genes do not match the fit's\n"); return(invisible(NULL)) }
  med <- aggregate(rms ~ test + band, ms, median)
  w <- reshape(med, idvar = "test", timevar = "band", direction = "wide")
  names(w) <- sub("^rms\\.", "", names(w))
  bcols <- intersect(paste0("Q", 1:5), names(w))
  out <- w[, c("test", bcols)]
  out$range <- sprintf("%.3f-%.3f", apply(out[bcols], 1, min, na.rm = TRUE), apply(out[bcols], 1, max, na.rm = TRUE))
  out$out_of_band <- apply(out[bcols], 1, function(v) sum(v < BAND[1] | v > BAND[2], na.rm = TRUE))
  out[bcols] <- lapply(out[bcols], function(v) sprintf("%.3f", v))
  print(out, row.names = FALSE)
  if (!is.null(scored_tests)) {
    for (tt in intersect(scored_tests, out$test)) {
      ok <- out$out_of_band[out$test == tt] == 0
      cat(sprintf("   %-9s bands in [%.2f, %.2f]: %s\n", tt, BAND[1], BAND[2], if (ok) "PASS" else "FAIL"))
    }
  }
  invisible(out)
}

dropout <- function(run) {
  fit <- run$fit; cond <- run$condition
  if (is.null(fit)) { cat("   skipped: a data.frame carries no patients (give the SpiDEResults)\n"); return(invisible()) }
  if (!length(cond)) { cat("   skipped: no condition was tested\n"); return(invisible()) }
  pt <- .a(fit, "patients")
  if (!cond %in% names(pt)) { cat(sprintf("   skipped: '%s' is not a patient-level column of the fit\n", cond)); return(invisible()) }
  lab <- setNames(as.character(pt[[cond]]), pt$patient)
  lv <- if (is.factor(pt[[cond]])) levels(droplevels(pt[[cond]])) else sort(unique(lab[!is.na(lab)]))
  idx <- .a(fit, "index"); eng <- .a(fit, "engine")
  cat(if (eng == "slopes") {
    "   usable = the patient enters the index type (>= min.cells cells) and has a finite slope\n   for at least half of its (gene, niche) pairs (slopes engine)\n"
  } else {
    "   usable = the patient enters the index type (>= min.cells cells); the sandwich engine\n   stores no per-patient slopes, so per-gene dropout is not visible here\n"
  })
  rows <- lapply(names(idx), function(k) {
    x <- idx[[k]]
    use <- setNames(pt$patient %in% x$patients, pt$patient)
    if (eng == "slopes" && !is.null(x$beta)) {
      fr <- apply(is.finite(x$beta), 2, mean)          # per patient, over genes x niches
      use[x$patients] <- fr >= 0.5
    }
    tt <- table(factor(lab, levels = lv), factor(use[names(lab)], levels = c(FALSE, TRUE)))
    p <- if (all(dim(tt) == 2) && all(rowSums(tt) > 0) && all(colSums(tt) > 0)) fisher.test(tt)$p.value else NA_real_
    data.frame(index = k, a = sprintf("%d/%d", tt[1, 2], sum(tt[1, ])),
               b = sprintf("%d/%d", tt[2, 2], sum(tt[2, ])), fisher_p = p, stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, rows)
  names(d)[2:3] <- paste("usable", lv)
  print(d, row.names = FALSE, digits = 3)
  bad <- d$index[!is.na(d$fisher_p) & d$fisher_p < 0.05]
  if (length(bad)) {
    cat(sprintf("   CONFOUNDED: %s -- which patients contribute depends on the condition;\n", paste(bad, collapse = ", ")))
    cat("   no threshold or variance correction fixes informative missingness.\n")
  } else cat("   no index type with inclusion associated with the condition (Fisher p < 0.05)\n")
}

describe <- function(r) {
  if (is.null(r$fit)) return(sprintf("a data.frame from results(), %d rows", nrow(r$tab)))
  sprintf("%s engine, condition %s%s, procedure %s, %d rows",
          if (is.na(r$engine)) "unknown" else r$engine,
          if (length(r$condition)) sprintf("'%s'", r$condition) else "none",
          if (length(r$contrast) && nzchar(r$contrast)) sprintf(" (%s)", r$contrast) else "",
          if (is.na(r$procedure)) "unknown" else r$procedure, nrow(r$tab))
}

scoreSet <- function(runs, kind, expr_bands) {
  if (kind == "perm") {
    # every permutation grid repeats the real pooled test (it never sees the
    # labels), so its rows are not null draws: score the condition test only
    runs <- lapply(runs, function(r) { r$tab <- r$tab[r$tab$test == "condition", , drop = FALSE]; r })
    cat("   (pooled rows left out: they are the real pooled test, identical in every permutation)\n\n")
  }
  tb <- do.call(rbind, lapply(runs, `[[`, "tab"))
  bands <- if (!is.null(runs[[1]]$fit)) bandsOf(runs[[1]]$fit) else expr_bands
  cat("(a) global spread of z per test\n")
  print(byGroup(tb, "test"), row.names = FALSE, digits = 3)
  if (kind == "real") cat("   (real data: signal inflates the spread; calibration is judged on the nulls)\n")
  scored <- switch(kind, block = c("pooled", "condition"), perm = "condition", NULL)
  cat("\n(b) per-band RMS of z (median over (gene, index) within expression quintiles of each index type)\n")
  bandTable(runs, bands, scored)
  cat("\n(c) per index\n")
  print(byGroup(tb, c("index", "test")), row.names = FALSE, digits = 3)
  invisible(tb)
}

# ---- run ----------------------------------------------------------------------
expr_bands <- NULL
if (!is.null(expr_file)) {
  e <- readRDS(expr_file)
  if (isS4(e) && .legacy(e)) .legacyStop(expr_file)
  if (isS4(e) && .cls(e) == "SpiDEResults") e <- .a(e, "fit")
  if (!isS4(e) || .cls(e) != "SpiDEFit") stop("--expr must name a SpiDEFit or SpiDEResults")
  expr_bands <- bandsOf(e)
}
if (is.null(expr_bands) && length(sets$real)) {
  first <- loadOne(sets$real[1]); expr_bands <- bandsOf(first$fit)
}

real_calls <- list()
for (f in sets$real) {
  r <- loadOne(f); r$tab <- addZ(r$tab)
  cat(sprintf("\n==== REAL: %s\n     %s\n\n", basename(f), describe(r)))
  scoreSet(list(r), "real", expr_bands)
  real_calls[[f]] <- families(r$tab, if (length(r$fdr)) r$fdr else ALPHA)
  checkFamily(real_calls[[f]], basename(f))
  cat(sprintf("\n   BH calls at %.2f: pooled %d | condition, every triplet %s | condition, filtered family %s\n",
              ALPHA, real_calls[[f]][1], real_calls[[f]][2], real_calls[[f]][3]))
  cat("\n(e) condition dropout (informative missingness)\n")
  dropout(r)
}

for (kind in c("block", "perm", "null")) {
  if (!length(sets[[kind]])) next
  runs <- lapply(sets[[kind]], function(f) { r <- loadOne(f); r$tab <- addZ(r$tab); r })
  cat(sprintf("\n==== NULL (%s): %d grid(s)\n     first: %s -- %s\n\n", kind, length(runs),
              basename(runs[[1]]$path), describe(runs[[1]])))
  scoreSet(runs, kind, expr_bands)
  fam <- lapply(runs, function(r) families(r$tab, if (length(r$fdr)) r$fdr else ALPHA))
  for (j in seq_along(fam)) checkFamily(fam[[j]], basename(runs[[j]]$path))
  calls <- do.call(rbind, fam)
  cat(sprintf("\n(d) tails: fraction of grids with >= 1 BH call at %.2f\n", ALPHA))
  lim <- if (kind %in% names(TAIL)) TAIL[[kind]] else NA
  fams <- colnames(calls)
  scored <- switch(kind, block = fams, perm = c("condition_all", "condition_filtered"), character())
  for (fm in fams) {
    v <- calls[, fm]
    if (all(is.na(v))) next
    if (kind == "perm" && fm == "pooled") {
      same <- length(unique(v)) == 1L
      cat(sprintf("   %-19s not scored: the pooled test never sees the labels%s\n", fm,
                  if (same) "" else " (BUT ITS CALLS DIFFER ACROSS PERMUTATIONS -- a label leak)"))
      next
    }
    frac <- mean(v > 0)
    verdict <- if (fm %in% scored && is.finite(lim)) sprintf("  limit %.2f: %s", lim, if (frac <= lim) "PASS" else "FAIL") else ""
    cat(sprintf("   %-19s %d/%d grids = %.3f  (mean calls %.2f, max %d)%s\n", fm, sum(v > 0), length(v),
                frac, mean(v), max(v), verdict))
  }
  if (length(runs) < 10L && kind %in% names(TAIL)) {
    cat(sprintf("   (only %d grids: a tail fraction this coarse cannot resolve the limit)\n", length(runs)))
  }
  if (length(real_calls)) {
    cat("\n   real vs this null (empirical FDP ~ mean null calls / real calls; meaningful only if\n   the grids are this real run's own nulls):\n")
    for (f in names(real_calls)) {
      rc <- real_calls[[f]]; nm <- colMeans(calls, na.rm = TRUE)
      for (fm in names(rc)) {
        if (is.na(rc[[fm]]) || (kind == "perm" && fm == "pooled")) next
        cat(sprintf("   %-19s real %4d | null mean %6.2f | FDP ~ %s\n", fm, rc[[fm]], nm[[fm]],
                    if (rc[[fm]] > 0) sprintf("%.3f", nm[[fm]] / rc[[fm]]) else "n/a"))
      }
    }
  }
}
cat("\nDecisions rest on the public cohorts; YTMA results enter at weight 0.25 (user directive 2026-09-29).\n")
