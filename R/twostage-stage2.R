# Stage-2 machinery: patient pooling of per-core slopes, the DL tau2, the
# weighted moderated limma contrast, and the dropout diagnostics.

#' Precision-weighted pooling of per-core slopes within patients
#'
#' Cores share their patient's slope, so fixed-effect (1/v) pooling is
#' correct here (the between-patient component enters later, in tau2).
#' @noRd
.poolPatientSlopes <- function(B, V, sample2patient) {
  pats <- unique(unname(sample2patient[dimnames(B)[[3]]]))
  db <- dimnames(B)
  Bp <- Vp <- array(NA_real_, c(dim(B)[1:2], length(pats)),
                    dimnames = list(db[[1]], db[[2]], pats))
  for (p in pats) {
    ss <- dimnames(B)[[3]][sample2patient[dimnames(B)[[3]]] == p]
    b <- B[, , ss, drop = FALSE]; v <- V[, , ss, drop = FALSE]
    wt <- 1 / v
    wt[is.na(b) | !is.finite(wt)] <- NA
    swt <- apply(wt, 1:2, sum, na.rm = TRUE)
    num <- apply(b * wt, 1:2, sum, na.rm = TRUE)
    est <- num / swt
    est[swt == 0] <- NA_real_
    Bp[, , p] <- est
    vp <- 1 / swt
    vp[swt == 0] <- NA_real_
    Vp[, , p] <- vp
  }
  list(beta = Bp, var = Vp)
}

#' Between-patient variance component, DerSimonian-Laird, pooled over genes
#'
#' Per gene: weighted (1/v) regression on the stage-2 design, Cochran's
#' Q from the weighted residuals, DL moment estimate
#' (Q - (S - p)) / (sum(w) - tr((X'WX)^{-1} X'W^2X)); then the median over
#' genes, floored at 0. One tau2 per (index, niche) -- matching how the
#' package shares strength across genes elsewhere.
#' @noRd
.tau2DL <- function(B, V, X) {
  p <- ncol(X)
  t2 <- apply(cbind(B, V), 1, function(row) {
    S <- length(row) / 2
    b <- row[seq_len(S)]; v <- row[S + seq_len(S)]
    ok <- is.finite(b) & is.finite(v) & v > 0
    if (sum(ok) < p + 2) return(NA_real_)
    w <- 1 / v[ok]; Xo <- X[ok, , drop = FALSE]; bo <- b[ok]
    XtW <- t(Xo * w)
    M <- tryCatch(solve(XtW %*% Xo), error = function(e) NULL)
    if (is.null(M)) return(NA_real_)
    res <- bo - Xo %*% (M %*% (XtW %*% bo))
    Q <- sum(w * res^2)
    denom <- sum(w) - sum(diag(M %*% (t(Xo * w^2) %*% Xo)))
    if (denom <= 0) return(NA_real_)
    max((Q - (sum(ok) - p)) / denom, 0)
  })
  med <- stats::median(t2, na.rm = TRUE)
  if (!is.finite(med)) 0 else med
}

#' Moderated condition contrast on the patient slope matrix
#'
#' Observation weights 1/(v + tau2): random-effects weighting that degrades
#' toward equal weights when between-patient variation dominates (the
#' pseudobulk pitfall guard) and does real work when subset precision varies.
#' NA slopes enter with weight 0 and value 0 (the standard limma idiom).
#' @importFrom limma lmFit eBayes
#' @noRd
.limmaStage2 <- function(B, V, tau2, X) {
  Wt <- 1 / (V + tau2)
  bad <- !is.finite(B) | !is.finite(Wt)
  Bf <- B; Bf[bad] <- 0; Wt[bad] <- 0
  fit <- limma::lmFit(Bf, design = X, weights = Wt)
  fit <- limma::eBayes(fit, robust = TRUE)
  cn <- colnames(X)[2]
  data.frame(gene = rownames(B),
             coef = fit$coefficients[, cn],
             t = fit$t[, cn],
             p = fit$p.value[, cn],
             row.names = NULL, stringsAsFactors = FALSE)
}

#' Per-index patient inclusion table, with a condition-association warning
#'
#' min.cells dropout is informative when cell-type abundance correlates with
#' condition; this makes it visible instead of silent.
#' @noRd
.inclusionDiagnostics <- function(ncells, sample2patient, pgrp, min.cells) {
  ncells$patient <- unname(sample2patient[ncells$sample])
  agg <- stats::aggregate(n ~ index + patient, ncells, max)
  agg$included <- agg$n >= min.cells
  agg$condition <- as.character(pgrp[agg$patient])
  out <- agg[, c("index", "patient", "condition", "included")]
  for (ix in unique(out$index)) {
    d <- out[out$index == ix, ]
    if (length(unique(d$condition)) < 2 || all(d$included) || !any(d$included)) next
    p <- stats::fisher.test(table(d$condition, d$included))$p.value
    if (p < 0.05) {
      warning("min.cells dropout for index '", ix, "' is associated with ",
              "the condition (Fisher p = ", signif(p, 2), "); the tested ",
              "patient subset is condition-biased", call. = FALSE)
    }
  }
  out
}
