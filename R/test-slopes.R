# The slopes engine's tests across patients.
#
# For every (index, tested niche) the per-patient slopes beta_s of each gene
# are combined across patients, weighted by 1 / (v_model_s x f_s + tau2_g):
#   the POOLED test   ~ 1                 (does the slope differ from zero?)
#   the CONDITION test ~ condition [+ strata] (does it differ between conditions?)
# The pooled test is a moderated limma test (below). The condition test is the
# same weighted least squares with an HC2 standard error across patients and
# Bell-McCaffrey df (.robustConditionTest(), 0.99.32). Its weights do not move
# with any one gene's departure from the pooled slope (the per-patient factor
# f_s is a median over genes), so the estimate is not pulled toward the pooled
# slope; and the robust SE lets a patient whose slope is noisier than its weight
# implies (a spatially patterned marker gene) widen the interval rather than
# make a call. To 0.99.31 the condition test was weighted by the gene's own tile
# variance, which attenuated the contrast (research/bench2/diag/FINDINGS.md,
# research/release/README.md).
# tau2_g is the gene's DerSimonian-Laird between-patient variance of the slope
# (label-free), the term the mixed model lacked. Two further choices make the
# pooled test calibrated on every null of the five-cohort study
# (research/simplify, 2026-09-29, "effective df + trended prior"):
#   * an expression-trended eBayes prior (eBayes(trend = TRUE) on the gene's log
#     mean expression in the index type), so bright genes are not squeezed toward
#     the typical gene; and
#   * the t statistic referred to min(limma's df.total, n_eff - k), where
#     n_eff = (sum w)^2 / sum w^2 is the Kish effective number of patients
#     behind the weights: when tau2_g truncates to 0 the weights can rest on a
#     few cell-rich patients, and limma's ~S df then overstate the evidence.
# Both are computed over the whole gene family of a column (the eBayes prior
# and trend are shared across genes), never within a gene block.

# Label-free DerSimonian-Laird moments per (gene, niche): the untruncated
# between-patient variance t = (Q - (m - 1)) / c, its numerator Qm and
# denominator c, and the usable patients m.
.dlMoments <- function(b, v) {
  w <- 1 / v
  ok <- is.finite(b) & is.finite(w) & w > 0
  w[!ok] <- 0; b[!ok] <- 0
  sw <- apply(w, c(1, 3), sum); sw2 <- apply(w^2, c(1, 3), sum)
  bbar <- apply(w * b, c(1, 3), sum) / sw
  m <- apply(ok, c(1, 3), sum)
  dev <- sweep(b, c(1, 3), bbar)
  Q <- apply(w * dev^2, c(1, 3), sum)
  cc <- sw - sw2 / sw
  list(t = (Q - (m - 1)) / cc, Qm = Q - (m - 1), c = cc, m = m)
}

# Label-free DerSimonian-Laird between-patient variance per (gene, niche).
.dlTau2 <- function(b, v) {
  # pmax keeps the attributes of its FIRST argument: the matrix goes first
  tau2 <- pmax(.dlMoments(b, v)$t, 0)
  tau2[!is.finite(tau2)] <- 0
  tau2
}

# The heterogeneity floor (branch feature/shared-tau2, research/sharedtau): per
# (gene, niche) the larger of the gene's own DerSimonian-Laird tau2 and its
# family's, a degree-1 loess of the untruncated DL moments on log mean
# expression over the column's genes, weighted by c and floored at 0 (locally
# the DL estimator pooled over genes, sum(Q - (m - 1)) / sum(c)). Under
# `min.rows` usable genes (spiGSEA's sets) the family value is that pooled
# ratio. A gene keeps its own heterogeneity where it is larger (markers,
# patterned genes) and its weights never fall back to 1 / v while its family is
# heterogeneous. Label-free, over the whole family of a column, never inside a
# gene block.
.heterogeneityFloor <- function(b, v, mean_expr, span = 0.5, min.rows = 20L) {
  mo <- .dlMoments(b, v)
  own <- pmax(mo$t, 0)
  own[!is.finite(own)] <- 0
  x <- log(mean_expr + 1e-3)
  for (j in seq_len(ncol(own))) {
    use <- which(is.finite(mo$t[, j]) & is.finite(mo$c[, j]) & mo$c[, j] > 0 & mo$m[, j] >= 2 & is.finite(x))
    own[, j] <- pmax(own[, j], .familyTrend(mo$t[, j], mo$c[, j], x, use, span, min.rows), 0)
  }
  own
}

# The family's heterogeneity at every gene's x: a c-weighted degree-1 loess of
# the untruncated moments t over the rows `use`, clamped to their x range (a
# gene outside it takes the edge value; a gene with unknown x the pooled
# ratio); under `min.rows` rows, or if the smoother fails, the pooled ratio
# sum(c t) / sum(c) everywhere.
.familyTrend <- function(t, cw, x, use, span = 0.5, min.rows = 20L) {
  fam <- rep(if (length(use)) sum(cw[use] * t[use]) / sum(cw[use]) else 0, length(t))
  if (length(use) >= min.rows && length(unique(x[use])) >= 5L) {
    d <- data.frame(t = t[use], x = x[use], w = cw[use])
    lo <- tryCatch(suppressWarnings(stats::loess(t ~ x, data = d, weights = w, span = span, degree = 1L,
                                                 control = stats::loess.control(surface = "direct"))),
                   error = function(e) NULL)
    if (!is.null(lo)) {
      pr <- tryCatch(suppressWarnings(stats::predict(lo, newdata = data.frame(x = pmin(pmax(x, min(d$x)), max(d$x))))),
                     error = function(e) NULL)
      if (!is.null(pr)) fam[is.finite(pr)] <- pr[is.finite(pr)]
    }
  }
  fam[!is.finite(fam)] <- 0
  fam
}

# The Kish effective number of patients behind each gene's weights, less the
# design's columns, floored at 1.
.kishDf <- function(wts, b, k) {
  ww <- wts
  ww[!is.finite(ww) | !is.finite(b)] <- NA
  neff <- rowSums(ww, na.rm = TRUE)^2 / rowSums(ww^2, na.rm = TRUE)
  pmax(neff - k, 1)
}

# The pooled test of one (index, niche) column: G x S slopes b and variances v,
# with tau2 (G): moderated limma on an intercept, the trended prior over the
# column's gene family (trend = FALSE: a flat prior, for spiGSEA's few sets),
# the df capped by the Kish effective patients (df_rule = "uncapped": limma's
# df.total, a research arm of feature/shared-tau2).
.pooledColumnTest <- function(b, v, tau2, mean_expr, min.pooled = 6L, trend = TRUE,
                              df_rule = c("capped", "uncapped")) {
  df_rule <- match.arg(df_rule)
  G <- nrow(b)
  wts <- 1 / (v + tau2)
  wts[!is.finite(wts) | wts <= 0 | !is.finite(b)] <- NA
  n_ok <- rowSums(is.finite(wts))
  est <- tt <- df <- rep(NA_real_, G)
  # rows with fewer than two usable patients carry no residual variance: limma
  # ignores them in the prior anyway, and in a small family (spiGSEA's sets)
  # an all-NA row can make the robust prior fail for every row
  fitted <- which(n_ok >= 2L)
  if (length(fitted)) {
    w <- wts[fitted, , drop = FALSE]
    w <- w / rowMeans(w, na.rm = TRUE)
    fit <- suppressWarnings(limma::lmFit(b[fitted, , drop = FALSE], matrix(1, ncol(b), 1), weights = w))
    fit$Amean <- log(mean_expr[fitted] + 1e-3)
    eb <- function(robust) tryCatch(suppressWarnings(limma::eBayes(fit, robust = robust, trend = trend)),
                                    error = function(e) NULL)
    fit <- eb(TRUE)
    if (is.null(fit)) fit <- eb(FALSE)
    if (!is.null(fit)) {
      est[fitted] <- fit$coefficients[, 1]
      tt[fitted] <- fit$t[, 1]
      df[fitted] <- if (df_rule == "capped") pmin(fit$df.total, .kishDf(w, b[fitted, , drop = FALSE], 1L)) else
        fit$df.total
    }
  }
  tt[n_ok < min.pooled | !is.finite(tt)] <- NA
  data.frame(estimate = est, se = est / tt, t = tt, df = df, p = 2 * stats::pt(-abs(tt), df),
             n_patients = unname(rowSums(is.finite(b))), stringsAsFactors = FALSE)
}

# The condition test of one (index, niche) column: per gene, weighted least
# squares of the patients' slopes b on [1, condition, strata] with weights
# 1 / (v + tau2), the condition coefficient's HC2 standard error across
# patients, and Bell-McCaffrey df under the working model var(b_s) ~ 1 / w_s.
# With K = X A X' (A = (X'WX)^-1) the residual covariance under that model is
# M W^-1 M' = W^-1 - K, so the df need only K. A patient with an unusable slope
# or weight drops out; a gene needs min.group effective (Kish) patients in each
# condition, so its fit cannot rest on one or two heavily weighted patients.
# Genes are independent and blocked over BPPARAM; every per-gene inversion is
# guarded and a gene that cannot be tested is NA. moderate = TRUE (a research
# arm of feature/shared-tau2) squeezes the HC2 variance across the column's
# genes (.conditionTable()).
.robustConditionTest <- function(b, v, tau2, trt, strata = NULL, min.group = 3L,
                                 BPPARAM = BiocParallel::SerialParam(), moderate = FALSE, mean_expr = NULL) {
  M <- .conditionPieces(b, v, tau2, trt, strata = strata, min.group = min.group, BPPARAM = BPPARAM)
  .conditionTable(M, rowSums(is.finite(b)), moderate = moderate, mean_expr = mean_expr)
}

# The per-gene pieces of the condition test: G x 4 (estimate, HC2 variance V,
# Bell-McCaffrey df nu, A22 = the condition coefficient's model variance
# [(X'WX)^-1]_22 under the UNNORMALISED weights 1 / (v + tau2)).
.conditionPieces <- function(b, v, tau2, trt, strata = NULL, min.group = 3L,
                             BPPARAM = BiocParallel::SerialParam()) {
  G <- nrow(b)
  X0 <- cbind(1, trt)
  if (!is.null(strata)) {
    sf <- droplevels(factor(strata))
    if (nlevels(sf) > 1L) {
      Sd <- stats::model.matrix(~ sf)[, -1, drop = FALSE]
      # confounded when the condition lies in the strata's span: pivoting would
      # otherwise drop a strata column and test the condition unadjusted
      S1 <- cbind(1, Sd)
      if (qr(cbind(S1, trt))$rank == qr(S1)$rank) {
        stop("the condition is fully confounded with 'strata'", call. = FALSE)
      }
      X0 <- cbind(X0, Sd)
    }
  }
  kish <- function(w) if (length(w)) sum(w)^2 / sum(w^2) else 0
  one <- function(gi) {
    out <- matrix(NA_real_, length(gi), 4L)
    for (a in seq_along(gi)) {
      g <- gi[a]
      if (!is.finite(tau2[g])) next
      w <- 1 / (v[g, ] + tau2[g])
      ok <- is.finite(b[g, ]) & is.finite(w) & w > 0
      if (kish(w[ok & trt == 1]) < min.group || kish(w[ok & trt == 0]) < min.group) next
      X <- X0[ok, , drop = FALSE]
      qx <- qr(X)
      keep <- sort(qx$pivot[seq_len(qx$rank)])
      if (!all(c(1L, 2L) %in% keep) || length(keep) >= nrow(X)) next
      X <- X[, keep, drop = FALSE]
      mw <- mean(w[ok])
      w <- w[ok] / mw
      y <- b[g, ok]
      A <- tryCatch(solve(crossprod(X * sqrt(w))), error = function(e) NULL)
      if (is.null(A) || !all(is.finite(A))) next
      XA <- X %*% A
      h <- rowSums(XA * X) * w                        # leverages
      if (any(h >= 1 - 1e-8)) next
      coef <- as.numeric(A %*% crossprod(X, w * y))
      e <- y - as.numeric(X %*% coef)
      d <- (XA[, 2] * w)^2 / (1 - h)
      V <- sum(d * e^2)
      Sig <- -tcrossprod(XA, X)                       # W^-1 - K
      diag(Sig) <- diag(Sig) + 1 / w
      nu <- sum(d * diag(Sig))^2 / sum(tcrossprod(d) * Sig^2)
      if (!is.finite(V) || V <= 0 || !is.finite(nu) || nu <= 0) next
      # A was computed with weights scaled to mean 1: (X' W X)^-1 = A / mw
      out[a, ] <- c(coef[2], V, nu, A[2, 2] / mw)
    }
    list(idx = gi, out = out)
  }
  res <- .bpGenes(G, one, BPPARAM)
  M <- matrix(NA_real_, G, 4L, dimnames = list(NULL, c("estimate", "V", "nu", "A22")))
  for (r in res) M[r$idx, ] <- r$out
  M
}

# The condition test's table from its pieces. HC2: se = sqrt(V) on nu df.
# moderate = TRUE: phi = V / A22 (about 1 under the working model; unitless) is
# squeezed with limma::squeezeVar(robust = TRUE, covariate = log mean
# expression) over the column's complete rows, after the gene blocks are
# reassembled; se = sqrt(phi_post * A22) on nu + df.prior df. Under `min.rows`
# complete rows (spiGSEA's sets) nothing is moderated.
.conditionTable <- function(M, n_patients, moderate = FALSE, mean_expr = NULL, min.rows = 20L) {
  est <- M[, 1]; se <- sqrt(M[, 2]); df <- M[, 3]
  done <- FALSE
  if (moderate) {
    stopifnot(length(mean_expr) == nrow(M))
    phi <- M[, 2] / M[, 4]
    x <- log(mean_expr + 1e-3)
    ok <- which(is.finite(phi) & phi > 0 & is.finite(M[, 3]) & M[, 3] > 0 & is.finite(x))
    if (length(ok) >= min.rows) {
      sq <- tryCatch(suppressWarnings(limma::squeezeVar(phi[ok], df = M[ok, 3], covariate = x[ok], robust = TRUE)),
                     error = function(e) NULL)
      if (!is.null(sq) && all(is.finite(sq$var.post)) && !anyNA(sq$df.prior)) {
        se[ok] <- sqrt(sq$var.post * M[ok, 4])
        df[ok] <- M[ok, 3] + rep_len(sq$df.prior, length(ok))
        done <- TRUE
      }
    }
  }
  tt <- est / se
  out <- data.frame(estimate = unname(est), se = unname(se), t = unname(tt), df = unname(df),
                    p = unname(2 * stats::pt(-abs(tt), df)), n_patients = n_patients, stringsAsFactors = FALSE)
  # whether a requested moderation happened (a small or failed family stays HC2)
  if (moderate) attr(out, "moderated") <- done
  out
}

# The test-stage arms of feature/shared-tau2 (research/sharedtau/README.md):
# the heterogeneity behind the weights ("own": each gene's DL tau2; "floor":
# .heterogeneityFloor(); "equal": every patient weight 1), the pooled test's df
# ("capped" by the Kish effective patients, or limma's "uncapped" df.total) and
# the condition test's variance ("hc2", or "moderated" across genes). The
# default is the shipped test.
.armSpec <- function(heterogeneity = c("own", "floor", "equal"), pooled_df = c("capped", "uncapped"),
                     condition_variance = c("hc2", "moderated")) {
  list(heterogeneity = match.arg(heterogeneity), pooled_df = match.arg(pooled_df),
       condition_variance = match.arg(condition_variance))
}

# The heterogeneity of an arm for one index type's stored slopes (G x nt), and
# the patient variances it goes with: "equal" sets v to 0 (unusable patients
# stay NA) and tau2 to 1, so every usable patient weighs 1.
.armWeights <- function(x, heterogeneity) {
  vpool <- sweep(x$v_model, c(2, 3), x$factor, "*")
  switch(heterogeneity,
         own = list(v = vpool, tau2 = .dlTau2(x$beta, vpool)),
         floor = list(v = vpool, tau2 = .heterogeneityFloor(x$beta, vpool, x$mean_expr)),
         equal = list(v = vpool * 0, tau2 = matrix(1, length(x$genes), length(x$niches))))
}

# The slopes engine's tests for one index type's stored slopes.
.slopesTests <- function(x, trt = NULL, strata = NULL, BPPARAM = BiocParallel::SerialParam(),
                         arm = .armSpec()) {
  stopifnot(is.list(arm), all(c("heterogeneity", "pooled_df", "condition_variance") %in% names(arm)))
  out <- list()
  aw <- .armWeights(x, arm$heterogeneity)
  vpool <- aw$v
  tau_pool <- aw$tau2
  for (j in seq_along(x$niches)) {
    b <- matrix(x$beta[, , j], nrow = length(x$genes))
    vj <- matrix(vpool[, , j], nrow = length(x$genes))
    p <- .pooledColumnTest(b, vj, tau_pool[, j], x$mean_expr, df_rule = arm$pooled_df)
    out[[length(out) + 1L]] <- data.frame(gene = x$genes, niche = x$niches[j], test = "pooled", p,
                                          stringsAsFactors = FALSE)
    if (!is.null(trt)) {
      cc <- .robustConditionTest(b, vj, tau_pool[, j], trt = trt, strata = strata, BPPARAM = BPPARAM,
                                 moderate = arm$condition_variance == "moderated", mean_expr = x$mean_expr)
      out[[length(out) + 1L]] <- data.frame(gene = x$genes, niche = x$niches[j], test = "condition", cc,
                                            stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

# Benjamini-Hochberg families. The pooled test: every (gene, index, niche).
# The condition test: every triplet ("all"), or only the triplets whose pooled
# test passes at `fdr` ("filtered", the default). The filter is label-free --
# it never looks at the condition -- so it is exact under any relabelling of
# patients, and it spends the multiplicity budget on triplets with a niche
# effect to modify.
.bhFamilies <- function(tab, procedure = c("filtered", "all"), fdr = 0.05) {
  procedure <- match.arg(procedure)
  tab$q <- NA_real_
  tab$in_family <- FALSE
  pl <- which(tab$test == "pooled" & is.finite(tab$p))
  tab$q[pl] <- stats::p.adjust(tab$p[pl], "BH")
  tab$in_family[pl] <- TRUE
  cl <- which(tab$test == "condition" & is.finite(tab$p))
  if (length(cl)) {
    if (procedure == "filtered") {
      key <- paste(tab$gene, tab$index, tab$niche, sep = "\r")
      pass <- key[pl][tab$q[pl] < fdr]
      cl <- cl[key[cl] %in% pass]
    }
    tab$q[cl] <- stats::p.adjust(tab$p[cl], "BH")
    tab$in_family[cl] <- TRUE
  }
  tab
}
