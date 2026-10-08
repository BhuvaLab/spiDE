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

# Label-free DerSimonian-Laird between-patient variance per (gene, niche).
.dlTau2 <- function(b, v) {
  w <- 1 / v
  ok <- is.finite(b) & is.finite(w) & w > 0
  w[!ok] <- 0; b[!ok] <- 0
  sw <- apply(w, c(1, 3), sum); sw2 <- apply(w^2, c(1, 3), sum)
  bbar <- apply(w * b, c(1, 3), sum) / sw
  m <- apply(ok, c(1, 3), sum)
  dev <- sweep(b, c(1, 3), bbar)
  Q <- apply(w * dev^2, c(1, 3), sum)
  # pmax keeps the attributes of its FIRST argument: the matrix goes first
  tau2 <- pmax((Q - (m - 1)) / (sw - sw2 / sw), 0)
  tau2[!is.finite(tau2)] <- 0
  tau2
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
# the df capped by the Kish effective patients. Research arms of
# feature/small-sample-df (research/smalldf/README.md) refer the same t to other
# df (df_rule): "proportional" = limma's df.total x n_eff / m (m usable
# patients); "kish_plus_prior" = min(df.total, n_eff - 1 + d0), the Kish
# residual df plus eBayes' prior df; "limma" = df.total (reference only).
.pooledColumnTest <- function(b, v, tau2, mean_expr, min.pooled = 6L, trend = TRUE,
                              df_rule = c("capped", "proportional", "kish_plus_prior", "limma")) {
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
      bf <- b[fitted, , drop = FALSE]
      df[fitted] <- switch(df_rule,
        capped = pmin(fit$df.total, .kishDf(w, bf, 1L)),
        proportional = {
          ww <- w; ww[!is.finite(ww) | !is.finite(bf)] <- NA
          neff <- rowSums(ww, na.rm = TRUE)^2 / rowSums(ww^2, na.rm = TRUE)
          fit$df.total * neff / rowSums(is.finite(ww))
        },
        kish_plus_prior = pmin(fit$df.total, .kishDf(w, bf, 1L) + fit$df.prior),
        limma = fit$df.total)
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
# guarded and a gene that cannot be tested is NA.
.robustConditionTest <- function(b, v, tau2, trt, strata = NULL, min.group = 3L,
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
    out <- matrix(NA_real_, length(gi), 3L)
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
      w <- w[ok] / mean(w[ok])
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
      out[a, ] <- c(coef[2], sqrt(V), nu)
    }
    list(idx = gi, out = out)
  }
  res <- .bpGenes(G, one, BPPARAM)
  M <- matrix(NA_real_, G, 3L)
  for (r in res) M[r$idx, ] <- r$out
  tt <- M[, 1] / M[, 2]
  data.frame(estimate = M[, 1], se = M[, 2], t = tt, df = M[, 3], p = 2 * stats::pt(-abs(tt), M[, 3]),
             n_patients = rowSums(is.finite(b)), stringsAsFactors = FALSE)
}

# The test-stage arm of feature/small-sample-df: the pooled test's df rule
# (.pooledColumnTest()). The default is the shipped test.
.armSpec <- function(pooled_df = c("capped", "proportional", "kish_plus_prior", "limma")) {
  list(pooled_df = match.arg(pooled_df))
}

# The slopes engine's tests for one index type's stored slopes.
.slopesTests <- function(x, trt = NULL, strata = NULL, BPPARAM = BiocParallel::SerialParam(),
                         arm = .armSpec()) {
  stopifnot(is.list(arm), "pooled_df" %in% names(arm))
  out <- list()
  vpool <- sweep(x$v_model, c(2, 3), x$factor, "*")
  tau_pool <- .dlTau2(x$beta, vpool)
  for (j in seq_along(x$niches)) {
    b <- matrix(x$beta[, , j], nrow = length(x$genes))
    vj <- matrix(vpool[, , j], nrow = length(x$genes))
    p <- .pooledColumnTest(b, vj, tau_pool[, j], x$mean_expr, df_rule = arm$pooled_df)
    out[[length(out) + 1L]] <- data.frame(gene = x$genes, niche = x$niches[j], test = "pooled", p,
                                          stringsAsFactors = FALSE)
    if (!is.null(trt)) {
      cc <- .robustConditionTest(b, vj, tau_pool[, j], trt = trt, strata = strata, BPPARAM = BPPARAM)
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
