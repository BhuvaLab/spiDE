# The slopes engine's tests across patients.
#
# For every (index, tested niche) the per-patient slopes beta_s of each gene
# are combined across patients with limma, weighted by 1 / (v_s + tau2_g):
#   the POOLED test   ~ 1                 (does the slope differ from zero?)
#   the CONDITION test ~ condition [+ strata] (does it differ between conditions?)
# tau2_g is the gene's DerSimonian-Laird between-patient variance of the slope
# (label-free), the term the mixed model lacked. Two further choices make the
# tests calibrated on every null of the five-cohort study (research/simplify,
# 2026-09-29, "effective df + trended prior"):
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

# One (index, niche) column: G x S slopes b and variances v, with tau2 (G).
# trt = NULL gives the pooled test; otherwise a 0/1 vector over the S patients
# (and strata an optional patient-level factor).
.slopeColumnTest <- function(b, v, tau2, mean_expr, trt = NULL, strata = NULL,
                             min.pooled = 6L, min.group = 3L) {
  G <- nrow(b)
  if (is.null(trt)) {
    des <- matrix(1, ncol(b), 1)
  } else {
    des <- cbind(1, trt)
    if (!is.null(strata)) {
      sf <- droplevels(factor(strata))
      if (nlevels(sf) > 1L) {
        des <- cbind(des, stats::model.matrix(~ sf)[, -1, drop = FALSE])
        qd <- qr(des)
        des <- des[, sort(qd$pivot[seq_len(qd$rank)]), drop = FALSE]
        if (ncol(des) < 2L || !isTRUE(all(des[, 2] == trt))) {
          stop("the condition is fully confounded with 'strata'", call. = FALSE)
        }
      }
    }
  }
  wts <- 1 / (v + tau2)
  wts[!is.finite(wts) | !is.finite(b)] <- NA
  wts <- wts / rowMeans(wts, na.rm = TRUE)
  fit <- suppressWarnings(limma::lmFit(b, des, weights = wts))
  fit$Amean <- log(mean_expr + 1e-3)
  fit <- suppressWarnings(limma::eBayes(fit, robust = TRUE, trend = TRUE))
  j <- if (is.null(trt)) 1L else 2L
  est <- fit$coefficients[, j]
  tt <- fit$t[, j]
  df <- pmin(fit$df.total, .kishDf(wts, b, ncol(des)))
  if (is.null(trt)) {
    ok <- rowSums(is.finite(b)) >= min.pooled
  } else {
    ok <- rowSums(is.finite(b[, trt == 1, drop = FALSE])) >= min.group &
      rowSums(is.finite(b[, trt == 0, drop = FALSE])) >= min.group
  }
  tt[!ok | !is.finite(tt)] <- NA
  data.frame(estimate = unname(est), se = unname(est / tt), t = unname(tt), df = unname(df),
             p = unname(2 * stats::pt(-abs(tt), df)), n_patients = unname(rowSums(is.finite(b))),
             stringsAsFactors = FALSE)
}

# The slopes engine's tests for one index type's stored slopes.
.slopesTests <- function(x, trt = NULL, strata = NULL) {
  out <- list()
  vpool <- sweep(x$v_model, c(2, 3), x$factor, "*")
  tau_pool <- .dlTau2(x$beta, vpool)
  tau_cond <- if (!is.null(trt)) .dlTau2(x$beta, x$v_tile) else NULL
  for (j in seq_along(x$niches)) {
    b <- matrix(x$beta[, , j], nrow = length(x$genes))
    p <- .slopeColumnTest(b, matrix(vpool[, , j], nrow = length(x$genes)), tau_pool[, j], x$mean_expr)
    out[[length(out) + 1L]] <- data.frame(gene = x$genes, niche = x$niches[j], test = "pooled", p,
                                          stringsAsFactors = FALSE)
    if (!is.null(trt)) {
      cc <- .slopeColumnTest(b, matrix(x$v_tile[, , j], nrow = length(x$genes)), tau_cond[, j],
                             x$mean_expr, trt = trt, strata = strata)
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
