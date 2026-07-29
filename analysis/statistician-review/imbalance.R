# =============================================================================
# Task 2: the reviewer's nested contr.sum fixed model vs the random-intercept
# model, under BALANCED and 3-fold cell-count IMBALANCE.
#
# Main-effect (resp) null: per-sample intercepts u ~ N(0, sd^2) i.i.d., no
# response effect. Arms (transparent Gaussian on log1p counts):
#   fixed-cell      : lm(y ~ cond), Wald at the cell residual df
#   nested-cell     : lm(y ~ cond + contr.sum(individual within group)), Wald,
#                     cell residual SE  -> the reviewer's nesting, tested wrongly
#   nested-between  : the same nesting, cond tested against the between-individual
#                     stratum = unweighted t-test on the S sample means (the
#                     balanced-ANOVA error term; assumes equal cell counts)
#   RE (lmerTest)   : lmer(y ~ cond + (1|sample)), REML + Satterthwaite (weights
#                     unequal cells correctly)
#
# Balanced: nested-between and RE should both be ~0.05 (equivalent). Under
# 3-fold imbalance the unweighted between-stratum test can drift (its equal-
# variance assumption breaks), whereas REML weights the samples by information
# and should stay calibrated -- the practical case for the mixed model over a
# hand-rolled nested-fixed ANOVA.
#
#   Rscript analysis/statistician-review/imbalance.R [--quick]
# =============================================================================
suppressPackageStartupMessages({ library(SpatialExperiment); library(S4Vectors) })
QUICK <- "--quick" %in% commandArgs(trailingOnly = TRUE)
n_samples <- 12L
n_genes   <- if (QUICK) 60L else 200L
seeds     <- if (QUICK) 1:2 else 1:6
n_per     <- 90L
SD_PAT    <- 0.7
have_lmerTest <- requireNamespace("lmerTest", quietly = TRUE)

# per-sample cell counts: 3-fold imbalance, balanced ACROSS conditions so it is
# not confounded with the tested contrast (half of each group large, half small)
cell_counts <- function(cond, imbalance) {
  n <- rep(n_per, length(cond))
  if (imbalance) for (g in unique(cond)) {
    gs <- which(cond == g)
    n[gs[seq_len(floor(length(gs) / 2))]] <- 3L * n_per
  }
  n
}

simulate_null <- function(seed, imbalance) {
  set.seed(seed)
  cts <- c("A", "B", "C"); gn <- sprintf("G%d", seq_len(n_genes))
  sids <- sprintf("S%d", seq_len(n_samples))
  cond <- rep(c("Responder", "Non-responder"), length.out = n_samples)
  nps <- cell_counts(cond, imbalance)
  cd <- do.call(rbind, lapply(seq_along(sids), function(i) {
    ct <- sample(cts, nps[i], replace = TRUE)
    data.frame(sample_id = sids[i], condition = cond[i], cell_type = ct,
               stringsAsFactors = FALSE)
  }))
  n <- nrow(cd)
  base <- matrix(rnorm(n_genes * 3, 1.5, 0.6), n_genes, dimnames = list(gn, cts))
  u <- matrix(rnorm(n_genes * n_samples, 0, SD_PAT), n_genes,
              dimnames = list(gn, sids))
  lmu <- base[, cd$cell_type, drop = FALSE] + u[, cd$sample_id, drop = FALSE]
  counts <- matrix(rnbinom(length(lmu), mu = as.vector(exp(lmu)), size = 5),
                   n_genes, dimnames = list(gn, sprintf("c%d", seq_len(n))))
  list(counts = counts, cd = cd)
}

# contr.sum coding of individuals nested within condition group (S-2 columns)
make_nested <- function(sample_id, condition) {
  Z <- NULL
  for (g in unique(condition)) {
    samps <- unique(sample_id[condition == g])
    Cm <- stats::contr.sum(length(samps)); rownames(Cm) <- samps
    Zg <- Cm[match(sample_id, rownames(Cm)), , drop = FALSE]; Zg[is.na(Zg)] <- 0
    colnames(Zg) <- paste0(make.names(g), "_c", seq_len(ncol(Zg)))
    Z <- if (is.null(Z)) Zg else cbind(Z, Zg)
  }
  Z
}

arms_type1 <- function(sim) {
  cd <- sim$cd; Y <- sim$counts
  condf <- factor(cd$condition); sidf <- factor(cd$sample_id)
  Znest <- make_nested(cd$sample_id, cd$condition)
  grp <- tapply(as.character(cd$condition), cd$sample_id, `[`, 1)
  p <- t(vapply(seq_len(nrow(Y)), function(gi) {
    y <- log1p(Y[gi, ])
    p_cell <- summary(stats::lm(y ~ condf))$coefficients[2, 4]
    sm <- summary(stats::lm(y ~ condf + Znest))$coefficients
    cr <- grep("^condf", rownames(sm))[1]
    p_ncell <- if (is.na(cr)) NA_real_ else sm[cr, 4]
    means <- tapply(y, cd$sample_id, mean)
    p_btw <- tryCatch(stats::t.test(means ~ grp[names(means)])$p.value,
                      error = function(e) NA_real_)
    p_re <- NA_real_
    if (have_lmerTest) {
      f <- try(suppressMessages(lmerTest::lmer(y ~ condf + (1 | sidf))),
               silent = TRUE)
      if (!inherits(f, "try-error")) {
        cc <- summary(f)$coefficients; rr <- grep("^condf", rownames(cc))[1]
        if (!is.na(rr)) p_re <- cc[rr, "Pr(>|t|)"]
      }
    }
    c(`fixed-cell` = p_cell, `nested-cell` = p_ncell,
      `nested-between` = p_btw, `RE-lmer` = p_re)
  }, numeric(4)))
  colMeans(p < 0.05, na.rm = TRUE)
}

cat(sprintf("Task 2: nested contr.sum fixed vs RE, balanced vs 3-fold imbalance\n"))
cat(sprintf("  S=%d, n_genes=%d, seeds=%s, n_per=%d (imbalanced: %d vs %d)\n\n",
            n_samples, n_genes, paste(range(seeds), collapse = "-"), n_per,
            n_per, 3L * n_per))
res <- list()
for (imb in c(FALSE, TRUE)) {
  acc <- NULL
  for (s in seeds) acc <- rbind(acc, arms_type1(simulate_null(s, imb)))
  m <- colMeans(acc)
  res[[if (imb) "3-fold imbalance" else "balanced"]] <- m
  cat(sprintf("%-16s  fixed-cell=%.3f  nested-cell=%.3f  nested-between=%.3f  RE-lmer=%.3f\n",
              if (imb) "3-fold imbal:" else "balanced:",
              m["fixed-cell"], m["nested-cell"], m["nested-between"], m["RE-lmer"]))
}
saveRDS(res, "analysis/statistician-review/imbalance_results.rds")
cat("\nwrote analysis/statistician-review/imbalance_results.rds\n")
