# =============================================================================
# Quantify the shrinkage/pooling cost of a fixed-nested inference vs the
# random-intercept model, on the response (between-sample) contrast.
#
# Both approaches are calibrated; the ONE statistical difference is how the
# between-sample variance is estimated for each gene's response test:
#   fixed-nested : each gene's OWN between-sample variance   (S-2 df, noisy)
#   RE / pooled  : a between-sample variance SHARED across genes (spiDE's tau2),
#                  here the eBayes-moderated variance (limma::squeezeVar), the
#                  canonical realisation of "pool the variance across genes"
# We measure power (planted response effect) and type-I (null genes) by gene
# sparsity tier. The pooled arm should gain the most power on the SPARSE tier,
# where a single gene's own between-sample variance is least reliable -- that
# gap is the concrete cost of dropping shrinkage for the fixed-nested engine.
#
#   Rscript analysis/statistician-review/shrinkage.R [--quick]
# =============================================================================
QUICK  <- "--quick" %in% commandArgs(trailingOnly = TRUE)
S_grid <- if (QUICK) c(12L, 55L) else c(12L, 24L, 55L)  # samples (real ~55)
n_per  <- 60L                       # cells per sample
tau2   <- 0.5                       # between-sample (log-mean) variance, planted
eta    <- 0.8                       # planted responder log-fold-change (power)
reps   <- if (QUICK) 3L else 20L
G_tier <- if (QUICK) 120L else 400L # genes per sparsity tier (half planted)
# sparsity tiers = base log-mean; exp() ~ mean counts/cell in a baseline sample
tiers  <- c(sparse = 0.0, mid = 1.5, dense = 3.0)   # ~1, ~4.5, ~20 counts
have_limma <- requireNamespace("limma", quietly = TRUE)

rej <- function(p) mean(p < 0.05, na.rm = TRUE)

# per-gene response contrast + own between-sample variance (S-2 df), on log1p
gene_stats <- function(base, seed, S, resp, sample_of, col_idx) {
  set.seed(seed)
  G <- G_tier
  planted <- rep(c(TRUE, FALSE), each = G / 2)
  u <- matrix(rnorm(G * S, 0, sqrt(tau2)), G)
  lm_mat <- base + u + outer(as.numeric(planted) * eta, as.numeric(resp))
  Mu <- exp(lm_mat)[, sample_of]                          # G x (S*n_per)
  Y <- matrix(rnbinom(length(Mu), mu = as.vector(Mu), size = 5), nrow = G)
  ylog <- log1p(Y)
  sm <- vapply(col_idx, function(cc) rowMeans(ylog[, cc, drop = FALSE]),
               numeric(G))                                # G x S sample means
  eta_hat <- rowMeans(sm[, resp, drop = FALSE]) -
             rowMeans(sm[, !resp, drop = FALSE])
  ss <- function(M) rowSums((M - rowMeans(M))^2)
  s2 <- (ss(sm[, resp, drop = FALSE]) + ss(sm[, !resp, drop = FALSE])) / (S - 2)
  list(planted = planted, eta_hat = eta_hat, s2 = s2)
}

one_rep <- function(seed, S, resp, sample_of, col_idx, se_factor) {
  do.call(rbind, lapply(names(tiers), function(tn) {
    gs <- gene_stats(tiers[[tn]], seed + match(tn, names(tiers)) * 10000L,
                     S, resp, sample_of, col_idx)
    df_fn <- S - 2
    # fixed-nested: each gene's own between-sample variance, S-2 df
    p_fn <- 2 * stats::pt(-abs(gs$eta_hat / (sqrt(gs$s2) * se_factor)), df_fn)
    # RE / pooled: eBayes-moderated variance shared across genes (per tier, where
    # the mean-variance level is homogeneous), with the moderated df
    if (have_limma) {
      sq <- limma::squeezeVar(gs$s2, df = df_fn)
      s2_mod <- sq$var.post
      df_re <- df_fn + sq$df.prior
    } else {
      s2_mod <- mean(gs$s2)                 # fully-pooled fallback
      df_re <- df_fn * length(gs$s2)
    }
    p_re <- 2 * stats::pt(-abs(gs$eta_hat / (sqrt(s2_mod) * se_factor)), df_re)
    data.frame(
      tier = tn,
      power_fixednested = rej(p_fn[gs$planted]),
      power_RE          = rej(p_re[gs$planted]),
      type1_fixednested = rej(p_fn[!gs$planted]),
      type1_RE          = rej(p_re[!gs$planted]))
  }))
}

cat("Shrinkage cost: pooled (RE) vs per-gene (fixed-nested) between-sample variance\n")
cat(sprintf("  n_per=%d, tau2=%.2f, eta=%.2f, %d genes/tier, %d reps, S sweep=%s, limma=%s\n\n",
            n_per, tau2, eta, G_tier, reps, paste(S_grid, collapse = "/"),
            have_limma))
out <- do.call(rbind, lapply(S_grid, function(S) {
  cond <- rep(c("R", "N"), length.out = S)
  resp <- cond == "R"
  sample_of <- rep(seq_len(S), each = n_per)
  col_idx <- split(seq_len(S * n_per), sample_of)
  se_factor <- sqrt(1 / sum(resp) + 1 / sum(!resp))
  acc <- do.call(rbind, lapply(seq_len(reps), function(r)
    one_rep(r, S, resp, sample_of, col_idx, se_factor)))
  agg <- aggregate(cbind(power_fixednested, power_RE,
                         type1_fixednested, type1_RE) ~ tier, acc, mean)
  agg <- agg[match(names(tiers), agg$tier), ]
  agg$power_gap <- agg$power_RE - agg$power_fixednested
  cbind(S = S, agg)
}))
print(out, row.names = FALSE, digits = 3)
saveRDS(out, "analysis/statistician-review/shrinkage_results.rds")
cat("\nwrote analysis/statistician-review/shrinkage_results.rds\n")
