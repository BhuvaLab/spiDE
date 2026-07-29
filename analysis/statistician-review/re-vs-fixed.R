# =============================================================================
# Responding to the statistician's review of the spiDE mixed-effects model.
#
# Investigation (not a package feature): does the page-5 "pseudo-replication"
# bias in the spiDE-model vignette survive the experiments the reviewer asked
# for, and is the random-effects (RE) model actually needed vs a properly
# *nested* fixed-effects model? Uses the vignette's own `sim_clustered` design,
# varying only HOW the per-sample effects are generated.
#
#   Run (background; slow, tens of minutes):
#     Rscript analysis/statistician-review/re-vs-fixed.R          # full
#     Rscript analysis/statistician-review/re-vs-fixed.R --quick  # fast smoke
#
# Writes analysis/statistician-review/results.rds and prints tables.
# =============================================================================

suppressPackageStartupMessages({
  # load THIS branch's spiDE (df.method); library(spiDE) would use a stale
  # installed build that predates the argument and error.
  devtools::load_all(quiet = TRUE)
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
})

QUICK <- "--quick" %in% commandArgs(trailingOnly = TRUE)
cfg <- if (QUICK) {
  list(n_samples = 8L, n_per = 60L, n_genes = 60L, seeds = 1:2,
       demo_genes = 30L, mle_genes = 6L)
} else {
  list(n_samples = 12L, n_per = 90L, n_genes = 200L, seeds = 1:5,
       demo_genes = 150L, mle_genes = 20L)
}
SIGMA <- 30
SD_PATIENT <- 0.7
set.seed(1)
cat(sprintf("mode: %s | n_samples=%d n_per=%d n_genes=%d seeds=%s\n",
            if (QUICK) "QUICK" else "FULL", cfg$n_samples, cfg$n_per,
            cfg$n_genes, paste(range(cfg$seeds), collapse = "-")))

# ---------------------------------------------------------------------------
# Simulator: the vignette's sim_clustered, with a `u_mode` controlling how the
# per-(gene, sample) intercepts u[g,s] relate to the responder/non-responder
# grouping. `iid` is exactly the vignette (u ~ N(0, sd^2), condition-independent).
#   centered : subtract the per-group mean per gene  -> zero mean within group
#   equal    : centered, then + m to BOTH groups     -> equal non-zero group means
#   unequal  : centered, then + delta to responders  -> a real group difference
# ---------------------------------------------------------------------------
simulate_clustered <- function(seed, u_mode = "iid", delta = 0.6, m = 0.5,
                               n_samples = cfg$n_samples, n_per = cfg$n_per,
                               n_genes = cfg$n_genes, field = 500,
                               sd_patient = SD_PATIENT) {
  set.seed(seed)
  cts <- c("A", "B", "C")
  gn <- sprintf("G%d", seq_len(n_genes))
  sids <- sprintf("S%d", seq_len(n_samples))
  cond <- rep(c("Responder", "Non-responder"), length.out = n_samples)
  names(cond) <- sids
  cells <- lapply(sids, function(sid) {
    x <- runif(n_per, 0, field)
    y <- runif(n_per, 0, field)
    ct <- ifelse(x > 0.6 * field & runif(n_per) < 0.7, "B",
                 sample(c("A", "C"), n_per, replace = TRUE))
    data.frame(sample_id = sid, condition = cond[[sid]], x = x, y = y,
               cell_type = ct, stringsAsFactors = FALSE)
  })
  cd <- do.call(rbind, cells)
  n <- nrow(cd)
  cd$cell_id <- sprintf("c%d", seq_len(n))
  base <- matrix(rnorm(n_genes * 3, 1.5, 0.6), n_genes,
                 dimnames = list(gn, cts))
  u <- matrix(rnorm(n_genes * n_samples, 0, sd_patient), n_genes,
              dimnames = list(gn, sids))
  resp_s <- which(cond == "Responder")
  nonr_s <- which(cond == "Non-responder")
  if (u_mode != "iid") {
    u[, resp_s] <- u[, resp_s] - rowMeans(u[, resp_s, drop = FALSE])
    u[, nonr_s] <- u[, nonr_s] - rowMeans(u[, nonr_s, drop = FALSE])
    if (u_mode == "equal") {
      u[, resp_s] <- u[, resp_s] + m
      u[, nonr_s] <- u[, nonr_s] + m
    } else if (u_mode == "unequal") {
      u[, resp_s] <- u[, resp_s] + delta
    }
  }
  lmu <- base[, cd$cell_type, drop = FALSE] + u[, cd$sample_id, drop = FALSE]
  counts <- matrix(rnbinom(length(lmu), mu = as.vector(exp(lmu)), size = 5),
                   n_genes, dimnames = list(gn, cd$cell_id))
  spe <- SpatialExperiment(assays = list(counts = counts),
                           colData = DataFrame(cd),
                           spatialCoords = as.matrix(cd[, c("x", "y")]))
  # per-gene realised group-mean difference in u (what 'resp' can pick up)
  delta_g <- rowMeans(u[, resp_s, drop = FALSE]) -
             rowMeans(u[, nonr_s, drop = FALSE])
  list(spe = spe, u = u, cond = cond, delta_g = delta_g)
}

# Response-effect null type-I for a fitted spiDE mode. Returns per-gene p-values.
resp_pvals <- function(spe, Y, random) {
  f <- spiDE::fitSpiDE(spe, "condition", sigma = SIGMA, random = random,
                       df.method = "between", verbose = FALSE)
  ff <- spiDE:::.blockedInference(fits(f)[[1]], Y)
  ct <- as.character(ff@covtype)
  rc <- ff@coefmap$covariate[ct == "Response"]
  dfv <- ff@df
  df_resp <- if (is.null(dfv)) Inf else if (length(dfv) == 1L) dfv else dfv[[rc]]
  list(p = 2 * stats::pt(-abs(ff@t_stat[, rc]), df_resp),
       tau2 = if (is.null(ff@tau2)) NA_real_ else ff@tau2[["SampleInt"]])
}

# ===========================================================================
# Experiments 1-5: intercept u_mode sweep (fixed vs mixed Response type-I)
# ===========================================================================
cat("\n== Experiments 1-5: per-sample intercept generation ==\n")
modes <- c("iid", "centered", "equal", "unequal")
rows <- list()
delta_by_gene <- list()   # for the decomposition (experiment 5)
for (um in modes) {
  for (s in cfg$seeds) {
    sim <- simulate_clustered(seed = s, u_mode = um)
    spe <- buildNiches(sim$spe, sigma = SIGMA)
    Y <- as.matrix(assay(spe, "counts"))
    fx <- resp_pvals(spe, Y, "none")
    mx <- resp_pvals(spe, Y, "intercept")
    rows[[length(rows) + 1L]] <- data.frame(
      u_mode = um, seed = s,
      fixed_type1 = mean(fx$p < 0.05, na.rm = TRUE),
      mixed_type1 = mean(mx$p < 0.05, na.rm = TRUE),
      mean_abs_delta = mean(abs(sim$delta_g)),
      tau2 = mx$tau2)
    if (um == "iid") {
      delta_by_gene[[length(delta_by_gene) + 1L]] <- data.frame(
        seed = s, delta_g = sim$delta_g,
        z_fixed = stats::qnorm(pmin(pmax(fx$p / 2, 1e-12), 1 - 1e-12),
                               lower.tail = FALSE))
    }
    cat(sprintf("  %-9s seed %d: fixed=%.3f mixed=%.3f mean|delta|=%.3f\n",
                um, s, mean(fx$p < 0.05), mean(mx$p < 0.05),
                mean(abs(sim$delta_g))))
  }
}
exp15 <- do.call(rbind, rows)
exp5 <- do.call(rbind, delta_by_gene)
cat("\n-- Experiments 1-4 summary (mean over seeds) --\n")
agg15 <- aggregate(cbind(fixed_type1, mixed_type1, mean_abs_delta) ~ u_mode,
                   exp15, mean)
agg15 <- agg15[match(modes, agg15$u_mode), ]
print(agg15, row.names = FALSE)
# experiment 5: does the fixed test track the per-gene group-mean difference?
exp5$abs_z <- abs(exp5$z_fixed)
cor_zd <- suppressWarnings(cor(abs(exp5$delta_g), exp5$abs_z,
                               use = "complete.obs"))
cat(sprintf("\n-- Experiment 5: cor(|delta_g|, |fixed z|) under iid = %.3f --\n",
            cor_zd))

# ===========================================================================
# Experiment 6: his nested (contr.sum) fixed model vs the error stratum.
# Per-gene demonstration on one iid dataset, Gaussian on log1p(counts):
#   (a) nested-fixed Wald with CELL residual   -> expect anti-conservative
#   (b) between-sample test (sample-mean t / aov Error stratum, S-2 df)
#   (c) RE via lme4::lmer                       -> expect calibrated, ~ (b)
# ===========================================================================
cat("\n== Experiment 6: nested fixed (contr.sum) vs error stratum ==\n")
make_nested <- function(sample_id, condition) {
  Z <- NULL
  for (g in unique(condition)) {
    samps <- unique(sample_id[condition == g])
    Cm <- stats::contr.sum(length(samps))
    rownames(Cm) <- samps
    Zg <- Cm[match(sample_id, rownames(Cm)), , drop = FALSE]
    Zg[is.na(Zg)] <- 0
    colnames(Zg) <- paste0("g", which(unique(condition) == g), "_c",
                           seq_len(ncol(Zg)))
    Z <- if (is.null(Z)) Zg else cbind(Z, Zg)
  }
  Z
}
have_lme4 <- requireNamespace("lme4", quietly = TRUE)
have_lmerTest <- requireNamespace("lmerTest", quietly = TRUE)
sim6 <- simulate_clustered(seed = 101, u_mode = "iid")
cd6 <- as.data.frame(colData(sim6$spe))
Y6 <- as.matrix(assay(sim6$spe, "counts"))
Znest <- make_nested(cd6$sample_id, cd6$condition)
condf <- factor(cd6$condition)
p6 <- t(vapply(seq_len(min(cfg$demo_genes, nrow(Y6))), function(gi) {
  y <- log1p(Y6[gi, ])
  # (a) nested-fixed, condition tested at the CELL residual
  m_nest <- stats::lm(y ~ condf + Znest)
  sm <- summary(m_nest)$coefficients
  crow <- grep("^condf", rownames(sm))[1]
  p_nested_cell <- if (is.na(crow)) NA_real_ else sm[crow, 4]
  # (b) between-sample: t-test on the S sample means (S-2 df)
  sm_means <- tapply(y, cd6$sample_id, mean)
  grp <- tapply(as.character(cd6$condition), cd6$sample_id,
                function(z) z[1])[names(sm_means)]
  p_between <- tryCatch(stats::t.test(sm_means ~ grp)$p.value,
                        error = function(e) NA_real_)
  # (c) RE
  p_re <- NA_real_
  if (have_lmerTest) {
    fit <- try(suppressMessages(lmerTest::lmer(
      y ~ condf + (1 | cd6$sample_id))), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      cc <- summary(fit)$coefficients
      rr <- grep("^condf", rownames(cc))[1]
      if (!is.na(rr)) p_re <- cc[rr, "Pr(>|t|)"]
    }
  }
  c(nested_cell = p_nested_cell, between = p_between, re = p_re)
}, numeric(3)))
exp6 <- data.frame(
  approach = c("nested-fixed (cell SE)", "between-sample (S-2 df)",
               "RE lmer (Satterthwaite)"),
  type1 = c(mean(p6[, "nested_cell"] < 0.05, na.rm = TRUE),
            mean(p6[, "between"] < 0.05, na.rm = TRUE),
            mean(p6[, "re"] < 0.05, na.rm = TRUE)))
cat("-- Experiment 6: Response-effect null type-I by approach --\n")
print(exp6, row.names = FALSE)

# ===========================================================================
# Experiment 7: the SLOPE case (the real spiDE target). Plant between-sample
# variation in the celltype-A x niche-B slope (true ResponseNiche beta = 0);
# compare random="none" vs random="slope" ResponseNiche (A x B) null type-I.
# ===========================================================================
cat("\n== Experiment 7: between-sample niche-slope null (why RE for slopes) ==\n")
simulate_slope <- function(seed, sd_slope = 0.6, n_samples = cfg$n_samples,
                           n_per = cfg$n_per, n_genes = cfg$n_genes,
                           field = 500, sd_patient = SD_PATIENT) {
  set.seed(seed)
  cts <- c("A", "B", "C")
  gn <- sprintf("G%d", seq_len(n_genes))
  sids <- sprintf("S%d", seq_len(n_samples))
  cond <- rep(c("Responder", "Non-responder"), length.out = n_samples)
  names(cond) <- sids
  cells <- lapply(sids, function(sid) {
    x <- runif(n_per, 0, field)
    y <- runif(n_per, 0, field)
    ct <- ifelse(x > 0.6 * field & runif(n_per) < 0.7, "B",
                 sample(c("A", "C"), n_per, replace = TRUE))
    data.frame(sample_id = sid, condition = cond[[sid]], x = x, y = y,
               cell_type = ct, stringsAsFactors = FALSE)
  })
  cd <- do.call(rbind, cells)
  n <- nrow(cd)
  cd$cell_id <- sprintf("c%d", seq_len(n))
  # build niches on placeholder counts to get the B-niche covariate for planting
  tmp <- SpatialExperiment(
    assays = list(counts = matrix(1L, n_genes, n,
                                  dimnames = list(gn, cd$cell_id))),
    colData = DataFrame(cd), spatialCoords = as.matrix(cd[, c("x", "y")]))
  tmp <- buildNiches(tmp, sigma = SIGMA)
  nb <- log1p(as.matrix(reducedDim(tmp, paste0("Niche", SIGMA)))[, "B"])
  nb <- (nb - mean(nb)) / stats::sd(nb)                 # standardise
  base <- matrix(rnorm(n_genes * 3, 1.5, 0.6), n_genes,
                 dimnames = list(gn, cts))
  u <- matrix(rnorm(n_genes * n_samples, 0, sd_patient), n_genes,
              dimnames = list(gn, sids))                # intercept clustering
  v <- matrix(rnorm(n_genes * n_samples, 0, sd_slope), n_genes,
              dimnames = list(gn, sids))                # per-sample SLOPE, mean 0
  # A-cells get a per-sample niche-B slope (condition-independent -> beta = 0)
  is_A <- cd$cell_type == "A"
  slope_mat <- v[, cd$sample_id, drop = FALSE]          # genes x cells
  slope_mat[, !is_A] <- 0                                # only A cells slope on niche
  lmu <- base[, cd$cell_type, drop = FALSE] +
         u[, cd$sample_id, drop = FALSE] +
         sweep(slope_mat, 2, nb, `*`)
  counts <- matrix(rnbinom(length(lmu), mu = as.vector(exp(lmu)), size = 5),
                   n_genes, dimnames = list(gn, cd$cell_id))
  spe <- SpatialExperiment(assays = list(counts = counts),
                           colData = DataFrame(cd),
                           spatialCoords = as.matrix(cd[, c("x", "y")]))
  reducedDim(spe, paste0("Niche", SIGMA)) <-
    reducedDim(tmp, paste0("Niche", SIGMA))              # reuse the niche
  spe
}
rn_type1 <- function(spe, random) {
  Y <- as.matrix(assay(spe, "counts"))
  f <- spiDE::fitSpiDE(spe, "condition", sigma = SIGMA, random = random,
                       df.method = "between", verbose = FALSE)
  ff <- spiDE:::.blockedInference(fits(f)[[1]], Y)
  ct <- as.character(ff@covtype)
  cm <- ff@coefmap
  # the A x B ResponseNiche column(s)
  sel <- ct == "ResponseNiche" & cm$index == "A" & cm$niche == "B"
  if (!any(sel)) return(NA_real_)
  cols <- cm$covariate[sel]
  tt <- ff@t_stat[, cols, drop = FALSE]
  dfv <- ff@df
  dfc <- if (is.null(dfv)) Inf else if (length(dfv) == 1L) dfv else dfv[cols]
  dfm <- matrix(dfc, nrow(tt), ncol(tt), byrow = TRUE)
  mean(2 * stats::pt(-abs(tt), dfm) < 0.05, na.rm = TRUE)
}
rows7 <- list()
for (s in cfg$seeds) {
  spe <- simulate_slope(seed = 200 + s)
  rows7[[length(rows7) + 1L]] <- data.frame(
    seed = s,
    fixed = rn_type1(spe, "none"),
    slope = rn_type1(spe, "slope"))
  cat(sprintf("  slope-null seed %d: fixed(A:B)=%.3f  RE-slope(A:B)=%.3f\n",
              s, rows7[[length(rows7)]]$fixed, rows7[[length(rows7)]]$slope))
}
exp7 <- do.call(rbind, rows7)
cat(sprintf("-- Experiment 7 (mean): fixed=%.3f  RE-slope=%.3f --\n",
            mean(exp7$fixed, na.rm = TRUE), mean(exp7$slope, na.rm = TRUE)))

# ===========================================================================
# Experiment 8 (light): spiDE pooled Schall/MoM tau2 vs per-gene REML tau2.
# ===========================================================================
cat("\n== Experiment 8: MoM (pooled) vs REML (per-gene) tau2 ==\n")
sim8 <- simulate_clustered(seed = 301, u_mode = "iid")
spe8 <- buildNiches(sim8$spe, sigma = SIGMA)
t0 <- proc.time()[["elapsed"]]
f8 <- spiDE::fitSpiDE(spe8, "condition", sigma = SIGMA, random = "intercept",
                      df.method = "between", verbose = FALSE)
mom_time <- proc.time()[["elapsed"]] - t0
tau2_mom <- fits(f8)[[1]]@tau2[["SampleInt"]]
# per-gene REML/MLE on the MATCHED NB model (glmmPQL, theta = 5 = the sim's NB
# size) so the tau2 is on the same log-mean scale as spiDE's Schall/MoM estimate.
reml_tau2 <- rep(NA_real_, cfg$mle_genes)
Y8 <- as.matrix(assay(spe8, "counts"))
cd8 <- as.data.frame(colData(spe8))
condf8 <- factor(cd8$condition)
sidf8 <- factor(cd8$sample_id)
have_pql <- requireNamespace("MASS", quietly = TRUE) &&
  requireNamespace("nlme", quietly = TRUE)
t0 <- proc.time()[["elapsed"]]
if (have_pql) {
  for (gi in seq_len(cfg$mle_genes)) {
    dat <- data.frame(count = Y8[gi, ], cond = condf8, sid = sidf8)
    fit <- try(suppressWarnings(suppressMessages(MASS::glmmPQL(
      count ~ cond, random = ~ 1 | sid,
      family = MASS::negative.binomial(theta = 5), data = dat,
      verbose = FALSE))), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      sd_int <- as.numeric(nlme::VarCorr(fit)["(Intercept)", "StdDev"])
      reml_tau2[gi] <- sd_int^2
    }
  }
}
reml_time <- proc.time()[["elapsed"]] - t0
cat(sprintf("  pooled MoM tau2 = %.3f (%.1fs, all %d genes)\n",
            tau2_mom, mom_time, cfg$n_genes))
cat(sprintf("  REML tau2 per-gene: mean=%.3f median=%.3f sd=%.3f (%.1fs, %d genes)\n",
            mean(reml_tau2, na.rm = TRUE), stats::median(reml_tau2, na.rm = TRUE),
            stats::sd(reml_tau2, na.rm = TRUE), reml_time, cfg$mle_genes))
cat(sprintf("  true simulated sd_patient^2 = %.3f\n", SD_PATIENT^2))
exp8 <- list(tau2_mom = tau2_mom, mom_time = mom_time,
             reml_tau2 = reml_tau2, reml_time = reml_time,
             truth = SD_PATIENT^2)

# ---------------------------------------------------------------------------
out <- list(config = cfg, exp15 = exp15, agg15 = agg15, exp5 = exp5,
            cor_zd = cor_zd, exp6 = exp6, exp7 = exp7, exp8 = exp8)
saveRDS(out, "analysis/statistician-review/results.rds")
cat("\nwrote analysis/statistician-review/results.rds\n")
