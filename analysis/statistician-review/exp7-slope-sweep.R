# Focused check: does the RE-slope calibration of the ResponseNiche(A:B) null
# improve with sample size S? (Experiment 7 at S in {12, 24}, milder slope.)
suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(SpatialExperiment); library(SummarizedExperiment); library(S4Vectors)
})
SIGMA <- 30
simulate_slope <- function(seed, n_samples, sd_slope = 0.4, n_per = 90,
                           n_genes = 150, field = 500, sd_patient = 0.7) {
  set.seed(seed)
  cts <- c("A", "B", "C"); gn <- sprintf("G%d", seq_len(n_genes))
  sids <- sprintf("S%d", seq_len(n_samples))
  cond <- rep(c("Responder", "Non-responder"), length.out = n_samples)
  names(cond) <- sids
  cells <- lapply(sids, function(sid) {
    x <- runif(n_per, 0, field); y <- runif(n_per, 0, field)
    ct <- ifelse(x > 0.6 * field & runif(n_per) < 0.7, "B",
                 sample(c("A", "C"), n_per, replace = TRUE))
    data.frame(sample_id = sid, condition = cond[[sid]], x = x, y = y,
               cell_type = ct, stringsAsFactors = FALSE)
  })
  cd <- do.call(rbind, cells); n <- nrow(cd); cd$cell_id <- sprintf("c%d", seq_len(n))
  tmp <- SpatialExperiment(
    assays = list(counts = matrix(1L, n_genes, n, dimnames = list(gn, cd$cell_id))),
    colData = DataFrame(cd), spatialCoords = as.matrix(cd[, c("x", "y")]))
  tmp <- buildNiches(tmp, sigma = SIGMA)
  nb <- log1p(as.matrix(reducedDim(tmp, paste0("Niche", SIGMA)))[, "B"])
  nb <- (nb - mean(nb)) / stats::sd(nb)
  base <- matrix(rnorm(n_genes * 3, 1.5, 0.6), n_genes, dimnames = list(gn, cts))
  u <- matrix(rnorm(n_genes * n_samples, 0, sd_patient), n_genes, dimnames = list(gn, sids))
  v <- matrix(rnorm(n_genes * n_samples, 0, sd_slope), n_genes, dimnames = list(gn, sids))
  is_A <- cd$cell_type == "A"
  slope_mat <- v[, cd$sample_id, drop = FALSE]; slope_mat[, !is_A] <- 0
  lmu <- base[, cd$cell_type, drop = FALSE] + u[, cd$sample_id, drop = FALSE] +
    sweep(slope_mat, 2, nb, `*`)
  counts <- matrix(rnbinom(length(lmu), mu = as.vector(exp(lmu)), size = 5),
                   n_genes, dimnames = list(gn, cd$cell_id))
  spe <- SpatialExperiment(assays = list(counts = counts), colData = DataFrame(cd),
                           spatialCoords = as.matrix(cd[, c("x", "y")]))
  reducedDim(spe, paste0("Niche", SIGMA)) <- reducedDim(tmp, paste0("Niche", SIGMA))
  spe
}
rn_type1 <- function(spe, random, dfm) {
  Y <- as.matrix(assay(spe, "counts"))
  f <- fitSpiDE(spe, "condition", sigma = SIGMA, random = random,
                df.method = dfm, verbose = FALSE)
  ff <- spiDE:::.blockedInference(fits(f)[[1]], Y)
  ct <- as.character(ff@covtype); cm <- ff@coefmap
  sel <- ct == "ResponseNiche" & cm$index == "A" & cm$niche == "B"
  if (!any(sel)) return(NA_real_)
  cols <- cm$covariate[sel]; tt <- ff@t_stat[, cols, drop = FALSE]
  dfv <- ff@df
  dfc <- if (is.null(dfv)) Inf else if (length(dfv) == 1L) dfv else dfv[cols]
  dfm2 <- matrix(dfc, nrow(tt), ncol(tt), byrow = TRUE)
  mean(2 * stats::pt(-abs(tt), dfm2) < 0.05, na.rm = TRUE)
}
for (S in c(12L, 24L)) {
  for (s in 1:4) {
    spe <- simulate_slope(seed = 400 + s, n_samples = S, sd_slope = 0.4)
    fx <- rn_type1(spe, "none", "between")
    rb <- rn_type1(spe, "slope", "between")
    rs <- rn_type1(spe, "slope", "satterthwaite")
    cat(sprintf("S=%2d seed %d: fixed=%.3f  RE-slope/between=%.3f  RE-slope/satt=%.3f\n",
                S, s, fx, rb, rs))
  }
}
