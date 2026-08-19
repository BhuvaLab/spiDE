#!/usr/bin/env Rscript
# Calibration + validity diagnostic for a spiDE results table.
#
#   Rscript calibration_check.R <results.rds> [fit.rds]
#
# Prints, in one pass:
#   1. global calibration      sd(t), frac|t|>1.96, raw p<0.05
#   2. per-index calibration   ranked, with cells/subset if the fit is given
#   3. dropout ~ condition     Fisher test per index (needs the fit object)
#   4. the valid subset        indexes passing BOTH criteria, re-FDR'd alone
#
# Reference values measured on the YTMA cohort, 2026-08-13 (sigma = 30):
#   calibrated   Tumor      sd(t) 1.05, frac 0.061 at 388 cells/subset
#   inflated     Mast       sd(t) 1.65, frac 0.226 at  38 cells/subset
#   driver       corr(sd_t, median cells) = -0.838; corr(sd_t, patients) = -0.659
#   NOT drivers  bandwidth (sd 1.363 at bw30 vs 1.364 at bw70), patient count
suppressPackageStartupMessages({library(stats)})
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: calibration_check.R <results.rds> [fit.rds]")
tb <- readRDS(args[1])
stopifnot(all(c("t", "p.niche", "ct_index") %in% names(tb)))

cat("== 1. global ==\n")
cat(sprintf("   triplets %d | sd(t) %.3f | frac|t|>1.96 %.4f | raw p<.05 %.4f\n",
            nrow(tb), sd(tb$t, na.rm = TRUE),
            mean(abs(tb$t) > 1.96, na.rm = TRUE),
            mean(tb$p.niche < 0.05, na.rm = TRUE)))
cat("   (sd(t) ~ 1.0 and frac ~ 0.05 if calibrated; >1.2 means FDR is not valid)\n")

cat("\n== 2. per index ==\n")
d <- do.call(rbind, lapply(split(tb, tb$ct_index), function(z) data.frame(
  index = z$ct_index[1], n = nrow(z), sd_t = sd(z$t, na.rm = TRUE),
  frac_sig = mean(abs(z$t) > 1.96, na.rm = TRUE))))

inc <- NULL
if (length(args) > 1 && file.exists(args[2])) {
  fit <- readRDS(args[2])
  inc <- tryCatch(fit@diagnostics$inclusion, error = function(e) NULL)
}
if (!is.null(inc)) {
  ds <- do.call(rbind, lapply(split(inc, inc$index), function(z) {
    tt <- table(z$condition, z$included)
    p <- if (all(dim(tt) == c(2, 2))) tryCatch(fisher.test(tt)$p.value,
                                               error = function(e) NA_real_) else NA_real_
    data.frame(index = z$index[1], included = sum(z$included), of = nrow(z),
               dropout_p = p)
  }))
  d <- merge(d, ds, by = "index", all.x = TRUE)
}
print(d[order(d$sd_t), ], row.names = FALSE, digits = 3)

if (!is.null(inc)) {
  cat("\n== 3. dropout ~ condition (informative missingness) ==\n")
  bad <- d$index[!is.na(d$dropout_p) & d$dropout_p < 0.05]
  if (length(bad)) {
    cat(sprintf("   CONFOUNDED: %s\n", paste(bad, collapse = ", ")))
    cat("   Which patients contribute depends on their condition. No threshold\n")
    cat("   or variance correction fixes this -- the missingness is informative.\n")
  } else cat("   none detected\n")

  cat("\n== 4. valid subset (complete inclusion AND sd(t) < 1.2) ==\n")
  ok <- d$index[d$included == d$of & d$sd_t < 1.2]
  ok <- ok[!is.na(ok)]
  if (!length(ok)) {
    cat("   NO index type passes both criteria.\n")
  } else {
    cat(sprintf("   %s\n", paste(ok, collapse = ", ")))
    z <- tb[tb$ct_index %in% ok, ]
    z$q <- p.adjust(z$p.niche, "BH")     # re-FDR within the valid subset ONLY
    cat(sprintf("   %d triplets | sd(t) %.3f | raw p<.05 %.4f\n", nrow(z),
                sd(z$t, na.rm = TRUE), mean(z$p.niche < 0.05, na.rm = TRUE)))
    for (thr in c(0.05, 0.10, 0.20))
      cat(sprintf("   q<=%.2f: %4d triplets / %4d genes\n", thr,
                  sum(z$q <= thr, na.rm = TRUE),
                  length(unique(z$gene[z$q <= thr & !is.na(z$q)]))))
    cat("\n   NOTE: this restriction is post-hoc. It is defensible only because\n")
    cat("   both criteria (complete inclusion, calibrated variance) are\n")
    cat("   independent of the outcome -- state it that way in any write-up.\n")
  }
}
