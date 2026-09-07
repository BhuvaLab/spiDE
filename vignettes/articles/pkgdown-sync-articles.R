#!/usr/bin/env Rscript
# Sync the validation-report articles from the research submodule. See README.
src <- "research/reports/benchmarks"
dst <- "vignettes/articles"
if (!file.exists(file.path(src, "build_site.R")))
  stop("research submodule not checked out (", src, " is empty) -- ",
       "in CI this means the deploy-key step failed")
rmd <- Sys.glob(file.path(src, "spiDE-*.Rmd"))
stopifnot(length(rmd) >= 6)
# a report retired from the submodule must not linger here as a stale copy:
# pkgdown lists every article it finds and refuses one missing from the index
stale <- setdiff(Sys.glob(file.path(dst, "spiDE-*.Rmd")), file.path(dst, basename(rmd)))
if (length(stale)) { unlink(stale); cat("removed stale:", paste(basename(stale), collapse = ", "), "\n") }
ok <- file.copy(rmd, dst, overwrite = TRUE, copy.date = TRUE)
stopifnot(all(ok))
# the cohort and rejected reports read tracked figures and summaries from the
# research tree; copy them to the sibling names those reports fall back to
sync_dir <- function(from, to) {
  unlink(to, recursive = TRUE); dir.create(to, showWarnings = FALSE)
  f <- list.files(from, full.names = TRUE, recursive = FALSE)
  f <- f[!dir.exists(f)]
  stopifnot(length(f) > 0, all(file.copy(f, to)))
  length(f)
}
n_fig <- sync_dir("research/fdr-ordering/figures", file.path(dst, "fdr-figures"))
n_sum <- sync_dir("research/fdr-ordering/summary", file.path(dst, "fdr-summary"))
n_str <- sync_dir("research/reports/data", file.path(dst, "stratum-data"))
n_nt  <- sync_dir("research/niche-transform/summary", file.path(dst, "niche-transform-summary"))
cat("synced supporting files:", n_fig, "figures,", n_sum + n_str + n_nt, "summaries\n")
# tables/ read at knit time relative to the article; copy for the same
# type-filter reason as the Rmds
unlink(file.path(dst, "tables"), recursive = TRUE)
dir.create(file.path(dst, "tables"), showWarnings = FALSE)
tb <- Sys.glob(file.path(src, "tables", "*.rds"))
stopifnot(length(tb) > 0, all(file.copy(tb, file.path(dst, "tables"))))
cat("synced", length(rmd), "articles and", length(tb), "tables from", src, "\n")
