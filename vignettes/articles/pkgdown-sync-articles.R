#!/usr/bin/env Rscript
# Sync the validation-report articles from the research submodule. See README.
src <- "research/reports/benchmarks"
dst <- "vignettes/articles"
if (!file.exists(file.path(src, "build_site.R")))
  stop("research submodule not checked out (", src, " is empty) -- ",
       "in CI this means the deploy-key step failed")
rmd <- Sys.glob(file.path(src, "spiDE-*.Rmd"))
stopifnot(length(rmd) >= 6)
ok <- file.copy(rmd, dst, overwrite = TRUE, copy.date = TRUE)
stopifnot(all(ok))
# tables/ read at knit time relative to the article; copy for the same
# type-filter reason as the Rmds
unlink(file.path(dst, "tables"), recursive = TRUE)
dir.create(file.path(dst, "tables"), showWarnings = FALSE)
tb <- Sys.glob(file.path(src, "tables", "*.rds"))
stopifnot(length(tb) > 0, all(file.copy(tb, file.path(dst, "tables"))))
cat("synced", length(rmd), "articles and", length(tb), "tables from", src, "\n")
