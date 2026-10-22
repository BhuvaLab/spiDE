# Convert the package's Rd files to Quarto pages under reference/r/ with altdoc's Rd converter
# (altdoc is installed in the git-ignored _rlib), plus an index grouped by topic.
# Examples are shown, and run only for the plotting pages (run_pages below).
.libPaths(c(normalizePath("_rlib"), .libPaths()))
suppressPackageStartupMessages(library(altdoc))

man_dir <- "../man"
out_dir <- "reference/r"
unlink(out_dir, recursive = TRUE)
dir.create(out_dir, recursive = TRUE)

# altdoc's converter looks for its settings next to the package: give it a stub.
stub <- tempfile("altdoc-stub")
dir.create(file.path(stub, "altdoc"), recursive = TRUE)
file.create(file.path(stub, "altdoc", "quarto_website.yml"))
file.copy("../DESCRIPTION", file.path(stub, "DESCRIPTION"))

rd_files <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE)

# altdoc:::.rd2qmd is private API (altdoc is pinned in build_reference.sh). If it is ever
# missing, fall back to base R's Rd2HTML, embedded as raw HTML in a qmd page.
have_rd2qmd <- exists(".rd2qmd", envir = asNamespace("altdoc"), inherits = FALSE)
if (!have_rd2qmd) message("altdoc:::.rd2qmd not found; using tools::Rd2HTML fallback")
for (f in rd_files) {
  if (have_rd2qmd) {
    altdoc:::.rd2qmd(f, target_dir = out_dir, path = stub)
  } else {
    html <- tempfile(fileext = ".html")
    tools::Rd2HTML(f, out = html)
    h <- readLines(html)
    h <- h[seq(grep("<body", h)[1] + 1, tail(grep("</body>", h), 1) - 1)]
    writeLines(c("```{=html}", h, "```"),
               file.path(out_dir, sub("\\.Rd$", ".qmd", basename(f))))
  }
}

# Pages whose examples are RUN at build time, so the figures appear under the code. Every
# other page shows its examples unexecuted (spiDE(), fitSpiDE() and the rest are slow). The
# examples are self-contained (toy data), so no setup chunk is needed.
run_pages <- c("plotNicheMap", "plotPatientSlopes", "plotNicheResponse", "plotCallMap",
               "plotTopTriplets", "plotSpillover", "plotTripletHeatmap", "plotPValues",
               "plotPatientEffects", "theme_spiDE", "spiDEColours")
missing_run <- setdiff(run_pages, sub("\\.Rd$", "", basename(rd_files)))
if (length(missing_run)) stop("run_pages without an Rd file: ", paste(missing_run, collapse = ", "), call. = FALSE)

# front matter (page title = Rd name); examples executed only for run_pages
for (f in rd_files) {
  nm <- sub("\\.Rd$", "", basename(f))
  qmd <- file.path(out_dir, paste0(nm, ".qmd"))
  x <- readLines(qmd)
  run <- nm %in% run_pages
  x <- gsub("eval\\s*=\\s*(TRUE|FALSE)", if (run) "eval=TRUE" else "eval=FALSE", x)
  if (run) x <- gsub("(```\\{r, warning=FALSE, message=FALSE, eval=TRUE)\\}", "\\1, fig.width=7, fig.height=4.5}", x)
  writeLines(c("---", sprintf("title: \"%s\"", nm), "toc: false",
               "execute:", sprintf("  eval: %s", tolower(run)), "---", "", x), qmd)
}

# explicit allow-list: only run_pages may carry an executable chunk, so a slow example
# is never run by accident
for (f in list.files(out_dir, pattern = "\\.qmd$", full.names = TRUE)) {
  nm <- sub("\\.qmd$", "", basename(f))
  if (!(nm %in% run_pages) && any(grepl("eval\\s*[=:]\\s*(TRUE|true)", readLines(f)))) {
    stop("executable chunk left in ", f, " (not in run_pages)", call. = FALSE)
  }
}

# index grouped by topic, as in _pkgdown.yml; Rd title taken from the Rd files.
# Pages not listed in a group land in "Other", so a new Rd file is never dropped.
groups <- list(
  "Analysis" = c("spiDE", "buildNiches", "mergeNiches", "fitSpiDE", "testSpiDE", "results",
                 "spiGSEA", "testNicheAbundance"),
  "Per-patient estimates" = c("patientSlopes", "patientIntercepts"),
  "Plots" = c("plotNicheMap", "plotPatientSlopes", "plotNicheResponse", "plotCallMap",
              "plotTopTriplets", "plotSpillover", "plotTripletHeatmap", "plotPValues",
              "plotPatientEffects", "theme_spiDE", "spiDEColours"),
  "Classes and data" = c("SpiDEFit-class", "SpiDEResults-class", "toySpiDE", "spiDE-package")
)
all_pages <- sub("\\.Rd$", "", basename(rd_files))
other <- setdiff(all_pages, unlist(groups))
if (length(other)) groups[["Other"]] <- other
groups <- lapply(groups, function(g) intersect(g, all_pages))

rd_title <- function(nm) {
  rd <- tools::parse_Rd(file.path(man_dir, paste0(nm, ".Rd")))
  tag <- vapply(rd, function(e) attr(e, "Rd_tag"), "")
  trimws(gsub("\\s+", " ", paste(unlist(rd[tag == "\\title"]), collapse = "")))
}

idx <- c("---", "title: \"R reference\"", "toc: false", "---", "")
for (g in names(groups)) {
  idx <- c(idx, paste("##", g), "", "| | |", "| --- | --- |")
  for (nm in groups[[g]]) {
    idx <- c(idx, sprintf("| [`%s`](%s.qmd) | %s |", nm, nm, rd_title(nm)))
  }
  idx <- c(idx, "")
}
writeLines(idx, file.path(out_dir, "index.qmd"))
message("wrote ", length(rd_files), " R reference pages")
