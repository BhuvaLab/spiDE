# Convert the package's Rd files to Quarto pages under reference/r/ with altdoc's Rd converter
# (altdoc is installed in the git-ignored _rlib), plus an index grouped by topic.
# Examples are shown and run on every page except skip_pages (below).
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

# Every page's examples are RUN at build time, so the output appears under the code, except the
# pages in skip_pages (examples too slow to run on every build). The examples are
# self-contained (toy data, a second or two each: about half a minute for the whole
# reference), so no shared fit is needed. Figure pages get a wider device.
skip_pages <- character()   # none: nothing is too slow on the toy data
fig_pages <- c("plotNicheMap", "plotPatientSlopes", "plotNicheResponse", "plotCallMap",
               "plotTopTriplets", "plotSpillover", "plotTripletHeatmap", "plotPValues",
               "plotPatientEffects", "theme_spiDE", "spiDEColours")
missing_skip <- setdiff(c(skip_pages, fig_pages), sub("\\.Rd$", "", basename(rd_files)))
if (length(missing_skip)) stop("page lists name an absent Rd file: ", paste(missing_skip, collapse = ", "), call. = FALSE)

# front matter (page title = Rd name); examples executed only for run_pages
for (f in rd_files) {
  nm <- sub("\\.Rd$", "", basename(f))
  qmd <- file.path(out_dir, paste0(nm, ".qmd"))
  x <- readLines(qmd)
  run <- !(nm %in% skip_pages)
  x <- gsub("eval\\s*=\\s*(TRUE|FALSE)", if (run) "eval=TRUE" else "eval=FALSE", x)
  if (nm %in% fig_pages) x <- gsub("(```\\{r, warning=FALSE, message=FALSE, eval=TRUE)\\}", "\\1, fig.width=7, fig.height=4.5}", x)
  writeLines(c("---", sprintf("title: \"%s\"", nm), "toc: false",
               "execute:", sprintf("  eval: %s", tolower(run)), "---", "", x), qmd)
}

# explicit exclusion list: a page in skip_pages may not carry an executable chunk
for (nm in skip_pages) {
  f <- file.path(out_dir, paste0(nm, ".qmd"))
  if (any(grepl("eval\\s*[=:]\\s*(TRUE|true)", readLines(f)))) {
    stop("executable chunk left in ", f, " (listed in skip_pages)", call. = FALSE)
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
