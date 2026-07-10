# Coarsen fine cell-type niche columns into named groups by summing their
# density columns. Cell types not named in any group are carried through
# unchanged, matching the niche_merge logic in the original analysis scripts.

#' Merge niche cell-type columns into groups
#'
#' Sums the density columns of grouped cell types in each niche
#' \code{reducedDim}, producing coarser niche covariates (e.g. merging
#' Macrophage subtypes into a single "Macrophage" niche). Cell types absent from
#' \code{groups} are retained as their own (singleton) niche.
#'
#' @param spe a SpatialExperiment with niche reducedDims (see [buildNiches()]).
#' @param groups a named list mapping each merged niche name to a character
#'   vector of cell-type column names to sum.
#' @param sigma a numeric vector of bandwidths to update; \code{NULL} (default)
#'   updates every \code{name<sigma>} reducedDim present.
#' @param name a character, the reducedDim name prefix (default "Niche").
#' @param ... ignored.
#'
#' @return the input \code{spe} with merged niche reducedDims.
#'
#' @examples
#' data(toySpiDE)
#' spe <- toySpiDE
#' spe <- buildNiches(spe, sigma = 20)
#' spe <- mergeNiches(spe, groups = list(AC = c("A", "C")), sigma = 20)
#' colnames(SingleCellExperiment::reducedDim(spe, "Niche20"))
#'
#' @rdname mergeNiches
#' @importFrom SingleCellExperiment reducedDim reducedDim<- reducedDimNames
#' @export
setMethod(
  "mergeNiches",
  signature = "ANY",
  definition = function(spe, groups, sigma = NULL, name = "Niche", ...) {
    stopifnot(is.list(groups), !is.null(names(groups)))

    rdn <- SingleCellExperiment::reducedDimNames(spe)
    targets <- if (is.null(sigma)) grep(sprintf("^%s", name), rdn, value = TRUE) else paste0(name, sigma)
    missing <- setdiff(targets, rdn)
    if (length(missing) > 0) {
      stop(sprintf("niche reducedDim(s) not found: %s", paste(missing, collapse = ", ")))
    }

    for (rd in targets) {
      mat <- SingleCellExperiment::reducedDim(spe, rd)
      # carry through any cell type not named in a group as its own niche
      leftover <- setdiff(colnames(mat), unlist(groups, use.names = FALSE))
      full_groups <- c(groups, stats::setNames(as.list(leftover), leftover))

      merged <- vapply(full_groups, function(cts) {
        cts <- intersect(cts, colnames(mat))
        if (length(cts) == 0) {
          return(rep(0, nrow(mat)))
        }
        rowSums(mat[, cts, drop = FALSE])
      }, numeric(nrow(mat)))
      rownames(merged) <- rownames(mat)
      SingleCellExperiment::reducedDim(spe, rd) <- merged
    }

    return(spe)
  }
)
