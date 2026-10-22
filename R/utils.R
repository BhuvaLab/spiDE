# Shared helpers used across spiDE: cell-type name sanitising, the self-niche
# rule and the sparse per-group column sums.

# sanitise cell-type / niche names the same way as the analysis scripts
.sanitise <- function(x) gsub(" |-", ".", x)

#' Which niche columns are the index type's own niche (to be dropped)
#'
#' A covariate pairing an index cell type with its own niche density is
#' meaningless and is dropped. With merged niches, "its own" means the index is
#' a \emph{member} of the group the niche column represents, so all member
#' indices of a merged niche are dropped, not only the exact-name match. A
#' group's own NAME is also one of its members: mergeNiches() maps a merged
#' column to its FINE sub-labels, but colData$cell_type may already carry the
#' MERGED labels, in which case the index label IS the group name.
#'
#' @param index a single index cell type.
#' @param niches character vector of niche column names.
#' @param group_map \code{NULL}, or a named list mapping each (raw) merged niche
#'   column name to a character vector of its (raw) fine member cell types.
#'   Names and members are sanitised here before matching. A niche column absent
#'   from \code{group_map} has member set \code{{itself}}, so the test reduces to
#'   \code{index == niche} (the pre-merge behaviour).
#' @return a logical vector over \code{niches}, \code{TRUE} where the column is
#'   the index type's own niche.
#' @noRd
.isSelfNiche <- function(index, niches, group_map = NULL) {
  index <- .sanitise(index)
  nm <- .sanitise(niches)
  members <- if (is.null(group_map)) list() else
    stats::setNames(lapply(group_map, .sanitise), .sanitise(names(group_map)))
  vapply(nm, function(n) index == n || index %in% members[[n]], logical(1), USE.NAMES = FALSE)
}

# Per-group sums of the columns of a genes x cells matrix, as a product with a
# sparse group indicator, so a sparse or DelayedArray count matrix is never
# densified whole (invariant 5): only genes x groups is dense. `group` gives
# each column's group as an integer code (NA = the column is left out),
# `weights` a scalar or per-column multiplier. Returns genes x ngroups, with
# Y's rownames and no column names.
.groupSums <- function(Y, group, weights = 1, ngroups = max(group, na.rm = TRUE)) {
  ok <- which(!is.na(group))
  ind <- Matrix::sparseMatrix(i = ok, j = as.integer(group[ok]), x = rep_len(weights, length(group))[ok],
                              dims = c(length(group), ngroups))
  as.matrix(Y %*% ind)
}
