# Shared helpers used across spiDE: cell-type name sanitising and the
# self-niche rule. Both moved here from R/design.R when the mixed-effects
# engine was archived (0.99.30); the bodies are unchanged.

# sanitise cell-type / niche names the same way as the analysis scripts
.sanitise <- function(x) gsub(" |-", ".", x)

#' Which interaction columns are index-vs-own-niche (to be dropped)
#'
#' A covariate pairing an index cell type with its own niche density is
#' meaningless and is dropped. With merged niches, "its own" means the index is
#' a \emph{member} of the group the niche column represents, so all member
#' indices of a merged niche are dropped, not only the exact-name match.
#'
#' @param index,niche equal-length character vectors of the sanitised parsed
#'   index / niche cell type of each design column (\code{NA} for non
#'   interaction columns).
#' @param group_map \code{NULL}, or a named list mapping each (raw) merged niche
#'   column name to a character vector of its (raw) fine member cell types.
#'   Names and members are sanitised here before matching. A niche column absent
#'   from \code{group_map} has member set \code{{itself}}, so the test reduces to
#'   \code{index == niche} — the pre-merge behaviour.
#' @return a logical vector, \code{TRUE} where the column should be dropped.
#' @noRd
.isSelfNiche <- function(index, niche, group_map = NULL) {
  san_map <- if (is.null(group_map)) {
    list()
  } else {
    stats::setNames(lapply(group_map, .sanitise), .sanitise(names(group_map)))
  }
  # The group's own NAME is one of its members: mergeNiches() maps a merged
  # column to its FINE sub-labels, but colData$cell_type may already carry the
  # MERGED labels, in which case the index label IS the group name and would
  # otherwise never match its own members.
  members <- function(nch) unique(c(nch, san_map[[nch]]))

  ok <- !is.na(index) & !is.na(niche)
  res <- logical(length(index))
  res[ok] <- vapply(which(ok), function(k) index[k] %in% members(niche[k]),
                    logical(1))
  res
}
