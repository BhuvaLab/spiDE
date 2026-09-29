#' Gene-set tests on spiDE's niche results (experimental)
#'
#' For every (index cell type, niche cell type) column of a test, asks whether
#' the genes of a set change with the niche more than the other tested genes
#' (a competitive test, as in \code{limma::camera()}): per-gene signed
#' \eqn{z}-scores from the test's p-values are averaged over the set and
#' compared with the rest, with the set's variance inflated by
#' \eqn{1 + (m - 1)\rho} for inter-gene correlation \eqn{\rho}. With the slopes
#' engine \eqn{\rho} is the mean correlation between genes of the per-patient
#' slopes in that column; otherwise it is \code{rho} (default 0).
#'
#' \strong{Experimental.} This replaces the mixed model's \code{spiGSEA()}
#' (archived in spiDEmixed); it has not yet passed a gene-set null benchmark
#' on the new engines, so its calls should be treated as exploratory.
#'
#' @param object a [SpiDEResults-class].
#' @param genesets a named list of character vectors of gene names.
#' @param test \code{NULL} (the condition-specific test if run, else pooled),
#'   \code{"condition"} or \code{"pooled"}.
#' @param min.size,max.size set size limits after intersecting with the tested
#'   genes of each index type.
#' @param rho \code{NULL} (estimate it, slopes engine) or a fixed inter-gene
#'   correlation.
#' @param ... unused.
#' @return a data.frame: \code{set}, \code{index}, \code{niche}, \code{test},
#'   \code{size}, \code{mean_z}, \code{t}, \code{p}, \code{q} (BH within each
#'   (index, niche) column), \code{direction} and \code{rho}.
#' @examples
#' data(toySpiDE)
#' res <- spiDE(toySpiDE, condition = "condition", sigma = 20, index = "A",
#'              min.patients = 6, procedure = "all")
#' sets <- list(first = paste0("G", 1:5), second = paste0("G", 6:12))
#' spiGSEA(res, sets, test = "pooled", min.size = 3)
#' @rdname spiGSEA
#' @export
setMethod("spiGSEA", "SpiDEResults", function(object, genesets, test = NULL, min.size = 5L,
                                             max.size = 500L, rho = NULL, ...) {
  .assertCurrent(object)
  if (!is.list(genesets) || is.null(names(genesets))) stop("'genesets' must be a named list", call. = FALSE)
  tb <- results(object, test = test)
  tt <- unique(tb$test)
  fit <- object@fit
  out <- list()
  for (k in unique(tb$index)) {
    for (n in unique(tb$niche[tb$index == k])) {
      x <- tb[tb$index == k & tb$niche == n & is.finite(tb$p), , drop = FALSE]
      if (nrow(x) < 3L) next
      z <- sign(x$t) * stats::qnorm(pmax(x$p, 1e-300) / 2, lower.tail = FALSE)
      names(z) <- x$gene
      r <- if (!is.null(rho)) rho else .slopeCorrelation(fit, k, n, names(z))
      G <- length(z)
      for (s in names(genesets)) {
        idx <- which(names(z) %in% genesets[[s]])
        m <- length(idx)
        if (m < min.size || m > max.size || m >= G) next
        vz <- stats::var(z[-idx])
        if (!is.finite(vz) || vz <= 0) vz <- 1
        vif <- 1 + (m - 1) * r
        d <- mean(z[idx]) - mean(z[-idx])
        tstat <- d / sqrt(vz * (vif / m + 1 / (G - m)))
        out[[length(out) + 1L]] <- data.frame(set = s, index = k, niche = n, test = tt, size = m,
                                              mean_z = mean(z[idx]), t = tstat,
                                              p = 2 * stats::pt(-abs(tstat), G - 2),
                                              direction = if (tstat > 0) "up" else "down", rho = r,
                                              stringsAsFactors = FALSE)
      }
    }
  }
  d <- do.call(rbind, out)
  if (is.null(d)) return(data.frame())
  d$q <- stats::ave(d$p, d$index, d$niche, FUN = function(p) stats::p.adjust(p, "BH"))
  d <- d[order(d$p), c("set", "index", "niche", "test", "size", "mean_z", "t", "p", "q", "direction", "rho")]
  rownames(d) <- NULL
  d
})

# Mean inter-gene correlation of the per-patient slopes in one (index, niche)
# column (slopes engine); 0 for the sandwich engine. Genes are standardised
# across patients and rho = (G var(mean z) - 1) / (G - 1), the mean
# off-diagonal correlation.
.slopeCorrelation <- function(fit, k, n, genes) {
  if (fit@engine != "slopes") return(0)
  x <- fit@index[[k]]
  j <- match(n, x$niches)
  if (is.na(j)) return(0)
  b <- matrix(x$beta[, , j], nrow = length(x$genes), dimnames = list(x$genes, NULL))
  b <- b[intersect(genes, rownames(b)), , drop = FALSE]
  keep <- colSums(!is.finite(b)) == 0
  b <- b[, keep, drop = FALSE]
  if (nrow(b) < 3L || ncol(b) < 3L) return(0)
  zs <- t(scale(t(b)))
  zs <- zs[stats::complete.cases(zs), , drop = FALSE]
  G <- nrow(zs)
  if (G < 3L) return(0)
  r <- (G * stats::var(colMeans(zs)) - 1) / (G - 1)
  max(0, min(r, 1))
}
