#' Gene-set tests on spiDE's niche results
#'
#' Asks, for every gene set and every (index cell type, niche cell type)
#' column, whether the set's genes as a group change with the niche (the
#' pooled test) or change with it differently between the two conditions (the
#' condition-specific test).
#'
#' A gene set is tested exactly like a gene. Every patient gets a slope for the
#' set: the average of its genes' niche slopes in that patient, each gene's
#' slope first put on a common scale (divided by the gene's typical spread of
#' slopes between patients, so a highly expressed gene does not dominate). By
#' default the set's slope is \emph{competitive}: the same average over the
#' column's other tested genes is subtracted, so the test asks whether the set
#' moves \emph{more} than the typical gene. This guards against a shift that
#' moves every gene, e.g. from segmentation spillover. With \code{type =
#' "self-contained"} nothing is subtracted.
#'
#' The patients' set slopes then go through the same pooled and condition
#' tests as the genes' slopes (see [testSpiDE()] and the model vignette).
#' Because the patient is the unit of replication, the correlation between the
#' genes of a set is accounted for without being estimated: it is already in
#' how much the set's slope varies from patient to patient.
#'
#' Needs the slopes engine, whose fit holds each patient's slopes. For gene sets
#' in the between-patient test, use [testNicheAbundance()] with \code{genesets}.
#'
#' @param object a [SpiDEResults-class] from the slopes engine.
#' @param genesets a named list of character vectors of gene names.
#' @param test \code{"pooled"}, \code{"condition"} or both (the default: every
#'   test the results hold).
#' @param type \code{"competitive"} (default) or \code{"self-contained"}.
#' @param min.size,max.size set size limits after intersecting with the genes
#'   tested in each index type.
#' @param BPPARAM a BiocParallelParam; the condition-specific test is run in
#'   parallel over sets.
#' @param ... unused.
#' @return a data.frame, one row per (set, index, niche, test):
#'   \code{size} (genes of the set tested there), \code{estimate} (the mean set
#'   slope, pooled, or its difference between conditions, in units of the
#'   genes' between-patient spread), \code{se}, \code{t}, \code{df}, \code{p},
#'   \code{q} (Benjamini-Hochberg within the (index, niche, test) column),
#'   \code{q.global} (over every row of that test), \code{n_patients} and
#'   \code{direction}.
#' @examples
#' data(toySpiDE)
#' res <- spiDE(toySpiDE, condition = "condition", sigma = 30, index = "A",
#'              procedure = "all")
#' sets <- list(first = paste0("G", 1:5), second = paste0("G", 6:12))
#' spiGSEA(res, sets, min.size = 3)
#' @seealso [testSpiDE()]
#' @rdname spiGSEA
#' @export
setMethod("spiGSEA", "SpiDEResults", function(object, genesets, test = NULL,
                                             type = c("competitive", "self-contained"),
                                             min.size = 5L, max.size = 500L,
                                             BPPARAM = BiocParallel::SerialParam(), ...) {
  .assertCurrent(object)
  type <- match.arg(type)
  checkGenesets(genesets)
  fit <- object@fit
  if (fit@engine != "slopes") {
    stop("spiGSEA() needs the slopes engine's per-patient slopes: fit with engine = \"slopes\"",
         call. = FALSE)
  }
  have <- if (length(object@condition)) c("pooled", "condition") else "pooled"
  if (is.null(test)) test <- have
  test <- match.arg(test, c("pooled", "condition"), several.ok = TRUE)
  if (!all(test %in% have)) {
    stop("the condition-specific set test needs results tested with a condition", call. = FALSE)
  }
  pt <- fit@patients
  trt_of <- if ("condition" %in% test) .conditionCoding(pt, object@condition) else NULL
  strata <- fit@params$strata
  out <- list()
  for (k in names(fit@index)) {
    x <- fit@index[[k]]
    trt <- if (!is.null(trt_of)) unname(trt_of[x$patients]) else NULL
    st <- if (!is.null(trt) && !is.null(strata)) {
      stats::setNames(as.character(pt[[strata]]), pt$patient)[x$patients]
    } else NULL
    member <- lapply(genesets, function(g) which(x$genes %in% g))
    m <- lengths(member)
    member <- member[m >= min.size & m <= max.size & m < length(x$genes)]
    if (!length(member)) next
    ss <- .setSlopes(x, member, type)
    for (j in seq_along(x$niches)) {
      b <- ss$b[, , j, drop = FALSE]
      v <- ss$v[, , j, drop = FALSE]
      tau <- .dlTau2(b, v)[, 1]
      b <- matrix(b, nrow = length(member)); v <- matrix(v, nrow = length(member))
      if ("pooled" %in% test) {
        r <- .pooledColumnTest(b, v, tau, ss$mean_expr, trend = length(member) >= 20L)
        out[[length(out) + 1L]] <- data.frame(set = names(member), index = k, niche = x$niches[j],
                                              test = "pooled", size = lengths(member), r,
                                              stringsAsFactors = FALSE)
      }
      if ("condition" %in% test) {
        r <- .robustConditionTest(b, v, tau, trt = trt, strata = st, BPPARAM = BPPARAM)
        out[[length(out) + 1L]] <- data.frame(set = names(member), index = k, niche = x$niches[j],
                                              test = "condition", size = lengths(member), r,
                                              stringsAsFactors = FALSE)
      }
    }
  }
  d <- do.call(rbind, out)
  if (is.null(d)) return(data.frame())
  d$q <- NA_real_; d$q.global <- NA_real_
  ok <- is.finite(d$p)
  col <- paste(d$index, d$niche, d$test, sep = "\r")
  d$q[ok] <- stats::ave(d$p[ok], col[ok], FUN = function(p) stats::p.adjust(p, "BH"))
  d$q.global[ok] <- stats::ave(d$p[ok], d$test[ok], FUN = function(p) stats::p.adjust(p, "BH"))
  d$direction <- ifelse(is.finite(d$t), ifelse(d$t > 0, "up", "down"), NA_character_)
  d <- d[order(d$test != "pooled", d$p), c("set", "index", "niche", "test", "size", "estimate", "se", "t", "df",
                                          "p", "q", "q.global", "n_patients", "direction")]
  rownames(d) <- NULL
  d
})

# Each set's per-patient slope for one index type, for every niche column.
# A gene's slopes are put on a common scale by sigma_g = sqrt(tau2_g + median_s
# v_gs), its typical spread of slopes between patients (the DerSimonian-Laird
# between-patient variance plus a typical within-patient variance, both from
# the pooled weights and computed without the condition). The set's slope in a
# patient is the mean of its genes' scaled slopes there (genes missing in that
# patient left out; the patient drops when under half the set is present),
# minus, when competitive, the mean over the column's other genes. Its
# variance treats the genes' sampling errors as independent within the
# patient; that only sets the patients' relative weights, since both tests
# take their error from the spread of the set slopes between patients.
.setSlopes <- function(x, member, type) {
  G <- length(x$genes); S <- length(x$patients); nt <- length(x$niches); K <- length(member)
  vpool <- sweep(x$v_model, c(2, 3), x$factor, "*")
  tau <- .dlTau2(x$beta, vpool)
  B <- V <- array(NA_real_, c(K, S, nt))
  for (j in seq_len(nt)) {
    bj <- matrix(x$beta[, , j], G, S); vj <- matrix(vpool[, , j], G, S)
    ok <- is.finite(bj) & is.finite(vj)
    bj[!ok] <- NA; vj[!ok] <- NA
    sig <- sqrt(tau[, j] + apply(vj, 1, stats::median, na.rm = TRUE))
    sig[!is.finite(sig) | sig <= 0] <- NA
    u <- bj / sig; vu <- vj / sig^2
    fin <- is.finite(u) & is.finite(vu)
    u[!fin] <- 0; vu[!fin] <- 0
    tot_n <- colSums(fin); tot_b <- colSums(u); tot_v <- colSums(vu)
    for (s in seq_len(K)) {
      g <- member[[s]]
      n_in <- colSums(fin[g, , drop = FALSE])
      sb <- colSums(u[g, , drop = FALSE]); sv <- colSums(vu[g, , drop = FALSE])
      mb <- sb / n_in
      mv <- sv / n_in^2
      if (type == "competitive") {
        n_out <- tot_n - n_in
        mb <- mb - (tot_b - sb) / n_out
        mv <- mv + (tot_v - sv) / n_out^2
      }
      drop <- n_in < length(g) / 2
      if (type == "competitive") drop <- drop | n_out < 1
      mb[drop] <- NA; mv[drop] <- NA
      B[s, , j] <- mb; V[s, , j] <- mv
    }
  }
  me <- vapply(member, function(g) exp(mean(log(x$mean_expr[g] + 1e-3))), numeric(1))
  list(b = B, v = V, mean_expr = me)
}
