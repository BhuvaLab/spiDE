#' Extract spiDE's test results
#'
#' @param object a [SpiDEResults-class].
#' @param test \code{NULL} (the condition-specific test if one was run, else
#'   the pooled test), \code{"condition"}, \code{"pooled"} or \code{"both"}.
#' @param fdr \code{NULL} (every tested triplet) or a level: keep rows with
#'   \code{q <= fdr}.
#' @param ... unused.
#' @return a data.frame, one row per (gene, index, niche, test), ordered by
#'   p-value: \code{estimate} (the pooled slope, or the difference in slope
#'   between conditions, on the log scale per unit log1p niche density),
#'   \code{se}, \code{t}, \code{df}, \code{p}, \code{q} (BH within the test's
#'   family), \code{in_family} and, for the slopes engine,
#'   \code{p.heterogeneity} (the label-free test that the triplet's slopes
#'   vary between patients, the filter of \code{procedure = "heterogeneity"}).
#'   Condition-specific rows outside the family have \code{q = NA}.
#' @examples
#' data(toySpiDE)
#' res <- spiDE(toySpiDE, condition = "condition", sigma = 30, index = "A",
#'              procedure = "all")
#' results(res, test = "pooled")
#' @rdname results
#' @export
setMethod("results", "SpiDEResults", function(object, test = NULL, fdr = NULL, ...) {
  .assertCurrent(object)
  tb <- object@table
  if (is.null(test)) test <- if (any(tb$test == "condition")) "condition" else "pooled"
  test <- match.arg(test, c("condition", "pooled", "both"))
  if (test != "both") tb <- tb[tb$test == test, , drop = FALSE]
  if (!is.null(fdr)) tb <- tb[!is.na(tb$q) & tb$q <= fdr, , drop = FALSE]
  tb <- tb[order(tb$p, na.last = TRUE), , drop = FALSE]
  rownames(tb) <- NULL
  tb
})

#' Each patient's niche slopes
#'
#' The slopes engine estimates every patient's own niche slopes (a one-step
#' negative binomial estimate from the shared fit); this returns them in long
#' form, e.g. to plot how a gene's slope varies between patients and
#' conditions.
#'
#' @param object a [SpiDEFit-class] or [SpiDEResults-class] from the slopes
#'   engine.
#' @param gene,index,niche optional character vectors to subset.
#' @param ... unused.
#' @return a data.frame: \code{gene}, \code{index}, \code{niche},
#'   \code{patient}, \code{slope}, \code{var_spatial} (the within-patient
#'   spatial sandwich variance), \code{var_model}, \code{ncells}, plus the
#'   patient-level colData columns recorded by [fitSpiDE()].
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = 30)
#' fit <- fitSpiDE(spe, index = "A", sigma = 30)
#' head(patientSlopes(fit, gene = "G1"))
#' @rdname patientSlopes
#' @export
setMethod("patientSlopes", "SpiDEFit", function(object, gene = NULL, index = NULL, niche = NULL, ...) {
  .assertCurrent(object)
  if (object@engine != "slopes") {
    stop("per-patient slopes are estimated by the slopes engine (fitSpiDE(engine = \"slopes\"))", call. = FALSE)
  }
  out <- list()
  for (k in intersect(if (is.null(index)) names(object@index) else index, names(object@index))) {
    x <- object@index[[k]]
    gi <- if (is.null(gene)) seq_along(x$genes) else which(x$genes %in% gene)
    ni <- if (is.null(niche)) seq_along(x$niches) else which(x$niches %in% niche)
    if (!length(gi) || !length(ni)) next
    g <- expand.grid(g = gi, s = seq_along(x$patients), j = ni)
    ix <- cbind(g$g, g$s, g$j)
    out[[k]] <- data.frame(gene = x$genes[g$g], index = k, niche = x$niches[g$j],
                           patient = x$patients[g$s], slope = x$beta[ix], var_spatial = x$v_tile[ix],
                           var_model = x$v_model[ix], ncells = unname(x$ncells[g$s]),
                           stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, out)
  if (is.null(d)) return(data.frame())
  pt <- object@patients
  d <- merge(d, pt[, setdiff(colnames(pt), "ncells"), drop = FALSE], by = "patient", all.x = TRUE, sort = FALSE)
  d <- d[order(d$index, d$gene, d$niche, d$patient), c("gene", "index", "niche", "patient",
                                                       setdiff(colnames(d), c("gene", "index", "niche", "patient")))]
  rownames(d) <- NULL
  d
})

#' @rdname patientSlopes
#' @export
setMethod("patientSlopes", "SpiDEResults", function(object, ...) {
  .assertCurrent(object)
  patientSlopes(object@fit, ...)
})

#' Each patient's intercept
#'
#' Every index type's shared fit has an intercept per patient and gene: the
#' patient's expression level of the gene in that type, at zero niche density
#' and a reference depth: the index type's mean log library size under
#' \code{depth = "loglib"} or \code{"nonlinear"}, each patient's own mean
#' under \code{"spatial_spline"}. The intercepts absorb
#' every difference between patients -- the condition's main effect,
#' composition, batch -- which is why the niche slopes are within-patient
#' slopes. [plotPatientEffects()] shows what they capture. Kept by
#' [fitSpiDE()] from spiDE 0.99.34 on, for both engines.
#'
#' @param object a [SpiDEFit-class] or [SpiDEResults-class].
#' @param gene,index optional character vectors to subset.
#' @param ... unused.
#' @return a data.frame: \code{gene}, \code{index}, \code{patient},
#'   \code{intercept} (log scale), \code{ncells} (the patient's cells of the
#'   index type), plus the patient-level colData columns recorded by
#'   [fitSpiDE()].
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = 30)
#' fit <- fitSpiDE(spe, index = "A", sigma = 30, verbose = FALSE)
#' head(patientIntercepts(fit, gene = "G1"))
#' @rdname patientIntercepts
#' @export
setMethod("patientIntercepts", "SpiDEFit", function(object, gene = NULL, index = NULL, ...) {
  .assertCurrent(object)
  out <- list()
  for (k in intersect(if (is.null(index)) names(object@index) else index, names(object@index))) {
    x <- object@index[[k]]
    checkIntercepts(x, "patientIntercepts()")
    gi <- if (is.null(gene)) seq_along(x$genes) else which(x$genes %in% gene)
    if (!length(gi)) next
    g <- expand.grid(g = gi, s = seq_along(x$patients))
    out[[k]] <- data.frame(gene = x$genes[g$g], index = k, patient = x$patients[g$s],
                           intercept = x$intercept[cbind(g$g, g$s)], ncells = unname(x$ncells[g$s]),
                           stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, out)
  if (is.null(d)) return(data.frame())
  pt <- object@patients
  d <- merge(d, pt[, setdiff(colnames(pt), "ncells"), drop = FALSE], by = "patient", all.x = TRUE, sort = FALSE)
  d <- d[order(d$index, d$gene, d$patient),
         c("gene", "index", "patient", setdiff(colnames(d), c("gene", "index", "patient")))]
  rownames(d) <- NULL
  d
})

#' @rdname patientIntercepts
#' @export
setMethod("patientIntercepts", "SpiDEResults", function(object, ...) {
  .assertCurrent(object)
  patientIntercepts(object@fit, ...)
})
