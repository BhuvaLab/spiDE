#' Test niche-dependent differential expression
#'
#' Always runs the \strong{pooled} test: for every (gene, index cell type, niche
#' cell type), whether the gene's expression in the index type changes with the
#' density of the niche type, within patients, consistently across patients.
#' If a \code{condition} is given (here, or recorded by [fitSpiDE()]), it also
#' runs the \strong{condition-specific} test: whether that slope differs
#' between the two conditions -- a between-patient contrast of within-patient
#' slopes.
#'
#' By default (\code{procedure = "filtered"}) the condition-specific test is
#' carried out only among triplets whose pooled test passes at \code{fdr}: the
#' filter never looks at the condition, so it is exact under any relabelling of
#' patients, and it spends the multiplicity budget on triplets with a niche
#' effect to modify. \code{procedure = "all"} tests every triplet.
#'
#' @param object a [SpiDEFit-class] from [fitSpiDE()].
#' @param condition \code{NULL} (use the condition recorded by [fitSpiDE()],
#'   if any) or a character naming a two-level, patient-level colData column.
#'   With the sandwich engine it must be the condition the model was fitted
#'   with. The contrast is second level minus first (factor levels, else
#'   sorted values).
#' @param procedure \code{"filtered"} or \code{"all"}, see Details.
#' @param fdr the FDR level of the filter.
#' @param strata \code{NULL} (the fit's) or a patient-level column adjusted for
#'   in the slopes engine's condition test (e.g. slide, where the condition is
#'   confounded with it).
#' @param ... unused.
#' @return a [SpiDEResults-class]; read it with [results()].
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = 30)
#' fit <- fitSpiDE(spe, index = "A", sigma = 30)
#' res <- testSpiDE(fit, condition = "condition", procedure = "all")
#' res
#' head(results(res))
#' @rdname testSpiDE
#' @export
setMethod(
  "testSpiDE", "SpiDEFit",
  function(object, condition = NULL, procedure = c("filtered", "all"), fdr = 0.05,
           strata = NULL, ...) {
    .assertCurrent(object)
    procedure <- match.arg(procedure)
    checkFdr(fdr)
    if (is.null(condition) && length(object@condition)) condition <- object@condition
    if (object@engine == "sandwich" && !is.null(condition) &&
        !identical(condition, object@condition)) {
      stop("the sandwich engine tests the condition it was fitted with; refit with ",
           sprintf("fitSpiDE(condition = \"%s\", engine = \"sandwich\")", condition), call. = FALSE)
    }
    if (is.null(strata)) strata <- object@params$strata
    pt <- object@patients
    trt_of <- if (!is.null(condition)) .conditionCoding(pt, condition) else NULL
    contrast <- if (!is.null(trt_of)) paste(rev(attr(trt_of, "levels")), collapse = " - ") else character()
    tabs <- list()
    for (k in names(object@index)) {
      x <- object@index[[k]]
      if (object@engine == "slopes") {
        trt <- if (!is.null(trt_of)) unname(trt_of[x$patients]) else NULL
        st <- if (!is.null(trt) && !is.null(strata)) {
          stats::setNames(as.character(pt[[strata]]), pt$patient)[x$patients]
        } else NULL
        tk <- .slopesTests(x, trt = trt, strata = st)
      } else {
        tk <- x$coef
        if (is.null(condition)) tk <- tk[tk$test == "pooled", , drop = FALSE]
        tk$n_patients <- length(x$patients)
        tk <- tk[, c("gene", "niche", "test", "estimate", "se", "t", "df", "p", "n_patients")]
      }
      tabs[[k]] <- cbind(index = k, tk, stringsAsFactors = FALSE)
    }
    tab <- do.call(rbind, tabs)
    tab <- tab[, c("gene", "index", "niche", setdiff(colnames(tab), c("gene", "index", "niche")))]
    rownames(tab) <- NULL
    tab <- .bhFamilies(tab, procedure = procedure, fdr = fdr)
    methods::new("SpiDEResults", table = tab,
                 condition = if (is.null(condition)) character() else condition,
                 contrast = contrast, procedure = procedure, fdr = fdr, fit = object)
  }
)
