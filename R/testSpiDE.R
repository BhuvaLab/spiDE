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
#' The condition-specific test is carried out within a family of triplets
#' chosen without looking at the condition, so the choice is exact under any
#' relabelling of patients and the multiplicity budget is spent on triplets
#' with something for the condition to explain:
#' \itemize{
#'   \item \code{procedure = "pooled_or_heterogeneity"} (the slopes engine's
#'     default): the triplets that pass either filter below.
#'   \item \code{procedure = "heterogeneity"}: the triplets whose slopes vary
#'     between patients more than their sampling variances allow (Cochran's Q
#'     across patients, the statistic behind the between-patient variance;
#'     Benjamini-Hochberg over every triplet at \code{fdr}). A condition effect
#'     is such variation, including one whose two conditions' slopes cancel in
#'     the pooled slope. With few patients the pooled test passes little, and
#'     this filter carries the condition test.
#'   \item \code{procedure = "filtered"} (the sandwich engine's default, and
#'     the slopes engine's to spiDE 0.99.36): the triplets whose pooled test
#'     passes at \code{fdr}. It keeps a modest condition difference on a strong
#'     pooled slope, which the heterogeneity test, spread over every patient,
#'     can miss.
#'   \item \code{procedure = "all"}: every triplet.
#' }
#' The heterogeneity filter needs each patient's slopes, so it is the slopes
#' engine's only. The calibration vignette has the measurements behind the
#' default.
#'
#' The slopes engine's pooled test refers its moderated t statistic to limma's
#' degrees of freedom scaled by the Kish effective share of the patients (the
#' proportional rule), so a column whose weights rest on a few patients is not
#' credited with the others. To spiDE 0.99.36 the df were the smaller of limma's
#' df and the Kish effective number of patients less one, which is conservative
#' with few patients.
#'
#' @param object a [SpiDEFit-class] from [fitSpiDE()].
#' @param condition \code{NULL} (use the condition recorded by [fitSpiDE()],
#'   if any) or a character naming a two-level, patient-level colData column.
#'   With the sandwich engine it must be the condition the model was fitted
#'   with. The contrast is second level minus first (factor levels, else
#'   sorted values).
#' @param procedure \code{"pooled_or_heterogeneity"},
#'   \code{"heterogeneity"}, \code{"filtered"} or \code{"all"}, see Details.
#'   The default is \code{"pooled_or_heterogeneity"} for the slopes engine and
#'   \code{"filtered"} for the sandwich engine.
#' @param fdr the FDR level of the filter.
#' @param strata \code{NULL} (the fit's) or a patient-level column adjusted for
#'   in the slopes engine's condition test (e.g. slide, where the condition is
#'   confounded with it).
#' @param BPPARAM a BiocParallelParam; the slopes engine's condition-specific
#'   test is run in parallel over genes.
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
  function(object, condition = NULL, procedure = c("pooled_or_heterogeneity", "heterogeneity", "filtered", "all"),
           fdr = 0.05,
           strata = NULL, BPPARAM = BiocParallel::SerialParam(), ...) {
    .assertCurrent(object)
    procedure <- if (missing(procedure) || is.null(procedure)) {
      if (object@engine == "slopes") "pooled_or_heterogeneity" else "filtered"
    } else match.arg(procedure)
    checkProcedure(procedure, object@engine)
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
        ci <- .conditionInputs(pt, x$patients, trt_of, strata)
        tk <- .slopesTests(x, trt = ci$trt, strata = ci$st, BPPARAM = BPPARAM)
      } else {
        tk <- x$coef
        if (is.null(condition)) tk <- tk[tk$test == "pooled", , drop = FALSE]
        tk$n_patients <- length(x$patients)
        tk$p.heterogeneity <- NA_real_
        tk <- tk[, c("gene", "niche", "test", "estimate", "se", "t", "df", "p", "n_patients", "p.heterogeneity")]
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
