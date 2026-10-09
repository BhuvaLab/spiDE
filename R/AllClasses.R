# The two result classes of the per-patient engines (0.99.30).
#
# SpiDEFit holds, per index cell type, what fitSpiDE() estimated: the genes
# tested, the tested niche columns, and either each patient's niche slopes with
# their variances (the slopes engine) or the fixed-effects fit's tested
# coefficients with their cluster-robust standard errors (the sandwich
# engine). SpiDEResults holds testSpiDE()'s tidy table and the fit it came from.
#
# The class NAMES are the ones spiDE <= 0.99.22 used for its mixed model, whose
# objects saved with saveRDS() record package "spiDE" and therefore resolve to
# these new definitions when read. They have none of the new slots; .isLegacy()
# recognises them and every entry point stops with a pointer to
# spiDEmixed::readSpiDE(), the archived engine's reader.

#' The fitted per-patient niche models of spiDE
#'
#' Returned by [fitSpiDE()]; tested with [testSpiDE()].
#'
#' @slot engine a character, \code{"slopes"} or \code{"sandwich"}.
#' @slot sigma a numeric, the niche bandwidth the fit used.
#' @slot condition a character (length 0 or 1), the condition the fit was
#'   given (the sandwich engine fits condition-specific slopes only if it was
#'   given one; the slopes engine only records it).
#' @slot patients a data.frame, one row per patient (sample): its id, the
#'   number of cells it contributes, and every colData column that is constant
#'   within patients (so a condition can be named at [testSpiDE()] time).
#' @slot index a named list, one element per index cell type, holding the
#'   engine's estimates (see Details).
#' @slot params a list of the settings the fit used.
#'
#' @details For the slopes engine each element of \code{index} holds
#'   \code{genes}, \code{niches} (tested niche columns), \code{patients},
#'   \code{beta} (genes x patients x niches: each patient's slope),
#'   \code{v_tile} (the within-patient spatial sandwich variance of each slope),
#'   \code{v_model} (the model variance), \code{factor} (patients x niches: the
#'   per-patient variance factor), \code{pooled} (genes x niches: the shared
#'   fit's slope), \code{mean_expr} and \code{ncells}. For the sandwich engine
#'   it holds \code{genes}, \code{niches} and a data.frame \code{coef} of the
#'   tested coefficients (\code{test} = "pooled" or "condition") with their
#'   CR2 standard errors and Bell-McCaffrey degrees of freedom. Both hold
#'   \code{status} (genes fitted; genes left out because the solver could not
#'   fit them or their coefficients ran off, the latter also counted as
#'   \code{runaway}; and genes whose dispersion is on its search bound) and, under
#'   \code{depth = "spatial_spline"}, \code{depth_r2} (patients x niches: the
#'   \eqn{R^2} of each niche column on the patient's library-size spline, a
#'   diagnostic of how much of the niche covariate the spline could absorb).
#'   From spiDE 0.99.34 both also hold \code{intercept} (genes x patients:
#'   each patient's intercept in the shared, condition-free fit; see
#'   [patientIntercepts()]), \code{niche_mean} (patients x niches: the mean
#'   log1p density of each tested niche over the patient's index cells),
#'   \code{loglib_mean} (the patients' mean log library size) and
#'   \code{ncells} (the patients' index cells).
#' @return An object of class \code{SpiDEFit}, created by [fitSpiDE()]; its
#'   \code{$} accessor returns a slot, and \code{show()} prints a summary.
#' @exportClass SpiDEFit
setClass("SpiDEFit", representation(
  engine = "character", sigma = "numeric", condition = "character",
  patients = "data.frame", index = "list", params = "list"
))

#' The results of spiDE's niche tests
#'
#' Returned by [testSpiDE()] and [spiDE()]; read it with [results()].
#'
#' @slot table a data.frame, one row per (gene, index, niche, test).
#' @slot condition a character (length 0 or 1), the condition tested.
#' @slot contrast a character, e.g. \code{"Responder - Non-responder"}.
#' @slot procedure a character, the condition test's family:
#'   \code{"pooled_or_heterogeneity"}, \code{"heterogeneity"},
#'   \code{"filtered"} or \code{"all"}.
#' @slot fdr a numeric, the FDR level of the filter.
#' @slot pooled.df a character, the slopes engine's pooled-test degrees of
#'   freedom (\code{"proportional"} or \code{"capped"}; empty for the
#'   sandwich engine). Results saved by spiDE <= 0.99.36 lack the slot; they
#'   were tested with \code{"capped"}.
#' @slot fit the [SpiDEFit-class] tested.
#' @return An object of class \code{SpiDEResults}, created by [testSpiDE()];
#'   its \code{$} accessor returns a slot, and \code{show()} prints a summary.
#' @exportClass SpiDEResults
setClass("SpiDEResults", representation(
  table = "data.frame", condition = "character", contrast = "character",
  procedure = "character", fdr = "numeric", pooled.df = "character", fit = "SpiDEFit"
))

# The pooled-test df rule of a results object; results saved by spiDE <= 0.99.36
# have no slot for it and were tested with the capped df.
.pooledDfOf <- function(object) {
  if (methods::.hasSlot(object, "pooled.df") && length(object@pooled.df)) object@pooled.df else "capped"
}

# A spiDE <= 0.99.22 mixed-model object read with readRDS(): it carries the
# old legacy slots (covtype, W, alpha, fits, ...) and none of the new ones.
.isLegacy <- function(object) {
  a <- names(attributes(object))
  any(c("covtype", "coefmap", "fits", "results.celltype") %in% a) ||  # legacy slots
    !any(c("engine", "table") %in% a)
}

.legacyMessage <- paste(
  "this is a spiDE <= 0.99.22 mixed-model object, from the engine archived as",
  "the research package spiDEmixed (research/mixed in the spiDE repository).",
  "Read it with spiDEmixed::readSpiDE(path) and use spiDEmixed:: functions on it."
)

.assertCurrent <- function(object) {
  if (.isLegacy(object)) stop(.legacyMessage, call. = FALSE)
  invisible(TRUE)
}

setValidity("SpiDEFit", function(object) {
  if (.isLegacy(object)) return(.legacyMessage)
  if (!length(object@engine) || !object@engine %in% c("slopes", "sandwich")) {
    return("'engine' must be \"slopes\" or \"sandwich\"")
  }
  TRUE
})

setValidity("SpiDEResults", function(object) {
  if (.isLegacy(object)) return(.legacyMessage)
  need <- c("gene", "index", "niche", "test", "estimate", "p", "q")
  miss <- setdiff(need, colnames(object@table))
  if (length(miss)) return(paste("results table lacks column(s):", paste(miss, collapse = ", ")))
  TRUE
})

#' @importFrom methods show
setMethod("show", "SpiDEFit", function(object) {
  if (.isLegacy(object)) {
    cat("<legacy spiDE mixed-model fit>\n", .legacyMessage, "\n", sep = "")
    return(invisible(NULL))
  }
  cat(sprintf("SpiDEFit (%s engine), bandwidth %s\n", object@engine, format(object@sigma)))
  cat(sprintf("  %d patients; condition: %s\n", nrow(object@patients),
              if (length(object@condition)) object@condition else "none"))
  for (k in names(object@index)) {
    x <- object@index[[k]]
    cat(sprintf("  index %s: %d genes x %d niches (%s)\n", k, length(x$genes),
                length(x$niches), paste(x$niches, collapse = ", ")))
  }
  invisible(NULL)
})

setMethod("show", "SpiDEResults", function(object) {
  if (.isLegacy(object)) {
    cat("<legacy spiDE mixed-model results>\n", .legacyMessage, "\n", sep = "")
    return(invisible(NULL))
  }
  tb <- object@table
  cat(sprintf("SpiDEResults (%s engine)\n", object@fit@engine))
  for (tt in unique(tb$test)) {
    x <- tb[tb$test == tt & !is.na(tb$q), , drop = FALSE]
    lab <- if (tt == "condition") {
      sprintf("condition-specific (%s; %s family)", object@contrast, object@procedure)
    } else if (object@fit@engine == "slopes") {
      sprintf("pooled across patients (%s df)", .pooledDfOf(object))
    } else "pooled across patients"
    cat(sprintf("  %s: %d tests, %d at FDR %s\n", lab, nrow(x),
                sum(x$q <= object@fdr), format(object@fdr)))
  }
  invisible(NULL)
})

#' @rdname SpiDEFit-class
#' @param x a SpiDEFit or SpiDEResults.
#' @param name a slot name.
#' @export
setMethod("$", "SpiDEFit", function(x, name) {
  .assertCurrent(x)
  methods::slot(x, name)
})

#' @rdname SpiDEResults-class
#' @param x a SpiDEResults.
#' @param name a slot name.
#' @export
setMethod("$", "SpiDEResults", function(x, name) {
  .assertCurrent(x)
  methods::slot(x, name)
})
