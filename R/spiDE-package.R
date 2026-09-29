#' spiDE: context-specific, neighbourhood-dependent differential expression
#'
#' spiDE asks, in spatial transcriptomics data, whether a gene's expression in
#' an \emph{index} cell type changes with the local density (the \emph{niche})
#' of another cell type, within patients, and whether that dependence differs
#' between two conditions. Build niche covariates with [buildNiches()], fit the
#' per-patient niche models with [fitSpiDE()] (the slopes or sandwich engine),
#' and test them with [testSpiDE()]: the pooled niche test, and the
#' condition-specific test when a condition is given. [spiDE()] chains the
#' three. [testNicheAbundance()] asks the different, between-patient question
#' of whether patients with more of a niche type express genes differently.
#'
#' The mixed-effects engine of spiDE <= 0.99.22 is archived as the research
#' package spiDEmixed (research/mixed in the spiDE repository); read objects
#' saved by it with \code{spiDEmixed::readSpiDE()}.
#'
#' @keywords internal
#' @name spiDE-package
#' @aliases spiDE-package
"_PACKAGE"

# Symbols used in non-standard evaluation (data.frame/model.matrix building)
# to keep R CMD check quiet about undefined globals.
utils::globalVariables(c(
  "sample_id", "cell_type", "CellType", "Response",
  "ct_index", "ct_niche", "gene", "bandwidth", "value"
))
