#' spiDE: context-specific, neighbourhood-dependent differential expression
#'
#' spiDE tests how gene expression within an *index* cell type changes with an
#' experimental *condition* as a function of the local density (*niche*) of
#' surrounding cell types, in spatial transcriptomics data. The workflow is:
#' build niche covariates ([buildNiches()]), fit a per-gene negative binomial
#' GLM over an interaction design ([fitSpiDE()]) using the SpaNorm engine, and
#' test neighbourhood effects with combined Wald statistics and a hierarchical
#' FDR ([testSpiDE()]). [spiDE()] wraps the three steps.
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
