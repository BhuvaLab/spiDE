# Design-matrix construction for the neighbourhood-dependent DE model. Builds
#   ~ 0 + <covariates> + CellType + <condition> * CellType:(niches) + niches
# then tags each column by type and drops symmetric self-interactions
# (an index cell type interacting with its own niche density). Reproduces the
# design in batch_nichede_v9.R; column tagging is done by parsing tokens so it
# is independent of the order R assigns to interaction labels.

# sanitise cell-type / niche names the same way as the analysis scripts
.sanitise <- function(x) gsub(" |-", ".", x)

#' Tag each design-matrix column by covariate type and parse index/niche cells
#'
#' @param cols character vector of column names.
#' @param niche_cols character vector of (sanitised) niche column names.
#' @param response_coef the condition main-effect coefficient name.
#' @return a data.frame with columns covariate, type, index, niche.
#' @noRd
.tagCovtype <- function(cols, niche_cols, response_coef) {
  parse_one <- function(col) {
    tokens <- strsplit(col, ":", fixed = TRUE)[[1]]
    ct_tok <- tokens[grepl("^CellType", tokens)]
    niche_tok <- intersect(tokens, niche_cols)
    has_resp <- response_coef %in% tokens

    index <- if (length(ct_tok) == 1) sub("^CellType", "", ct_tok) else NA_character_
    niche <- if (length(niche_tok) == 1) niche_tok else NA_character_

    if (length(tokens) == 1) {
      type <- if (col == response_coef) {
        "Response"
      } else if (grepl("^CellType", col)) {
        "CellType"
      } else {
        "Other"
      }
      return(c(type = type, index = NA_character_, niche = NA_character_))
    }
    # interaction columns
    if (length(ct_tok) == 1 && length(niche_tok) == 1) {
      type <- if (has_resp) "ResponseNiche" else "Niche"
      return(c(type = type, index = index, niche = niche))
    }
    c(type = "Other", index = NA_character_, niche = NA_character_)
  }

  parsed <- t(vapply(cols, parse_one, character(3)))
  data.frame(
    covariate = cols,
    type = parsed[, "type"],
    index = parsed[, "index"],
    niche = parsed[, "niche"],
    stringsAsFactors = FALSE
  )
}

#' Build the neighbourhood-interaction design for one bandwidth
#'
#' @param spe a SpatialExperiment with a niche reducedDim for \code{sigma}.
#' @param condition a character, the colData column of the tested condition.
#' @param sigma a numeric, the bandwidth (a single value).
#' @param index,niche character vectors restricting index / niche cell types
#'   (NULL = all).
#' @param covariates a character vector of nuisance colData columns.
#' @param cell_type a character, the colData column of cell type labels.
#' @param name a character, the niche reducedDim prefix.
#' @return a list with `W` (design matrix), `covtype` (factor), `coefmap`
#'   (data.frame), and `response_coef` (character).
#' @importFrom stats model.matrix as.formula
#' @importFrom SingleCellExperiment reducedDim
#' @noRd
.buildNicheDesign <- function(spe, condition, sigma, index = NULL, niche = NULL,
                              covariates = character(), cell_type = "cell_type",
                              name = "Niche") {
  cd <- SummarizedExperiment::colData(spe)

  # niche matrix -> log1p, sanitised column names
  nichemat <- SingleCellExperiment::reducedDim(spe, paste0(name, sigma))
  colnames(nichemat) <- .sanitise(colnames(nichemat))
  nichemat <- log1p(as.matrix(nichemat))

  niche_cols <- colnames(nichemat)
  if (!is.null(niche)) {
    niche_cols <- intersect(.sanitise(niche), niche_cols)
    if (length(niche_cols) == 0) {
      stop("no requested niche cell types found in the niche reducedDim")
    }
  }

  # assemble the model data.frame
  df <- as.data.frame(nichemat[, niche_cols, drop = FALSE])
  df[["CellType"]] <- factor(.sanitise(as.character(cd[[cell_type]])))
  df[[condition]] <- factor(cd[[condition]])
  for (cv in covariates) {
    df[[cv]] <- cd[[cv]]
  }

  # tested (non-reference) level of the condition -> coefficient name
  tested_level <- levels(df[[condition]])[2]
  response_coef <- paste0(condition, tested_level)

  # formula: 0 + covariates + CellType + condition*CellType:(niches) + niches
  niche_f <- paste(niche_cols, collapse = " + ")
  terms <- c(
    covariates,
    "CellType",
    sprintf("%s * CellType:(%s)", condition, niche_f),
    niche_f
  )
  f <- stats::as.formula(paste("~ 0 +", paste(terms, collapse = " + ")))
  W <- stats::model.matrix(f, df)

  # tag columns and drop symmetric self-interactions (index == niche)
  coefmap <- .tagCovtype(colnames(W), niche_cols, response_coef)
  is_self <- !is.na(coefmap$index) & !is.na(coefmap$niche) &
    coefmap$index == coefmap$niche
  keep <- !is_self

  # optionally restrict interaction columns to requested index cell types
  # (main effects and non-interaction columns are always kept)
  if (!is.null(index)) {
    idx_san <- .sanitise(index)
    is_interaction <- coefmap$type %in% c("Niche", "ResponseNiche")
    keep <- keep & (!is_interaction | coefmap$index %in% idx_san)
  }

  W <- W[, keep, drop = FALSE]
  coefmap <- coefmap[keep, , drop = FALSE]
  covtype <- factor(
    coefmap$type,
    levels = c("CellType", "Niche", "Response", "ResponseNiche", "Other")
  )

  list(W = W, covtype = covtype, coefmap = coefmap, response_coef = response_coef)
}

#' Build a spiDE design matrix
#'
#' Constructs the neighbourhood-interaction design matrix for a single niche
#' bandwidth and tags each covariate by type ("CellType", "Niche", "Response",
#' "ResponseNiche", or "Other"). The scientifically important covariates are the
#' three-way `CellType:condition:niche` interactions ("ResponseNiche"), which
#' capture how expression within an index cell type changes with the condition
#' as a function of a niche cell type's local density. Symmetric self
#' interactions (an index cell type against its own density) are dropped. This
#' is an escape hatch for custom fits; most users should call [fitSpiDE()].
#'
#' @inheritParams .buildNicheDesign
#' @param ... ignored.
#'
#' @return a list with `W` (the design matrix), `covtype` (a factor of column
#'   types), and `coefmap` (a data.frame mapping columns to index/niche cells).
#'
#' @examples
#' data(toySpiDE)
#' spe <- toySpiDE
#' spe <- buildNiches(spe, sigma = 20)
#' des <- nicheDesign(spe, condition = "condition", sigma = 20)
#' table(des$covtype)
#'
#' @rdname nicheDesign
#' @export
nicheDesign <- function(spe, condition, sigma, index = NULL, niche = NULL,
                        covariates = character(), cell_type = "cell_type",
                        name = "Niche", ...) {
  checkSPE(spe, cell_type = cell_type)
  checkCondition(spe, condition)
  checkCovariates(spe, covariates)
  checkNiche(spe, sigma, name = name)
  res <- .buildNicheDesign(spe, condition, sigma, index, niche, covariates,
                           cell_type, name)
  res[c("W", "covtype", "coefmap")]
}
