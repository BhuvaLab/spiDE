# Design-matrix construction for the neighbourhood-dependent DE model. Builds
#   ~ 0 + <covariates> + CellType + <condition> * CellType:(niches) + niches
# or, for a condition-free (niche-only) design (condition = NULL),
#   ~ 0 + <covariates> + CellType + CellType:(niches)
# then tags each column by type and drops symmetric self-interactions
# (an index cell type interacting with its own niche density). Reproduces the
# design in batch_nichede_v9.R; column tagging is done by parsing tokens so it
# is independent of the order R assigns to interaction labels.

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
  members <- function(nch) if (!is.null(san_map[[nch]])) san_map[[nch]] else nch

  ok <- !is.na(index) & !is.na(niche)
  res <- logical(length(index))
  res[ok] <- vapply(which(ok), function(k) index[k] %in% members(niche[k]),
                    logical(1))
  res
}

#' Build the patient-level random-effect design block
#'
#' Constructs the columns \code{Z} that carry the sample (patient) random
#' effects used to correct cell-level pseudo-replication (see the mixed-effects
#' vignette). Random effects are represented as ordinary design columns that are
#' L2-penalised at fit time (the ridge = random-effects equivalence).
#'
#' The random effects target the *response-related* fixed effects, which are the
#' only ones of inferential interest: for every response coefficient we add its
#' patient-level random counterpart, so the response contrast is tested against
#' between-patient variability. Concretely the random intercept (one indicator
#' per sample) is the patient-level counterpart of the response main effect, and
#' each random slope (a sample indicator times a fixed \code{CellType:niche} base
#' column) is the patient-level counterpart of the corresponding
#' response x \code{CellType:niche} coefficient.
#'
#' @param sample_vec a factor/character of sample ids, length = ncells.
#' @param slope_base a cells x k matrix of the fixed \code{CellType:niche} base
#'   columns (the non-response bases of the ResponseNiche terms); the random
#'   slopes are built from these when \code{random == "slope"}.
#' @param random one of "intercept" or "slope".
#' @return a list with \code{Z} (the random-effect design block) and
#'   \code{re_group} (a character label per column: "SampleInt" / "SampleSlope").
#' @noRd
.buildRandomEffects <- function(sample_vec, slope_base, random) {
  smp <- factor(sample_vec)
  Zint <- stats::model.matrix(~ 0 + smp)
  colnames(Zint) <- paste0("Sample", levels(smp))
  Z <- Zint
  re_group <- rep("SampleInt", ncol(Zint))

  if (random == "slope" && ncol(slope_base) > 0) {
    Zslope <- do.call(cbind, lapply(colnames(slope_base), function(bc) {
      M <- Zint * slope_base[, bc]
      colnames(M) <- paste0(colnames(Zint), ":", bc)
      M
    }))
    Z <- cbind(Zint, Zslope)
    re_group <- c(re_group, rep("SampleSlope", ncol(Zslope)))
  }
  list(Z = Z, re_group = re_group)
}

#' Tag each design-matrix column by covariate type and parse index/niche cells
#'
#' @param cols character vector of column names.
#' @param niche_cols character vector of (sanitised) niche column names.
#' @param response_coef the condition main-effect coefficient name.
#' @return a data.frame with columns covariate, type, index, niche.
#' @noRd
.tagCovtype <- function(cols, niche_cols, response_coef = NULL) {
  # response_coef = NULL is the condition-free design: there is no response
  # token, so has_resp is uniformly FALSE and the existing parse rules fall
  # through to "Niche" for the CellType:niche interactions -- which are the
  # tested terms in that mode.
  parse_one <- function(col) {
    tokens <- strsplit(col, ":", fixed = TRUE)[[1]]
    ct_tok <- tokens[grepl("^CellType", tokens)]
    niche_tok <- intersect(tokens, niche_cols)
    has_resp <- !is.null(response_coef) && response_coef %in% tokens

    index <- if (length(ct_tok) == 1) sub("^CellType", "", ct_tok) else NA_character_
    niche <- if (length(niche_tok) == 1) niche_tok else NA_character_

    if (length(tokens) == 1) {
      type <- if (!is.null(response_coef) && col == response_coef) {
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
    # E2 fix: CellType x condition with no niche token. Tagged "ResponseCellType"
    # deliberately -- fitSpiDE selects tested columns with
    # grepl("Response", covtype), so these DO receive standard errors and
    # t-statistics, while .nicheRecords() subsets on covtype == "ResponseNiche"
    # exactly, so they are excluded from the 3-level triplet FDR cascade and
    # leave it unchanged.
    if (length(ct_tok) == 1 && length(niche_tok) == 0 && has_resp) {
      return(c(type = "ResponseCellType", index = index, niche = NA_character_))
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
#' @param condition a character, the colData column of the tested condition, or
#'   \code{NULL} for a condition-free (niche-only) design.
#' @param sigma a numeric, the bandwidth (a single value).
#' @param index,niche character vectors restricting index / niche cell types
#'   (NULL = all).
#' @param covariates a character vector of nuisance colData columns.
#' @param cell_type a character, the colData column of cell type labels.
#' @param name a character, the niche reducedDim prefix.
#' @param sample_id a character, the colData column identifying samples
#'   (patients); used only when \code{random != "none"}.
#' @param random one of "none" (fixed-effects design, the default), "intercept"
#'   (add a per-sample random intercept) or "slope" (also add per-sample random
#'   slopes on the niche covariates). The random-effect columns are penalised at
#'   fit time to implement a mixed model via ridge (see the vignette).
#' @return a list with `W` (design matrix), `covtype` (factor), `coefmap`
#'   (data.frame), `response_coef` (character, `NULL` without a condition),
#'   `re_group` (per-column random-effect group label, `NA` for fixed columns),
#'   and `mode` ("condition" or "niche").
#' @importFrom stats model.matrix as.formula
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom S4Vectors metadata
#' @noRd
.buildNicheDesign <- function(spe, condition = NULL, sigma, index = NULL,
                              niche = NULL,
                              covariates = character(), cell_type = "cell_type",
                              name = "Niche", sample_id = "sample_id",
                              random = c("none", "intercept", "slope")) {
  random <- match.arg(random)
  has_cond <- !is.null(condition)
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
  if (has_cond) {
    df[[condition]] <- factor(cd[[condition]])
  }
  for (cv in covariates) {
    df[[cv]] <- cd[[cv]]
  }

  # tested (non-reference) level of the condition -> coefficient name
  response_coef <- if (has_cond) {
    paste0(condition, levels(df[[condition]])[2])
  } else {
    NULL
  }

  # formula: 0 + covariates + CellType + CellType:condition +
  #          CellType:(niches) + CellType:condition:(niches) + niches
  # Cell-means coding: no bare `condition` main effect -- the response
  # contrast lives in the CellType:condition columns, one per cell type.
  niche_f <- paste(niche_cols, collapse = " + ")
  # E2 fix: the released formula expands to condition + CellType:niches +
  # condition:CellType:niches but NO CellType:condition, so a cell-type-specific
  # response INTERCEPT shift aliases onto the three-way term (the niche
  # covariates are uncentred, so a coefficient there can mimic a constant).
  # Writing the terms explicitly and OMITTING the condition main effect makes the
  # CellType:condition block cell-means coded: one coefficient per cell type,
  # each the responder-vs-non-responder shift WITHIN that cell type, readable
  # without a contrast. Keeping the main effect would instead give k-1
  # treatment-coded columns requiring condition + CellType_x:condition.
  # Without a condition (condition = NULL) the condition groups drop out AND so
  # do the bare niche main effects, leaving
  #   ~ 0 + covariates + CellType + CellType:(niches)
  # so that the CellType:(niches) block -- tagged "Niche" either way -- becomes
  # the tested set under CELL-MEANS coding, exactly as CellType:condition is
  # cell-means coded above. c() drops the NULL entries.
  #
  # Dropping the niche main effects is required, not cosmetic. With them in,
  # niche_n = sum_c CellType_c:niche_n, so model.matrix() silently drops one
  # interaction per niche as aliased -- always the FIRST cell type's. In
  # condition mode that is harmless, because the tested terms are the three-way
  # CellType:condition:niche columns, which are not aliased with anything. Here
  # the two-way columns ARE the tested terms, so the alphabetically-first cell
  # type would lose every one of its niche slopes and could never be tested
  # (nicheDesign(index = "A") returned an empty design). Without the main
  # effects no aliasing occurs, every index cell type keeps a slope against
  # every non-self niche, and each coefficient is directly the within-cell-type
  # slope of expression on that niche's log density -- readable without a
  # contrast.
  terms <- c(
    covariates,
    "CellType",
    if (has_cond) sprintf("CellType:%s", condition),
    sprintf("CellType:(%s)", niche_f),
    if (has_cond) sprintf("CellType:%s:(%s)", condition, niche_f),
    if (has_cond) niche_f
  )
  f <- stats::as.formula(paste("~ 0 +", paste(terms, collapse = " + ")))
  W <- stats::model.matrix(f, df)

  # tag columns and drop self-interactions: an index cell type against its own
  # niche density. When niches have been merged (mergeNiches records the group
  # membership in metadata), "its own" means the index is a member of the merged
  # niche's group, so every member index of a merged niche is dropped. Absent a
  # stored mapping this reduces to the exact index == niche drop.
  coefmap <- .tagCovtype(colnames(W), niche_cols, response_coef)
  group_map <- S4Vectors::metadata(spe)[["spiDE_niche_groups"]][[paste0(name, sigma)]]
  is_self <- .isSelfNiche(coefmap$index, coefmap$niche, group_map)
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
  re_group <- rep(NA_character_, ncol(W))

  # append the patient random-effect block (penalised at fit time); these carry
  # the mixed-effects correction for cell-level pseudo-replication. The slopes
  # sit on the fixed CellType:niche columns (the non-response bases of the
  # ResponseNiche terms) so they target the response-related effects.
  if (random != "none") {
    if (!sample_id %in% colnames(cd)) {
      stop(sprintf("sample id column '%s' not found in colData(spe)", sample_id))
    }
    slope_base <- W[, coefmap$type == "Niche", drop = FALSE]
    re <- .buildRandomEffects(cd[[sample_id]], slope_base, random)
    W <- cbind(W, re$Z)
    coefmap <- rbind(coefmap, data.frame(
      covariate = colnames(re$Z), type = "Random",
      index = NA_character_, niche = NA_character_, stringsAsFactors = FALSE
    ))
    re_group <- c(re_group, re$re_group)
  }

  # E2 fix: "ResponseCellType" must be a declared level, otherwise factor()
  # silently maps the new columns to NA and cols_tested becomes NA for them.
  lvls <- c("CellType", "Niche", "Response", "ResponseNiche",
            "ResponseCellType", "Other")
  if (random != "none") lvls <- c(lvls, "Random")
  covtype <- factor(coefmap$type, levels = lvls)

  list(W = W, covtype = covtype, coefmap = coefmap,
       response_coef = response_coef, re_group = re_group,
       mode = if (has_cond) "condition" else "niche")
}

#' Build a spiDE design matrix
#'
#' Constructs the neighbourhood-interaction design matrix for a single niche
#' bandwidth and tags each covariate by type ("CellType", "Niche", "Response",
#' "ResponseNiche", or "Other"). The scientifically important covariates are the
#' three-way `CellType:condition:niche` interactions ("ResponseNiche"), which
#' capture how expression within an index cell type changes with the condition
#' as a function of a niche cell type's local density. Self interactions — an
#' index cell type against its own niche density — are dropped; when niches have
#' been merged with [mergeNiches()], this extends to any index that is a member
#' of the merged niche's group. This is an escape hatch for custom fits; most
#' users should call [fitSpiDE()].
#'
#' With `condition = NULL` the condition terms are omitted and the two-way
#' `CellType:niche` interactions ("Niche") become the tested effects. That
#' design also omits the bare niche main effects, so the interaction block is
#' cell-means coded: each coefficient is directly the slope of expression on
#' that niche cell type's log density *within* the index cell type, rather than
#' a contrast against a reference cell type.
#'
#' @param spe a SpatialExperiment with a niche reducedDim for \code{sigma}.
#' @param condition a character, the colData column of the tested condition, or
#'   \code{NULL} (default) for a condition-free design. With \code{NULL} the
#'   condition terms are omitted and the two-way \code{CellType:niche}
#'   interactions (tagged "Niche") become the tested effects; no \code{Response}
#'   columns are produced.
#' @param sigma a numeric, the bandwidth (a single value).
#' @param index,niche character vectors restricting index / niche cell types
#'   (NULL = all).
#' @param covariates a character vector of nuisance colData columns.
#' @param cell_type a character, the colData column of cell type labels.
#' @param name a character, the niche reducedDim prefix.
#' @param sample_id a character, the colData column identifying samples
#'   (patients); used only when \code{random != "none"}.
#' @param random one of "none" (fixed-effects design, the default), "intercept"
#'   (add a per-sample random intercept) or "slope" (also add per-sample random
#'   slopes on the niche covariates). The random-effect columns are penalised at
#'   fit time to implement a mixed model via ridge (see the vignette).
#' @param ... ignored.
#'
#' @return a list with `W` (the design matrix), `covtype` (a factor of column
#'   types), `coefmap` (a data.frame mapping columns to index/niche cells), and
#'   `mode` ("condition" or "niche").
#'
#' @examples
#' data(toySpiDE)
#' spe <- toySpiDE
#' spe <- buildNiches(spe, sigma = 20)
#' des <- nicheDesign(spe, condition = "condition", sigma = 20)
#' table(des$covtype)
#'
#' des0 <- nicheDesign(spe, condition = NULL, sigma = 20)
#' table(des0$covtype)
#'
#' @rdname nicheDesign
#' @export
nicheDesign <- function(spe, condition = NULL, sigma, index = NULL,
                        niche = NULL,
                        covariates = character(), cell_type = "cell_type",
                        name = "Niche", sample_id = "sample_id",
                        random = c("none", "intercept", "slope"), ...) {
  random <- match.arg(random)
  checkSPE(spe, cell_type = cell_type)
  if (!is.null(condition)) checkCondition(spe, condition)
  checkCovariates(spe, covariates)
  checkNiche(spe, sigma, name = name)
  res <- .buildNicheDesign(spe, condition, sigma, index, niche, covariates,
                           cell_type, name, sample_id, random)
  keep <- c("W", "covtype", "coefmap", "mode",
            if (random != "none") "re_group")
  res[keep]
}
