#' Fit spiDE's per-patient niche models
#'
#' Fits, for every index cell type, one negative binomial GLM per gene of its
#' expression on the local density of every other (niche) cell type, with an
#' intercept per patient so that each niche effect is a within-patient slope.
#' The slopes are then tested across patients by [testSpiDE()], pooled or
#' between conditions.
#'
#' Two engines estimate the between-patient error the tests need:
#' \describe{
#'   \item{\code{"slopes"} (default)}{fits the shared model once, without the
#'     condition, then gives every patient a one-step NB estimate of its own
#'     niche slopes with a within-patient spatial sandwich variance.
#'     [testSpiDE()] combines them across patients with weighted limma. It is
#'     the more powerful engine when effects vary between patients, and its
#'     fit never sees the condition, so any condition can be tested later.}
#'   \item{\code{"sandwich"}}{fits the condition-specific model directly (the
#'     condition x niche columns) and tests its coefficients with a
#'     patient-clustered CR2 sandwich and Bell-McCaffrey degrees of freedom. It
#'     is simpler, and strongest for effects shared by all patients, but needs
#'     the condition at fit time.}
#' }
#' Both are calibrated on the null grids of five spatial cohorts; see the
#' calibration vignette.
#'
#' @param spe a SpatialExperiment with raw counts and a niche reducedDim from
#'   [buildNiches()].
#' @param condition \code{NULL} or a character, the colData column holding a
#'   two-level, patient-level condition. The slopes engine only records it
#'   (it can also be given to [testSpiDE()]); the sandwich engine fits the
#'   condition-specific model only when it is given.
#' @param engine \code{"slopes"} or \code{"sandwich"}.
#' @param index \code{NULL} (every cell type with at least \code{min.cells}
#'   cells in at least \code{min.patients} patients) or a character vector of
#'   index cell types.
#' @param niche \code{NULL} (test every niche column) or a character vector of
#'   niche columns to test; the others stay in the model as adjustments. An
#'   index type is never tested against its own niche (see [mergeNiches()]).
#' @param sigma the niche bandwidth (needed only if \code{spe} carries several).
#' @param cell_type,sample_id the colData columns of the cell type and of the
#'   patient (the unit that is replicated: one sample per patient).
#' @param section \code{NULL} (= \code{sample_id}) or the colData column of the
#'   tissue section, used to lay the spatial tiles of the variance estimate.
#' @param covariates a character vector of cell-level colData columns to adjust
#'   for (centred; patient-level ones are absorbed by the patient intercepts).
#' @param depth how sequencing depth enters: \code{"loglib"} (default, the
#'   centred log library size as a covariate with a slope per gene) or
#'   \code{"none"}.
#' @param strata \code{NULL} or a patient-level colData column (e.g. slide or
#'   batch) whose niche slopes are adjusted for. Needed when the condition is
#'   confounded with slide; used by the sandwich engine here and by the slopes
#'   engine in [testSpiDE()].
#' @param assay the counts assay (raw integer counts).
#' @param name the niche reducedDim prefix.
#' @param min.cells minimum cells of an index type a patient must contribute.
#' @param min.patients minimum patients an index type needs.
#' @param min.detect minimum fraction of an index type's cells in which a gene
#'   is detected for it to be tested in that type.
#' @param tile side of the square tiles of the spatial variance estimate, in
#'   coordinate units (default 3 x the bandwidth).
#' @param min.tiles minimum tiles a patient needs for its own spatial variance.
#' @param BPPARAM a BiocParallelParam; genes are fitted in parallel.
#' @param verbose report progress.
#' @param ... unused.
#' @return a [SpiDEFit-class].
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = 20)
#' fit <- fitSpiDE(spe, index = "A", sigma = 20, min.patients = 6)
#' fit
#' @seealso [testSpiDE()], [spiDE()], [testNicheAbundance()]
#' @importFrom SummarizedExperiment assay colData
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom S4Vectors metadata
#' @rdname fitSpiDE
#' @export
setMethod(
  "fitSpiDE", "SpatialExperiment",
  function(spe, condition = NULL, engine = c("slopes", "sandwich"), index = NULL,
           niche = NULL, sigma = NULL, cell_type = "cell_type", sample_id = "sample_id",
           section = NULL, covariates = character(), depth = c("loglib", "none"),
           strata = NULL, assay = "counts", name = "Niche", min.cells = 10L,
           min.patients = 6L, min.detect = 0.1, tile = NULL, min.tiles = 5L,
           BPPARAM = BiocParallel::SerialParam(), verbose = TRUE, ...) {
    engine <- match.arg(engine)
    depth <- match.arg(depth)
    checkSPE(spe, assay = assay, cell_type = cell_type, sample_id = sample_id)
    if (!is.null(condition)) checkCondition(spe, condition)
    checkSample(spe, condition = condition, sample_id = sample_id, covariates = covariates)
    checkCovariates(spe, covariates, finite.only = TRUE)
    sigma <- .pickSigma(spe, sigma, name)
    checkNiche(spe, sigma, name)
    cd <- SummarizedExperiment::colData(spe)
    if (!is.null(strata) && !strata %in% colnames(cd)) {
      stop(sprintf("strata column '%s' not found in colData(spe)", strata), call. = FALSE)
    }
    if (!is.null(section) && !section %in% colnames(cd)) {
      stop(sprintf("section column '%s' not found in colData(spe)", section), call. = FALSE)
    }
    Y <- SummarizedExperiment::assay(spe, assay)
    checkCounts(Y, integer.only = TRUE)
    ct <- as.character(cd[[cell_type]])
    smp <- as.character(cd[[sample_id]])
    sec <- if (is.null(section)) smp else as.character(cd[[section]])
    NM <- as.matrix(SingleCellExperiment::reducedDim(spe, paste0(name, sigma)))
    xy <- SpatialExperiment::spatialCoords(spe)
    if (is.null(tile)) tile <- 3 * sigma
    tiles <- paste(sec, floor(xy[, 1] / tile), floor(xy[, 2] / tile), sep = "\r")
    cov <- .cellCovariates(cd, covariates, depth, Y)
    # a cell with no counts has no library size: it cannot enter a model with a
    # log-depth term, so it is left out (reported, not silently)
    usable <- if (depth == "loglib") is.finite(cov[, "loglib"]) else rep(TRUE, ncol(Y))
    if (!all(usable) && verbose) message(sprintf("fitSpiDE: %d cell(s) with no counts left out", sum(!usable)))
    patients <- .patientTable(cd, smp)
    # mergeNiches() records its groups per niche reducedDim
    group_map <- S4Vectors::metadata(spe)[["spiDE_niche_groups"]][[paste0(name, sigma)]]
    if (is.null(index)) {
      index <- names(which(vapply(split(smp, ct), function(s) sum(table(s) >= min.cells), numeric(1)) >= min.patients))
    }
    index <- intersect(index, unique(ct))
    if (!length(index)) stop("no index cell type has enough cells in enough patients", call. = FALSE)
    trt_of <- if (engine == "sandwich" && !is.null(condition)) .conditionCoding(patients, condition) else NULL
    str_of <- if (!is.null(strata)) stats::setNames(as.character(patients[[strata]]), patients$patient) else NULL
    fits <- list()
    for (k in index) {
      ik <- .indexCells(ifelse(usable, ct, NA_character_), smp, k, min.cells)
      pk <- unique(smp[ik])
      if (length(pk) < min.patients) {
        if (verbose) message(sprintf("fitSpiDE: skipping %s (%d patients)", k, length(pk)))
        next
      }
      nc <- .nicheColumns(NM, k, ik, niche, group_map)
      if (!length(nc$tested)) next
      det <- Matrix::rowMeans(Y[, ik, drop = FALSE] > 0)
      gk <- rownames(Y)[det >= min.detect]
      if (!length(gk)) next
      cols <- c(nc$tested, setdiff(nc$cols, nc$tested))
      L <- log1p(NM[ik, cols, drop = FALSE])
      pat <- factor(smp[ik])
      Yk <- Y[gk, ik, drop = FALSE]
      covk <- if (ncol(cov)) scale(cov[ik, , drop = FALSE], scale = FALSE) else NULL
      if (verbose) message(sprintf("fitSpiDE: %s engine, index %s: %d genes, %d cells, %d patients, %d niches",
                                   engine, k, length(gk), length(ik), nlevels(pat), length(nc$tested)))
      stk <- if (!is.null(str_of) && engine == "sandwich") factor(str_of[as.character(pat)]) else NULL
      des0 <- .indexDesign(L, covk, pat, strata = stk, tested = nc$tested)
      fit0 <- .fitIndexGLM(Yk, des0, BPPARAM = BPPARAM)
      if (engine == "slopes") {
        ps <- .patientSlopes(fit0, des0, Yk, L, tiles[ik], min.cells = min.cells,
                             min.tiles = min.tiles, BPPARAM = BPPARAM)
        fits[[k]] <- c(list(genes = gk, niches = nc$tested, adjusted = setdiff(nc$cols, nc$tested),
                            patients = des0$patients), ps,
                       list(factor = .patientFactor(ps$v_tile, ps$v_model),
                            mean_expr = Matrix::rowMeans(Yk), psi = fit0$psi))
      } else {
        coef <- cbind(test = "pooled", .sandwichCR2(fit0, des0, Yk, BPPARAM = BPPARAM))
        if (!is.null(trt_of)) {
          des1 <- .indexDesign(L, covk, pat, trt = unname(trt_of[as.character(pat)]), strata = stk,
                               tested = nc$tested)
          fit1 <- .fitIndexGLM(Yk, des1, BPPARAM = BPPARAM)
          coef <- rbind(coef, cbind(test = "condition", .sandwichCR2(fit1, des1, Yk, BPPARAM = BPPARAM)))
        }
        fits[[k]] <- list(genes = gk, niches = nc$tested, adjusted = setdiff(nc$cols, nc$tested),
                          patients = des0$patients, coef = coef, mean_expr = Matrix::rowMeans(Yk))
      }
    }
    if (!length(fits)) stop("no index cell type could be fitted", call. = FALSE)
    methods::new("SpiDEFit", engine = engine, sigma = sigma,
                 condition = if (is.null(condition)) character() else condition,
                 patients = patients, index = fits,
                 params = list(cell_type = cell_type, sample_id = sample_id, section = section,
                               covariates = covariates, depth = depth, strata = strata,
                               min.cells = min.cells, min.patients = min.patients,
                               min.detect = min.detect, tile = tile, min.tiles = min.tiles,
                               version = as.character(utils::packageVersion("spiDE"))))
  }
)

# The niche bandwidth: the only one present, or the one asked for.
.pickSigma <- function(spe, sigma, name) {
  rdn <- SingleCellExperiment::reducedDimNames(spe)
  have <- suppressWarnings(as.numeric(sub(sprintf("^%s", name), "", grep(sprintf("^%s[0-9.]+$", name), rdn, value = TRUE))))
  if (!length(have)) stop("no niche reducedDims found; run buildNiches() first", call. = FALSE)
  if (is.null(sigma)) {
    if (length(have) > 1L) {
      stop(sprintf("spe carries several niche bandwidths (%s): choose one with 'sigma'",
                   paste(sort(have), collapse = ", ")), call. = FALSE)
    }
    return(have)
  }
  if (length(sigma) != 1L) stop("'sigma' must be a single bandwidth", call. = FALSE)
  sigma
}

# Cell-level covariates: the centred log library size (depth = "loglib") and
# any user covariates, as a numeric matrix over all cells.
.cellCovariates <- function(cd, covariates, depth, Y) {
  cols <- list()
  if (depth == "loglib") {
    lib <- Matrix::colSums(Y)
    cols$loglib <- ifelse(lib > 0, log(lib), NA_real_)
  }
  for (cv in covariates) cols[[cv]] <- as.numeric(cd[[cv]])
  if (!length(cols)) return(matrix(0, nrow(cd), 0))
  do.call(cbind, cols)
}

# One row per patient: id, cell count, and every colData column constant
# within patients (so a condition or strata can be named later).
.patientTable <- function(cd, smp) {
  pats <- sort(unique(smp))
  out <- data.frame(patient = pats, ncells = as.integer(table(smp)[pats]), stringsAsFactors = FALSE)
  first <- match(pats, smp)
  for (cn in colnames(cd)) {
    x <- cd[[cn]]
    if (!is.atomic(x) || is.matrix(x)) next
    const <- tapply(as.character(x), smp, function(v) length(unique(v[!is.na(v)])) <= 1L)
    if (all(const[pats])) out[[cn]] <- if (is.factor(x)) droplevels(x[first]) else x[first]
  }
  out
}

# 0/1 per patient for a two-level condition: 1 = the second level (factor
# levels, else sorted values), so the contrast is "second - first".
.conditionCoding <- function(patients, condition) {
  x <- patients[[condition]]
  if (is.null(x)) stop(sprintf("condition '%s' is not a patient-level column", condition), call. = FALSE)
  lv <- if (is.factor(x)) levels(droplevels(x)) else sort(unique(as.character(x[!is.na(x)])))
  if (length(lv) != 2L) stop(sprintf("condition '%s' must have exactly two levels", condition), call. = FALSE)
  out <- stats::setNames(as.numeric(as.character(x) == lv[2]), patients$patient)
  attr(out, "levels") <- lv
  out
}
