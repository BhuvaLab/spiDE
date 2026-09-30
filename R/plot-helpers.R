# Data builders shared by the plots (0.99.34). Anything shared across genes is
# computed over the whole gene family, as the tests compute it (invariant 2);
# counts are read only as the plotted genes' rows or through a sparse product
# (invariant 5).

# The slopes engine's pooled-test variance of each patient's slope
# (v_model x the per-patient factor) and the DerSimonian-Laird tau2 of each
# (gene, niche) over the whole gene family: the weights 1 / (v + tau2).
.slopeWeights <- function(xi) {
  vpool <- sweep(xi$v_model, c(2, 3), xi$factor, "*")
  list(v = vpool, tau = .dlTau2(xi$beta, vpool))
}

# NULL = the condition-specific test if the results have one, else the pooled
# test (as results() does).
.resolveTest <- function(x, test) {
  have <- unique(x@table$test)
  if (is.null(test)) test <- if ("condition" %in% have) "condition" else "pooled"
  test <- match.arg(test, c("pooled", "condition"))
  if (!test %in% have) stop("these results have no condition-specific test", call. = FALSE)
  test
}

# Each patient's group: its condition level (factor levels in contrast order),
# or one group when there is no condition.
.groupLevels <- function(fit, condition) {
  pt <- fit@patients
  if (!length(condition)) {
    return(stats::setNames(factor(rep("all patients", nrow(pt))), pt$patient))
  }
  lv <- attr(.conditionCoding(pt, condition), "levels")
  stats::setNames(factor(as.character(pt[[condition]]), lv), pt$patient)
}

.directionLabels <- function(x, test) {
  if (test == "pooled") return(c(up = "rises with density", down = "falls with density"))
  lv <- attr(.conditionCoding(x@fit@patients, x@condition), "levels")
  c(up = sprintf("steeper in %s", lv[2]), down = sprintf("steeper in %s", lv[1]))
}

# The table rows of one test for these genes, in the genes' order.
.tripletRows <- function(tab, test, gene, index, niche) {
  r <- tab[tab$test == test & tab$index == index & tab$niche == niche, , drop = FALSE]
  r[match(gene, r$gene), , drop = FALSE]
}

# Each index type's genes in k bands by mean expression (the calibration
# vignette's expression fifths).
.expressionBands <- function(fit, k = 5L) {
  labs <- paste0("Q", seq_len(k))
  labs[1] <- paste(labs[1], "(dimmest)")
  labs[k] <- paste(labs[k], "(brightest)")
  do.call(rbind, lapply(names(fit@index), function(i) {
    xi <- fit@index[[i]]
    r <- rank(xi$mean_expr, ties.method = "first") / length(xi$genes)
    data.frame(gene = xi$genes, index = i,
               band = cut(r, seq(0, 1, length.out = k + 1L), labels = labs, include.lowest = TRUE),
               stringsAsFactors = FALSE)
  }))
}

# Mean expression per 10,000 counts of these genes in every cell type, by a
# product with a sparse, library-scaled indicator: only genes x cell types is
# dense. attr "ncells": the cells of each type.
.cellTypeMeans <- function(spe, genes, assay = "counts", cell_type = "cell_type") {
  Y <- SummarizedExperiment::assay(spe, assay)
  lib <- Matrix::colSums(Y)
  ct <- as.character(SummarizedExperiment::colData(spe)[[cell_type]])
  ok <- which(lib > 0 & !is.na(ct))
  types <- sort(unique(ct[ok]))
  j <- match(ct[ok], types)
  ind <- Matrix::sparseMatrix(i = ok, j = j, x = 1e4 / lib[ok], dims = c(ncol(Y), length(types)))
  n <- tabulate(j, length(types))
  M <- as.matrix(Y[genes, , drop = FALSE] %*% ind) / rep(n, each = length(genes))
  dimnames(M) <- list(genes, types)
  attr(M, "ncells") <- stats::setNames(n, types)
  M
}

# For each row of a results table: the gene's expression in the niche type
# over its expression in the index type (CP10k, + pseudo on both). A merged
# niche (mergeNiches()) pools its member cell types' cells.
.nicheFold <- function(tab, spe, fit, assay = "counts", name = "Niche", pseudo = 0.01) {
  genes <- unique(tab$gene)
  M <- .cellTypeMeans(spe, genes, assay, fit@params$cell_type)
  n <- attr(M, "ncells")
  groups <- S4Vectors::metadata(spe)[["spiDE_niche_groups"]][[paste0(name, fit@sigma)]]
  typeExpr <- function(type) {
    if (type %in% colnames(M)) return(M[, type])
    m <- intersect(groups[[type]], colnames(M))
    if (!length(m)) return(rep(NA_real_, nrow(M)))
    as.numeric(M[, m, drop = FALSE] %*% n[m]) / sum(n[m])
  }
  types <- unique(c(tab$index, tab$niche))
  E <- matrix(unlist(lapply(types, typeExpr)), nrow(M), length(types), dimnames = list(genes, types))
  (E[cbind(tab$gene, tab$niche)] + pseudo) / (E[cbind(tab$gene, tab$index)] + pseudo)
}

# One index type's cells as the fit used them (the fit's patients, cells with
# counts): these genes' counts, library sizes, patients and the log1p density
# of one niche column.
.indexCellData <- function(spe, fit, index, niche, genes, assay = "counts", name = "Niche") {
  cd <- SummarizedExperiment::colData(spe)
  ct <- as.character(cd[[fit@params$cell_type]])
  smp <- as.character(cd[[fit@params$sample_id]])
  xi <- fit@index[[index]]
  Y <- SummarizedExperiment::assay(spe, assay)
  miss <- setdiff(genes, rownames(Y))
  if (length(miss)) stop(sprintf("gene(s) not in 'spe': %s", paste(miss, collapse = ", ")), call. = FALSE)
  lib <- Matrix::colSums(Y)
  ik <- which(!is.na(ct) & ct == index & smp %in% xi$patients & lib > 0)
  if (!length(ik) || !all(xi$patients %in% smp[ik])) {
    stop(sprintf("the fit's %s cells are not all in 'spe': is it the object the fit was made from?", index),
         call. = FALSE)
  }
  checkNiche(spe, fit@sigma, name)
  NM <- SingleCellExperiment::reducedDim(spe, paste0(name, fit@sigma))
  if (!niche %in% colnames(NM)) {
    stop(sprintf("niche '%s' is not a column of %s%s", niche, name, fit@sigma), call. = FALSE)
  }
  list(counts = Y[genes, ik, drop = FALSE], lib = lib[ik], patient = smp[ik],
       density = log1p(as.numeric(NM[ik, niche])))
}
