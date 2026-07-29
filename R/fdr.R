# Hierarchical FDR. Three nested BH corrections gate the results:
#   1. gene level  - the gene-level niche p-value (all ResponseNiche combined),
#   2. index level - per index cell type, within surviving genes,
#   3. niche level - per (gene, index) across individual niche cell types.
# Up (positive) and Down (negative) directions are tested separately at
# fdr/2 and merged (a gene/pair significant in both directions is "Both").
# Reproduces the sig_niche cascade in YTMA_nicheDE_v9.md, in base R.

#' Gene- and index-level FDR gating, with direction merging
#'
#' @param p.pos,p.neg genes x (1 + n_index) Cauchy-combined p-value matrices.
#' @param fdr the FDR threshold.
#' @return a data.frame with gene, ct_index, fdr.gene, fdr.index,
#'   DirectionGene, DirectionIndex (or NULL if nothing passes).
#' @importFrom stats p.adjust
#' @noRd
.geneIndexFDR <- function(p.pos, p.neg, fdr) {
  per_dir <- function(P, dir) {
    gene_p <- P[, 1]
    # doubled, capped-at-1 BH-adjusted p-value: the two-sided-equivalent FDR
    # for this direction, comparable directly against the user's `fdr` (which
    # ranges over (0, 1]) rather than against an uncapped fdr/2 pre-gate.
    fdr_gene <- pmin(stats::p.adjust(gene_p, method = "BH") * 2, 1)
    pass <- !is.na(fdr_gene) & fdr_gene <= fdr
    if (!any(pass)) {
      return(NULL)
    }
    idxmat <- P[pass, -1, drop = FALSE]
    genes <- rownames(idxmat)
    fdr_index <- t(apply(idxmat, 1, stats::p.adjust, method = "BH"))
    dim(fdr_index) <- dim(idxmat)
    dimnames(fdr_index) <- dimnames(idxmat)
    fdr_index <- pmin(fdr_index * 2, 1)
    keep <- which(fdr_index <= fdr, arr.ind = TRUE)
    if (nrow(keep) == 0) {
      return(NULL)
    }
    data.frame(
      gene = genes[keep[, 1]],
      ct_index = colnames(idxmat)[keep[, 2]],
      fdr.gene = fdr_gene[pass][keep[, 1]],
      fdr.index = fdr_index[keep],
      Direction = dir,
      stringsAsFactors = FALSE
    )
  }

  tab <- rbind(per_dir(p.pos, "Up"), per_dir(p.neg, "Down"))
  if (is.null(tab) || nrow(tab) == 0) {
    return(NULL)
  }

  key <- paste(tab$gene, tab$ct_index, sep = "\r")
  out <- lapply(split(seq_len(nrow(tab)), key), function(ix) {
    sub <- tab[ix, , drop = FALSE]
    best <- which.min(sub$fdr.index)
    n <- nrow(sub)
    dir <- if (n == 2) "Both" else sub$Direction[best]
    # fdr.gene/fdr.index are already doubled and capped at 1 per direction;
    # combining the two directions is a min, not a further doubling.
    data.frame(
      gene = sub$gene[1], ct_index = sub$ct_index[1],
      fdr.gene = min(sub$fdr.gene),
      fdr.index = min(sub$fdr.index),
      DirectionGene = dir, DirectionIndex = dir,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

#' Long per-(gene, index, niche, bandwidth) record of niche coefficients
#' @noRd
.nicheRecords <- function(fits, genes) {
  recs <- lapply(fits, function(f) {
    ct <- as.character(f@covtype)
    respcols <- grepl("Response", ct)
    rn <- ct[respcols] == "ResponseNiche"
    tmat <- f@t_stat[genes, rn, drop = FALSE]
    cm <- f@coefmap[respcols, , drop = FALSE][rn, , drop = FALSE]
    cols <- cm$covariate
    # per-column df aligned to the ResponseNiche subset (rn) of the tested cols;
    # a scalar/NULL @df is used as-is (recycled / normal reference).
    dfn <- if (is.null(f@df) || length(f@df) == 1L) f@df else f@df[rn]
    lower <- .ptByCol(tmat, dfn)
    pmat <- pmin(lower, 1 - lower)
    do.call(rbind, lapply(seq_along(cols), function(j) {
      data.frame(
        gene = genes,
        ct_index = cm$index[j],
        ct_niche = cm$niche[j],
        bandwidth = f@sigma,
        t = tmat[, j],
        coef = f@alpha[genes, cols[j]],
        p = pmat[, j],
        stringsAsFactors = FALSE
      )
    }))
  })
  do.call(rbind, recs)
}

#' Niche-level FDR on the surviving (gene, index) pairs
#'
#' @param fits list of SpiDEFit.
#' @param gi_tab output of .geneIndexFDR().
#' @param gene.w gene x bandwidth weights.
#' @param fdr the FDR threshold.
#' @return gi_tab expanded to (gene, ct_index, ct_niche) with fdr.niche,
#'   DirectionNiche, coef, t, bandwidth.max.
#' @importFrom stats p.adjust
#' @noRd
.nicheLevelFDR <- function(fits, gi_tab, gene.w, fdr) {
  genes <- unique(gi_tab$gene)
  sigmas <- vapply(fits, function(f) f@sigma, numeric(1))
  recs <- .nicheRecords(fits, genes)
  # restrict to surviving (gene, ct_index) pairs
  pair <- paste(recs$gene, recs$ct_index, sep = "\r")
  keep_pair <- paste(gi_tab$gene, gi_tab$ct_index, sep = "\r")
  recs <- recs[pair %in% keep_pair, , drop = FALSE]
  if (nrow(recs) == 0) {
    return(NULL)
  }

  # combine p across bandwidths per (gene, ct_index, ct_niche)
  trip <- paste(recs$gene, recs$ct_index, recs$ct_niche, sep = "\r")
  combined <- lapply(split(seq_len(nrow(recs)), trip), function(ix) {
    sub <- recs[ix, , drop = FALSE]
    g <- sub$gene[1]
    w <- gene.w[g, match(sub$bandwidth, sigmas)]
    pc <- .cauchyCombine(matrix(sub$p, nrow = 1), matrix(w, nrow = 1))
    best <- which.min(sub$p)
    data.frame(
      gene = g, ct_index = sub$ct_index[1], ct_niche = sub$ct_niche[1],
      p.niche = pc, bandwidth.max = sub$bandwidth[best],
      t = sub$t[best], coef = sub$coef[best],
      DirectionNiche = if (sub$t[best] > 0) "Up" else "Down",
      stringsAsFactors = FALSE
    )
  })
  combined <- do.call(rbind, combined)

  # BH per (gene, ct_index) across ct_niche, then gate
  pair2 <- paste(combined$gene, combined$ct_index, sep = "\r")
  combined$fdr.niche <- NA_real_
  for (k in split(seq_len(nrow(combined)), pair2)) {
    combined$fdr.niche[k] <-
      pmin(stats::p.adjust(combined$p.niche[k], "BH") * 2, 1)
  }
  combined[combined$fdr.niche <= fdr, , drop = FALSE]
}

#' Full hierarchical FDR: assemble the tidy results table
#'
#' @param fits list of SpiDEFit.
#' @param p.pos,p.neg Cauchy-combined p-value matrices.
#' @param gene.w gene x bandwidth weights.
#' @param fdr the FDR threshold.
#' @return the tidy results data.frame (possibly with zero rows).
#' @noRd
.hierarchicalFDR <- function(fits, p.pos, p.neg, gene.w, fdr) {
  empty <- data.frame(
    gene = character(), ct_index = character(), ct_niche = character(),
    bandwidth.max = numeric(), coef = numeric(), t = numeric(),
    DirectionGene = character(), DirectionIndex = character(),
    DirectionNiche = character(), fdr.gene = numeric(),
    fdr.index = numeric(), fdr.niche = numeric(),
    stringsAsFactors = FALSE
  )

  gi <- .geneIndexFDR(p.pos, p.neg, fdr)
  if (is.null(gi) || nrow(gi) == 0) {
    return(empty)
  }
  niche <- .nicheLevelFDR(fits, gi, gene.w, fdr)
  if (is.null(niche) || nrow(niche) == 0) {
    return(empty)
  }

  out <- merge(niche, gi, by = c("gene", "ct_index"))
  out <- out[, c(
    "gene", "ct_index", "ct_niche", "bandwidth.max", "coef", "t",
    "DirectionGene", "DirectionIndex", "DirectionNiche",
    "fdr.gene", "fdr.index", "fdr.niche"
  )]
  out <- out[order(out$fdr.niche, -abs(out$t)), , drop = FALSE]
  rownames(out) <- NULL
  out
}
