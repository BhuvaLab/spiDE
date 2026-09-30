# The niche views (0.99.34): the niche covariate in tissue, and the index
# cells' expression against it. Design: design/specs/2026-10-01-plots.md.

#' The niche covariate in tissue
#'
#' One sample's cells, coloured by the log1p density of a niche type: the
#' covariate spiDE regresses on. With \code{index}, only the index type's
#' cells are coloured -- the covariate as the index type's model reads it --
#' and the others are light grey.
#'
#' Read it to judge the bandwidth: a speckled field that follows single cells
#' is too small; a field nearly flat within a section is too large (the
#' patient intercept then absorbs it). It also shows where the index cells sit
#' relative to the niche.
#'
#' @param spe a SpatialExperiment with niche reducedDims ([buildNiches()]).
#' @param niche a niche column (a cell type or a [mergeNiches()] group).
#' @param sample the sample (patient) to show.
#' @param index \code{NULL} or an index cell type to highlight.
#' @param sigma the bandwidth (needed only if \code{spe} carries several).
#' @param cell_type,sample_id the colData columns of cell type and sample.
#' @param section \code{NULL} or the colData column of the section; a sample
#'   with several sections gets one panel each.
#' @param name the niche reducedDim prefix.
#' @param point.size the point size (default: from the number of cells).
#' @param scale.bar draw a scale bar of about a fifth of the width.
#' @param unit the coordinates' unit, for the scale bar's label (e.g.
#'   \code{"µm"}).
#' @return a ggplot; \code{p$data} has one row per cell of the sample:
#'   \code{x}, \code{y}, \code{cell_type}, \code{density} (log1p),
#'   \code{section} and \code{focal} (coloured).
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = 30)
#' plotNicheMap(spe, niche = "B", sample = "S1", index = "A")
#' @export
plotNicheMap <- function(spe, niche, sample, index = NULL, sigma = NULL, cell_type = "cell_type",
                         sample_id = "sample_id", section = NULL, name = "Niche", point.size = NULL,
                         scale.bar = TRUE, unit = NULL) {
  checkSPE(spe, cell_type = cell_type, sample_id = sample_id)
  checkColumn(spe, section, "section")
  sigma <- .pickSigma(spe, sigma, name)
  checkNiche(spe, sigma, name)
  NM <- SingleCellExperiment::reducedDim(spe, paste0(name, sigma))
  if (!is.character(niche) || length(niche) != 1L || !niche %in% colnames(NM)) {
    stop(sprintf("niche '%s' is not a column of %s%s", paste(niche, collapse = ", "), name, sigma), call. = FALSE)
  }
  cd <- SummarizedExperiment::colData(spe)
  smp <- as.character(cd[[sample_id]])
  if (length(sample) != 1L || !sample %in% smp) {
    stop(sprintf("sample '%s' not found in colData(spe)$%s", paste(sample, collapse = ", "), sample_id),
         call. = FALSE)
  }
  i <- which(smp == sample)
  ct <- as.character(cd[[cell_type]])[i]
  if (!is.null(index) && !index %in% ct) stop(sprintf("no %s cells in sample %s", index, sample), call. = FALSE)
  xy <- SpatialExperiment::spatialCoords(spe)[i, , drop = FALSE]
  d <- data.frame(x = xy[, 1], y = xy[, 2], cell_type = ct, density = log1p(as.numeric(NM[i, niche])),
                  section = if (is.null(section)) sample else as.character(cd[[section]])[i],
                  stringsAsFactors = FALSE)
  d$focal <- if (is.null(index)) rep(TRUE, nrow(d)) else d$cell_type == index
  if (is.null(point.size)) point.size <- min(1.5, max(0.05, 50 / sqrt(nrow(d))))
  p <- ggplot(d, aes(.data$x, .data$y))
  if (!is.null(index)) {
    p <- p + geom_point(data = function(d) d[!d$focal, ], colour = .spideCols[["faint"]],
                        size = point.size * 0.7, shape = 16)
  }
  p <- p + geom_point(data = function(d) d[d$focal, ], aes(colour = .data$density), size = point.size, shape = 16) +
    scale_colour_gradientn(colours = .magnitudeRamp(256L, skip = 0.1), name = sprintf("%s density\n(log1p)", niche)) +
    coord_equal() + theme_spiDE() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
  if (length(unique(d$section)) > 1L) p <- p + facet_wrap(~ section)
  if (isTRUE(scale.bar)) p <- p + .scaleBar(d, unit)
  p
}

# A scale bar of about a fifth of the width, under the bottom-right corner,
# with room below the cells for its label.
.scaleBar <- function(d, unit) {
  xr <- range(d$x)
  yr <- range(d$y)
  len <- signif(diff(xr) / 5, 1)
  bar <- data.frame(x = xr[2] - len, xend = xr[2], y = yr[1] - 0.04 * diff(yr),
                    label = if (is.null(unit)) format(len) else paste(format(len), unit))
  list(geom_segment(data = bar, aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y),
                    linewidth = 0.8, inherit.aes = FALSE),
       geom_text(data = bar, aes(x = (.data$x + .data$xend) / 2, y = .data$y, label = .data$label),
                 vjust = 1.6, size = 3.2, inherit.aes = FALSE),
       scale_y_continuous(expand = expansion(mult = c(0.16, 0.04))))
}

#' Expression against niche density: the cells behind a slope
#'
#' The index type's cells' expression of each gene (counts per 10,000, log
#' axis) against the log1p density of the niche type, binned, one panel per
#' gene.
#'
#' With \code{adjust = TRUE} (the default) each cell's expression is divided by
#' its patient's fitted intercept, so the patients share one reference level
#' and overlay: faint points are each patient's mean in a density bin, points
#' and bars the mean across patients and its 95\% interval, per condition, and
#' the dashed line each group's fitted slope. This is what the tests compare.
#' With \code{adjust = FALSE} the data are raw: each patient's own slope is a
#' segment through its mean (open circle), and the dashed line through the
#' patients' means is the between-patient association that the intercepts
#' absorb -- [testNicheAbundance()]'s question, not a niche slope.
#'
#' Read it for the data behind a slope, whether the log1p scale is roughly
#' linear, and whether a few extreme cells drive it. The display normalises
#' depth by library size (the model fits a slope on it). Every layer carries
#' \code{patient}, so \code{+ facet_wrap(~ patient)} splits it by patient.
#'
#' @param spe the SpatialExperiment the fit was made from.
#' @param x a [SpiDEResults-class] or [SpiDEFit-class]; \code{adjust = TRUE}
#'   needs one from spiDE >= 0.99.34. Under the sandwich engine every
#'   patient's and group's line has the pooled slope.
#' @param gene one or more genes, one panel each.
#' @param index,niche the index and niche cell types.
#' @param adjust divide out each patient's intercept (see Details).
#' @param bins the number of density bins (\code{adjust = TRUE}; bins are
#'   quantiles over all patients). Raw data use fifths of each patient's own
#'   density range.
#' @param assay the counts assay.
#' @param name the niche reducedDim prefix.
#' @return a ggplot; \code{p$data} has one row per (gene, patient, bin):
#'   \code{gene}, \code{patient}, \code{bin}, \code{expr} (mean CP10k),
#'   \code{density} (mean log1p density), \code{condition} and
#'   \code{gene_lab}. Bins without counts are left out of that layer (not out
#'   of the means).
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = 30)
#' res <- spiDE(spe, condition = "condition", sigma = 30, index = "A", verbose = FALSE)
#' plotNicheResponse(spe, res, gene = "G1", index = "A", niche = "B")
#' plotNicheResponse(spe, res, gene = "G1", index = "A", niche = "B", adjust = FALSE)
#' @export
plotNicheResponse <- function(spe, x, gene, index, niche, adjust = TRUE, bins = 8L, assay = "counts",
                              name = "Niche") {
  checkResults(x, "plotNicheResponse()", fit.ok = TRUE)
  fit <- if (is(x, "SpiDEResults")) x@fit else x
  condition <- if (is(x, "SpiDEResults")) x@condition else fit@condition
  checkSPE(spe, assay = assay, cell_type = fit@params$cell_type, sample_id = fit@params$sample_id)
  gene <- unique(gene)
  checkTriplet(fit, index, niche, gene)
  xi <- fit@index[[index]]
  if (adjust) {
    checkIntercepts(xi, "plotNicheResponse(adjust = TRUE)")
    # a gene the fit could not fit has no intercept: it drops out, not the plot
    fitted <- rowSums(is.finite(xi$intercept[gene, , drop = FALSE])) > 0
    if (!all(fitted)) {
      message(sprintf("plotNicheResponse(): %s not fitted in %s, left out", paste(gene[!fitted], collapse = ", "), index))
    }
    gene <- gene[fitted]
    if (!length(gene)) stop("no gene left to plot: none of them was fitted", call. = FALSE)
  }
  cells <- .indexCellData(spe, fit, index, niche, gene, assay, name)
  grp <- .groupLevels(fit, condition)
  lv <- levels(grp)
  lab <- stats::setNames(.itGene(gene), gene)
  sl <- .responseSlopes(fit, xi, gene, niche)
  l <- cells$density
  pat <- cells$patient
  if (adjust) {
    br <- unique(stats::quantile(l, seq(0, 1, length.out = bins + 1L)))
    bin <- if (length(br) < 2L) rep(1L, length(l)) else cut(l, br, include.lowest = TRUE, labels = FALSE)
  } else {
    bin <- stats::ave(l, pat, FUN = function(v) {
      br <- unique(stats::quantile(v, seq(0, 1, length.out = 6L)))
      if (length(br) < 2L) rep(1, length(v)) else cut(v, br, include.lowest = TRUE, labels = FALSE)
    })
  }
  pb <- do.call(rbind, lapply(gene, function(g) {
    e <- as.numeric(cells$counts[g, ]) / cells$lib * 1e4
    if (adjust) {
      a <- xi$intercept[g, ]
      # the reference is the patients' median intercept: one patient pinned at an
      # extreme value by the intercepts' ridge would move a mean
      e <- e * exp(-(a - stats::median(a[is.finite(a)])))[pat]
    }
    d <- stats::aggregate(cbind(expr = e, density = l) ~ patient + bin,
                          data.frame(e = e, l = l, patient = pat, bin = bin), mean)
    d$gene <- g
    d
  }))
  pb$condition <- grp[pb$patient]
  pb$gene_lab <- factor(lab[pb$gene], unname(lab))
  p <- ggplot(pb, aes(.data$density, .data$expr))
  if (adjust) {
    sm <- do.call(rbind, lapply(split(pb, list(pb$gene, pb$condition, pb$bin), drop = TRUE), function(s) {
      n <- nrow(s)
      se <- if (n > 1L) stats::sd(s$expr) / sqrt(n) else NA_real_
      data.frame(gene = s$gene[1], condition = s$condition[1], bin = s$bin[1], density = mean(s$density),
                 mean = mean(s$expr), half = if (n > 1L) stats::qt(0.975, n - 1L) * se else NA_real_)
    }))
    sm <- sm[sm$mean > 0, , drop = FALSE]
    if (!nrow(sm)) stop("no counts of the chosen genes in these cells", call. = FALSE)
    floor_y <- min(sm$mean) / 2
    sm$lo <- pmax(sm$mean - sm$half, floor_y)
    sm$hi <- sm$mean + sm$half
    sm$gene_lab <- factor(lab[sm$gene], unname(lab))
    ln <- do.call(rbind, lapply(split(sm, list(sm$gene, sm$condition), drop = TRUE), function(s) {
      g <- s$gene[1]
      cv <- as.character(s$condition[1])
      b <- sl$group(g, names(grp)[grp == cv])
      x0 <- mean(s$density)
      y0 <- mean(log(s$mean))
      xs <- range(s$density)
      data.frame(gene_lab = s$gene_lab[1], condition = s$condition[1], x = xs[1], xend = xs[2],
                 y = exp(y0 + b * (xs[1] - x0)), yend = exp(y0 + b * (xs[2] - x0)))
    }))
    p <- p +
      geom_point(data = function(d) d[d$expr > 0, ], aes(colour = .data$condition), size = 0.6, alpha = 0.3,
                 shape = 16) +
      geom_segment(data = ln, aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend,
                                  colour = .data$condition, linetype = "fitted slope"),
                   linewidth = 0.8, inherit.aes = FALSE) +
      geom_line(data = sm, aes(.data$density, .data$mean, colour = .data$condition,
                               linetype = "mean across patients"), linewidth = 0.5) +
      geom_pointrange(data = sm, aes(.data$density, .data$mean, ymin = .data$lo, ymax = .data$hi,
                                     colour = .data$condition), size = 0.3, linewidth = 0.5,
                      inherit.aes = FALSE, na.rm = TRUE) +
      scale_linetype_manual(values = c(`mean across patients` = "solid", `fitted slope` = "22"), name = NULL)
    ylab <- "expression relative to patient (CP10k)"
  } else {
    cen <- do.call(rbind, lapply(gene, function(g) {
      e <- as.numeric(cells$counts[g, ]) / cells$lib * 1e4
      d <- stats::aggregate(cbind(expr = e, density = l) ~ patient, data.frame(e = e, l = l, patient = pat), mean)
      d$gene <- g
      d$beta <- sl$patient(g, d$patient)
      rng <- t(vapply(split(l, pat)[d$patient], stats::quantile, numeric(2), probs = c(0.05, 0.95)))
      d$x <- rng[, 1]
      d$xend <- rng[, 2]
      d$y <- exp(log(d$expr) + d$beta * (d$x - d$density))
      d$yend <- exp(log(d$expr) + d$beta * (d$xend - d$density))
      d
    }))
    cen <- cen[is.finite(cen$beta) & cen$expr > 0, , drop = FALSE]
    cen$condition <- grp[cen$patient]
    cen$gene_lab <- factor(lab[cen$gene], unname(lab))
    p <- p +
      geom_point(data = function(d) d[d$expr > 0, ], aes(colour = .data$condition), size = 0.6, alpha = 0.35,
                 shape = 16) +
      geom_segment(data = cen, aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend,
                                   colour = .data$condition, linetype = "within patient"),
                   linewidth = 0.35, inherit.aes = FALSE) +
      geom_smooth(data = cen, aes(.data$density, .data$expr, linetype = "between patients"), method = "lm",
                  formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.7, inherit.aes = FALSE) +
      geom_point(data = cen, aes(.data$density, .data$expr, colour = .data$condition), shape = 21, fill = "white",
                 size = 1.8, stroke = 0.6, inherit.aes = FALSE) +
      scale_linetype_manual(values = c(`within patient` = "solid", `between patients` = "22"), name = NULL)
    ylab <- "expression (CP10k)"
  }
  p + .scaleCondition(lv) +
    scale_y_log10(labels = scales::label_number(drop0trailing = TRUE, big.mark = ",")) +
    facet_wrap(~ gene_lab, scales = "free_y", labeller = label_parsed) +
    labs(x = sprintf("%s density (log1p)", niche), y = ylab, colour = NULL) +
    guides(colour = guide_legend(order = 1, override.aes = list(alpha = 1, linetype = 0,
                                                                size = if (adjust) 0.6 else 1.8)),
           linetype = guide_legend(order = 2, override.aes = list(colour = "black"))) +
    theme_spiDE()
}

# The slopes drawn by plotNicheResponse(): each patient's own slope and each
# group's weighted mean (the slopes engine's pooled-test weights), or, under
# the sandwich engine, the pooled slope for every patient and group.
.responseSlopes <- function(fit, xi, gene, niche) {
  if (fit@engine == "slopes") {
    j <- match(niche, xi$niches)
    sw <- .slopeWeights(xi)
    get <- function(g) {
      gi <- match(g, xi$genes)
      list(b = stats::setNames(xi$beta[gi, , j], xi$patients),
           w = stats::setNames(1 / (sw$v[gi, , j] + sw$tau[gi, j]), xi$patients))
    }
    return(list(
      patient = function(g, pats) get(g)$b[pats],
      group = function(g, pats) {
        o <- get(g)
        k <- intersect(pats, names(o$b))
        k <- k[is.finite(o$b[k]) & is.finite(o$w[k]) & o$w[k] > 0]
        if (!length(k)) NA_real_ else sum(o$w[k] * o$b[k]) / sum(o$w[k])
      }))
  }
  cf <- xi$coef[xi$coef$test == "pooled" & xi$coef$niche == niche, , drop = FALSE]
  pooled <- function(g) cf$estimate[match(g, cf$gene)]
  list(patient = function(g, pats) rep(pooled(g), length(pats)),
       group = function(g, pats) pooled(g))
}
