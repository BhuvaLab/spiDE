# The results views (0.99.34): where the calls are, the strongest triplets,
# the spillover signature, triplet heatmaps and p-value bands. Design:
# design/specs/2026-10-01-plots.md.

#' Where the calls are
#'
#' The number of genes called for every (index, niche) pair, split by
#' direction: for the pooled test, a slope that rises or falls with the niche's
#' density; for the condition-specific test, a slope steeper in one condition
#' or the other. The heatmap gives exact counts and marks untested pairs (an
#' index type against its own niche) in grey; the graph draws cell types on a
#' fixed circle with an arrow from niche to index type, its width the number
#' of genes -- it reads best when calls are few (the condition test), the
#' heatmap when every pair has some.
#'
#' Read it for which cell-type pairs carry the signal, and for pairs whose
#' calls all go one way toward one niche type (the shape segmentation spillover
#' gives; see [plotSpillover()]). An arrow is an association, not signalling.
#'
#' @param x a [SpiDEResults-class].
#' @param test \code{NULL} (the condition-specific test if run, else the
#'   pooled test), \code{"pooled"} or \code{"condition"}.
#' @param fdr the q-value threshold of a call.
#' @param style \code{"heatmap"} or \code{"graph"}.
#' @param min.calls the graph's smallest drawn edge.
#' @return a ggplot. For the heatmap, \code{p$data} has one row per (index,
#'   niche, direction): \code{calls} and \code{tested}.
#' @examples
#' data(toySpiDE)
#' res <- spiDE(toySpiDE, condition = "condition", sigma = 30, verbose = FALSE)
#' plotCallMap(res, test = "pooled")
#' plotCallMap(res, style = "graph")
#' @export
plotCallMap <- function(x, test = NULL, fdr = x@fdr, style = c("heatmap", "graph"), min.calls = 1L) {
  checkResults(x, "plotCallMap()")
  style <- match.arg(style)
  test <- .resolveTest(x, test)
  checkFdr(fdr)
  tab <- x@table[x@table$test == test, , drop = FALSE]
  tested <- unique(tab[is.finite(tab$p), c("index", "niche")])
  hit <- tab[!is.na(tab$q) & tab$q <= fdr, , drop = FALSE]
  n <- if (nrow(hit)) {
    stats::aggregate(list(calls = hit$gene), list(index = hit$index, niche = hit$niche,
                                                  direction = ifelse(hit$estimate > 0, "up", "down")), length)
  } else data.frame(index = character(), niche = character(), direction = character(), calls = integer())
  lab <- .directionLabels(x, test)
  if (style == "heatmap") .callHeatmap(n, tested, lab) else .callGraph(n, tested, lab, min.calls)
}

.callHeatmap <- function(n, tested, lab) {
  grid <- expand.grid(index = sort(unique(tested$index)), niche = sort(unique(tested$niche)),
                      direction = c("up", "down"), stringsAsFactors = FALSE)
  grid$tested <- paste(grid$index, grid$niche) %in% paste(tested$index, tested$niche)
  grid <- merge(grid, n, all.x = TRUE)
  grid$calls[is.na(grid$calls) & grid$tested] <- 0L
  grid$signed <- ifelse(grid$direction == "up", 1, -1) * log10(1 + grid$calls)
  grid$direction <- factor(grid$direction, c("up", "down"), unname(lab[c("up", "down")]))
  lim <- max(1, abs(grid$signed), na.rm = TRUE)
  ggplot(grid, aes(.data$niche, .data$index)) +
    geom_tile(data = function(d) d[!d$tested, ], fill = .spideCols[["faint"]], colour = "white", linewidth = 0.6) +
    geom_tile(data = function(d) d[d$tested, ], aes(fill = .data$signed), colour = "white", linewidth = 0.6) +
    geom_text(data = function(d) d[d$tested, ], aes(label = .data$calls, colour = abs(.data$signed) > 0.55 * lim),
              size = 3.4, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
    scale_fill_gradient2(low = .spideCols[["down"]], mid = .spideCols[["mid"]], high = .spideCols[["up"]],
                         midpoint = 0, limits = c(-lim, lim), guide = "none") +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    facet_wrap(~ direction) + coord_equal() +
    labs(x = "niche type", y = "index type") +
    theme_spiDE() + theme(axis.text.x = element_text(angle = 40, hjust = 1))
}

.callGraph <- function(n, tested, lab, min.calls) {
  types <- sort(unique(c(tested$index, tested$niche)))
  ang <- pi / 2 - 2 * pi * (seq_along(types) - 1) / length(types)
  pos <- data.frame(type = types, x = cos(ang), y = sin(ang), stringsAsFactors = FALSE)
  e <- n[n$calls >= min.calls, , drop = FALSE]
  e$x0 <- pos$x[match(e$niche, pos$type)]
  e$y0 <- pos$y[match(e$niche, pos$type)]
  e$x1 <- pos$x[match(e$index, pos$type)]
  e$y1 <- pos$y[match(e$index, pos$type)]
  r <- 0.13
  len <- sqrt((e$x1 - e$x0)^2 + (e$y1 - e$y0)^2)
  e$xs <- e$x0 + r * (e$x1 - e$x0) / len
  e$ys <- e$y0 + r * (e$y1 - e$y0) / len
  e$xe <- e$x1 - r * (e$x1 - e$x0) / len
  e$ye <- e$y1 - r * (e$y1 - e$y0) / len
  e$panel <- factor(unname(lab[e$direction]), unname(lab[c("up", "down")]))
  nodes <- merge(pos, data.frame(panel = factor(unname(lab[c("up", "down")]), unname(lab[c("up", "down")]))))
  # labels outside the nodes, anchored away from the circle's centre
  nodes$hj <- ifelse(abs(nodes$x) < 0.2, 0.5, ifelse(nodes$x > 0, 0, 1))
  nodes$vj <- ifelse(nodes$y > 0.3, 0, ifelse(nodes$y < -0.3, 1, 0.5))
  ggplot() +
    geom_curve(data = e, aes(x = .data$xs, y = .data$ys, xend = .data$xe, yend = .data$ye,
                             linewidth = .data$calls, colour = .data$direction),
               curvature = 0.22, lineend = "round",
               arrow = arrow(length = unit(2, "mm"), type = "closed")) +
    geom_point(data = nodes, aes(.data$x, .data$y), size = 6, shape = 21, fill = "white", colour = "black",
               stroke = 0.6) +
    geom_text(data = nodes, aes(1.3 * .data$x, 1.3 * .data$y, label = .data$type, hjust = .data$hj,
                                vjust = .data$vj), size = 3.6) +
    .scaleDirection(lab) +
    scale_linewidth(range = c(0.3, 2.6), name = "genes") +
    guides(colour = "none", linewidth = guide_legend(override.aes = list(colour = .spideCols[["grey"]]))) +
    facet_wrap(~ panel, drop = FALSE) + coord_equal(clip = "off") +
    scale_x_continuous(expand = expansion(add = 0.6)) + scale_y_continuous(expand = expansion(add = 0.3)) +
    theme_spiDE() + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
}

#' The strongest triplets, or those of chosen genes
#'
#' A forest plot: each triplet's estimate (the pooled slope, or the slope
#' difference between conditions) with its 95\% confidence interval, faded
#' where it is not called. \code{gene}, \code{index} and \code{niche} filter
#' the table before the \code{n} smallest p-values are taken, so
#' \code{gene = c("CD74", "HLA-DRA")} shows those genes' triplets. With
#' \code{spe}, a hollow triangle marks a gene at least \code{fold}-fold higher
#' in the niche type than in the index type (see [plotSpillover()]).
#'
#' @param x a [SpiDEResults-class].
#' @param test \code{NULL} (the condition-specific test if run, else the
#'   pooled test), \code{"pooled"} or \code{"condition"}.
#' @param n the number of triplets.
#' @param gene,index,niche optional filters.
#' @param spe \code{NULL} or the SpatialExperiment, for the niche-marker flag.
#' @param fold the niche-marker threshold.
#' @param fdr the q-value threshold of a call.
#' @param assay,name the counts assay and the niche reducedDim prefix.
#' @return a ggplot; \code{p$data} is the plotted rows of the results table
#'   with \code{lo}, \code{hi}, \code{called}, \code{direction},
#'   \code{marker} and \code{label} (plotmath).
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = 30)
#' res <- spiDE(spe, condition = "condition", sigma = 30, verbose = FALSE)
#' plotTopTriplets(res, test = "pooled", n = 15, spe = spe)
#' plotTopTriplets(res, gene = "G1")
#' @export
plotTopTriplets <- function(x, test = NULL, n = 30L, gene = NULL, index = NULL, niche = NULL, spe = NULL,
                            fold = 4, fdr = x@fdr, assay = "counts", name = "Niche") {
  checkResults(x, "plotTopTriplets()")
  test <- .resolveTest(x, test)
  checkFdr(fdr)
  tab <- results(x, test = test)
  if (!is.null(gene)) tab <- tab[tab$gene %in% gene, , drop = FALSE]
  if (!is.null(index)) tab <- tab[tab$index %in% index, , drop = FALSE]
  if (!is.null(niche)) tab <- tab[tab$niche %in% niche, , drop = FALSE]
  tab <- tab[is.finite(tab$p) & is.finite(tab$se), , drop = FALSE]
  if (!nrow(tab)) stop("no triplet of the chosen genes, index and niche types was tested", call. = FALSE)
  tab <- utils::head(tab[order(tab$p), , drop = FALSE], n)
  tq <- stats::qt(0.975, tab$df)
  tab$lo <- tab$estimate - tq * tab$se
  tab$hi <- tab$estimate + tq * tab$se
  tab$called <- factor(ifelse(!is.na(tab$q) & tab$q <= fdr, "yes", "no"), c("yes", "no"))
  tab$direction <- factor(ifelse(tab$estimate > 0, "up", "down"), c("up", "down"))
  tab$marker <- if (is.null(spe)) rep(FALSE, nrow(tab)) else .nicheFold(tab, spe, x@fit, assay, name) >= fold
  tab$label <- .itTriplet(tab$gene, tab$index, tab$niche)
  tab$label <- factor(tab$label, tab$label[order(tab$estimate)])
  ggplot(tab, aes(.data$estimate, .data$label, colour = .data$direction, alpha = .data$called)) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3) +
    geom_linerange(aes(xmin = .data$lo, xmax = .data$hi), linewidth = 0.6) +
    geom_point(aes(shape = .data$marker), size = 2, fill = "white", stroke = 0.8) +
    .scaleDirection(.directionLabels(x, test)) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 24), breaks = "TRUE",
                       labels = sprintf("≥ %g-fold higher\nin the niche type", fold), name = NULL,
                       drop = FALSE) +
    scale_alpha_manual(values = c(yes = 1, no = 0.35), name = sprintf("q ≤ %s", format(fdr)), drop = FALSE) +
    scale_y_discrete(labels = .parseLabels) +
    labs(x = if (test == "pooled") "slope on niche density (log1p)" else sprintf("slope difference (%s)", x@contrast),
         y = NULL, colour = NULL) +
    guides(colour = guide_legend(order = 1, ncol = 1),
           alpha = guide_legend(order = 2, ncol = 1, override.aes = list(colour = "black", shape = 16)),
           shape = guide_legend(order = 3)) +
    theme_spiDE() + theme(legend.position = "bottom", legend.justification = "left")
}

#' Are the calls the niche type's own markers?
#'
#' Every tested triplet, placed by how much higher the gene is in the niche
#' type than in the index type (log2, from each cell type's mean counts per
#' 10,000) against the signed z of its test; calls are coloured by direction.
#' A gene several-fold higher in the niche type whose measured level in the
#' index cells rises with that niche's density is what transcripts leaking
#' across segmentation boundaries would produce: the association is real in
#' the data, but it needs orthogonal evidence (nuclear counts, a
#' spillover-corrected segmentation) before it is read as biology.
#'
#' @param x a [SpiDEResults-class].
#' @param spe the SpatialExperiment the results came from.
#' @param test \code{"pooled"} or \code{"condition"}.
#' @param fold the niche-marker threshold (shaded).
#' @param label the number of markers labelled (strongest first).
#' @param fdr the q-value threshold of a call.
#' @param assay,name the counts assay and the niche reducedDim prefix.
#' @return a ggplot; \code{p$data} is the test's rows of the results table
#'   with \code{fold} (log2), \code{z} and \code{status}.
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = 30)
#' res <- spiDE(spe, condition = "condition", sigma = 30, verbose = FALSE)
#' plotSpillover(res, spe)
#' @export
plotSpillover <- function(x, spe, test = "pooled", fold = 4, label = 6L, fdr = x@fdr, assay = "counts",
                          name = "Niche") {
  checkResults(x, "plotSpillover()")
  test <- .resolveTest(x, test)
  checkFdr(fdr)
  checkSPE(spe, assay = assay, cell_type = x@fit@params$cell_type, sample_id = x@fit@params$sample_id)
  tab <- x@table[x@table$test == test & is.finite(x@table$p), , drop = FALSE]
  tab$fold <- log2(.nicheFold(tab, spe, x@fit, assay, name))
  tab$z <- sign(tab$t) * stats::qnorm(pmin(pmax(tab$p, 1e-300), 1) / 2, lower.tail = FALSE)
  called <- !is.na(tab$q) & tab$q <= fdr
  tab$status <- factor(ifelse(called, ifelse(tab$estimate > 0, "up", "down"), "not called"),
                       c("up", "down", "not called"))
  mk <- called & is.finite(tab$fold) & tab$fold >= log2(fold)
  lab <- tab[mk, , drop = FALSE]
  lab <- utils::head(lab[order(-abs(lab$z)), , drop = FALSE], label)
  lab$text <- .itTriplet(lab$gene, lab$index, lab$niche)
  labels <- if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(data = lab, aes(label = .data$text), parse = TRUE, size = 3, min.segment.length = 0,
                             segment.size = 0.3, box.padding = 0.35, max.overlaps = Inf, seed = 1,
                             xlim = c(-Inf, log2(fold) - 0.4), direction = "y", hjust = 1)
  } else {
    geom_text(data = lab, aes(label = .data$text), parse = TRUE, size = 3, hjust = 1.05, check_overlap = TRUE)
  }
  ggplot(tab, aes(.data$fold, .data$z)) +
    annotate("rect", xmin = log2(fold), xmax = Inf, ymin = -Inf, ymax = Inf, fill = .spideCols[["faint"]]) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5, size = 3.4,
             label = sprintf("%d of %d calls\n≥ %g-fold in niche type", sum(mk), sum(called), fold)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = log2(fold), linetype = "dashed", linewidth = 0.4) +
    geom_point(data = function(d) d[d$status == "not called", ], aes(colour = .data$status), size = 0.7, shape = 16) +
    geom_point(data = function(d) d[d$status != "not called", ], aes(colour = .data$status), size = 1.1, shape = 16) +
    labels +
    .scaleDirection(.directionLabels(x, test), extra = c(`not called` = .spideCols[["light"]])) +
    labs(x = "log2 fold (niche type / index type)", y = "signed z", colour = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 2.5))) +
    theme_spiDE()
}
