# The patient-level views (0.99.34): each patient's slopes for a triplet, and
# what the patients' intercepts (or slopes) capture. Design:
# design/specs/2026-10-01-plots.md.

#' Each patient's niche slope for a triplet
#'
#' The evidence the tests weigh: every patient's own slope of the gene's
#' expression in the index type on the niche type's density (log1p), with an
#' interval from its within-patient variance -- the variance the tests use --
#' and a point sized by the patient's weight in them. The pooled estimate and
#' its confidence interval are the grey band. With a condition the patients
#' are split by it, each group's weighted mean is a coloured line, and each
#' strip gives the condition-specific slope difference and its q-value.
#'
#' Read it for whether a call rests on patients that agree or on one or two
#' heavily weighted ones, and for why a condition test is missing (each group
#' needs three effective patients). Patients are ordered by group, then by
#' their mean slope over the genes shown, so a row is one patient across
#' genes.
#'
#' @param x a [SpiDEResults-class] from the slopes engine.
#' @param gene one or more genes, one column each.
#' @param index,niche the index and niche cell types.
#' @param ci the confidence level of the intervals.
#' @return a ggplot. \code{p$data} has one row per (gene, patient) with a
#'   slope: \code{gene}, \code{patient}, \code{slope}, \code{se} (within
#'   patient), \code{weight} (the patient's share of the pooled test's
#'   weights), \code{condition} and \code{gene_lab} (the strip label, plotmath).
#'   Change it like any ggplot, e.g. \code{+ theme(axis.text.y = element_blank())}.
#' @examples
#' data(toySpiDE)
#' res <- spiDE(toySpiDE, condition = "condition", sigma = 30, index = "A", verbose = FALSE)
#' plotPatientSlopes(res, gene = c("G1", "G2"), index = "A", niche = "B")
#' @seealso [patientSlopes()] for the numbers, [plotNicheResponse()] for the cells.
#' @export
plotPatientSlopes <- function(x, gene, index, niche, ci = 0.95) {
  checkResults(x, "plotPatientSlopes()")
  fit <- x@fit
  checkSlopesEngine(fit, "plotPatientSlopes()")
  gene <- unique(gene)
  checkTriplet(fit, index, niche, gene)
  if (!is.numeric(ci) || length(ci) != 1L || !is.finite(ci) || ci <= 0 || ci >= 1) {
    stop("'ci' must be a single number in (0, 1)", call. = FALSE)
  }
  xi <- fit@index[[index]]
  j <- match(niche, xi$niches)
  sw <- .slopeWeights(xi)
  z <- stats::qnorm(1 - (1 - ci) / 2)
  grp <- .groupLevels(fit, x@condition)
  d <- do.call(rbind, lapply(gene, function(g) {
    gi <- match(g, xi$genes)
    b <- xi$beta[gi, , j]
    v <- sw$v[gi, , j]
    w <- 1 / (v + sw$tau[gi, j])
    ok <- is.finite(b) & is.finite(w) & w > 0
    data.frame(gene = g, patient = xi$patients, slope = b, se = sqrt(v),
               weight = ifelse(ok, w / sum(w[ok]), NA_real_), stringsAsFactors = FALSE)[ok, , drop = FALSE]
  }))
  if (!nrow(d)) stop("no patient has a usable slope for these genes", call. = FALSE)
  d$condition <- grp[d$patient]
  ord <- stats::aggregate(slope ~ patient + condition, d, mean)
  ord <- ord[order(ord$condition, ord$slope), ]
  d$patient <- factor(d$patient, unique(ord$patient))
  lab <- .slopeStripLabels(x, gene, index, niche)
  d$gene_lab <- factor(lab[d$gene], unname(lab))
  pr <- .tripletRows(x@table, "pooled", gene, index, niche)
  tq <- stats::qt(1 - (1 - ci) / 2, pr$df)
  band <- data.frame(gene_lab = factor(unname(lab), unname(lab)), est = pr$estimate,
                     lo = pr$estimate - tq * pr$se, hi = pr$estimate + tq * pr$se)
  gm <- do.call(rbind, lapply(split(d, list(d$gene_lab, d$condition), drop = TRUE), function(s) {
    data.frame(gene_lab = s$gene_lab[1], condition = s$condition[1],
               mean = sum(s$weight * s$slope) / sum(s$weight))
  }))
  p <- ggplot(d, aes(.data$slope, .data$patient)) +
    geom_rect(data = band, aes(xmin = .data$lo, xmax = .data$hi, ymin = -Inf, ymax = Inf,
                               fill = sprintf("pooled estimate (%g%% CI)", 100 * ci)),
              inherit.aes = FALSE) +
    geom_vline(data = band, aes(xintercept = .data$est), colour = .spideCols[["grey"]], linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3) +
    geom_vline(data = gm, aes(xintercept = .data$mean, colour = .data$condition), linewidth = 0.8) +
    geom_linerange(aes(xmin = .data$slope - z * .data$se, xmax = .data$slope + z * .data$se,
                       colour = .data$condition), linewidth = 0.4) +
    geom_point(aes(colour = .data$condition, size = .data$weight)) +
    .scaleCondition(levels(grp)) +
    scale_fill_manual(values = .spideCols[["faint"]], name = NULL) +
    scale_size_area(max_size = 2.8, name = "weight", labels = scales::label_percent()) +
    facet_grid(condition ~ gene_lab, scales = "free", space = "free_y",
               labeller = labeller(gene_lab = label_parsed)) +
    labs(x = sprintf("slope on %s density (log1p)", niche), y = NULL) +
    guides(colour = "none", size = guide_legend(order = 1, override.aes = list(colour = .spideCols[["grey"]])),
           fill = guide_legend(order = 2)) +
    theme_spiDE() + theme(axis.text.y = element_text(size = rel(0.65)))
  if (length(x@condition)) p <- p + labs(caption = sprintf("\u0394: slope difference, %s", x@contrast))
  p
}

# Each gene's strip: the italic symbol over the condition test's slope
# difference and q (or, without a condition, the pooled slope and q).
.slopeStripLabels <- function(x, gene, index, niche) {
  cond <- length(x@condition) == 1L
  r <- .tripletRows(x@table, if (cond) "condition" else "pooled", gene, index, niche)
  est <- ifelse(is.finite(r$estimate), formatC(r$estimate, format = "f", digits = 2, flag = "+"), "NA")
  stats::setNames(sprintf('atop(%s, scriptstyle("%s = %s, q = %s"))', .itGene(gene),
                          if (cond) "\u0394" else "slope", est, .fmtq(r$q)), gene)
}

#' What the patient intercepts, or slopes, capture
#'
#' Every index type's shared fit gives each patient an intercept per gene
#' (see [patientIntercepts()]): they absorb every difference between patients
#' -- the condition's main effect, composition, batch -- which is what makes
#' the niche slopes within-patient slopes. This relates the principal
#' components of the patients' gene-centred intercepts to patient covariates:
#' every patient-level colData column, plus the patients' index cells (log10),
#' mean log library size, and mean log1p density of each tested niche type.
#' Each tile is the adjusted \eqn{R^2} of a component on a covariate.
#'
#' Structure in the intercepts is harmless to the tests: it is absorbed by
#' design. With \code{niche}, the same analysis runs on the patients' slopes on
#' that niche (each gene's slopes scaled by its between-patient spread, as in
#' [spiGSEA()]): a slide or batch that organises the slopes can reach the
#' condition test, and belongs in \code{strata}.
#'
#' @param x a [SpiDEResults-class] or [SpiDEFit-class] (spiDE >= 0.99.34).
#' @param index the index cell type.
#' @param niche \code{NULL} (the intercepts) or a niche type (its slopes;
#'   slopes engine).
#' @param covariates \code{NULL} (all) or the covariates to show, by the names
#'   in the plot.
#' @param type \code{"association"} (the \eqn{R^2} tiles) or \code{"pca"}
#'   (the first two components, coloured by \code{colour.by}).
#' @param colour.by for \code{type = "pca"}: a covariate (default the
#'   condition, else the first).
#' @param n.pcs the number of components shown.
#' @param ntop the number of genes, most expressed first, whose intercepts
#'   enter the components (a gene with any non-finite value is left out).
#' @return a ggplot. For \code{"association"}, \code{p$data} has \code{pc},
#'   \code{covariate}, \code{r2} and \code{group} (patient or derived); for
#'   \code{"pca"}, \code{patient}, \code{PC1}, \code{PC2} and the covariate.
#' @examples
#' data(toySpiDE)
#' res <- spiDE(toySpiDE, condition = "condition", sigma = 30, index = "A", verbose = FALSE)
#' plotPatientEffects(res, index = "A")
#' plotPatientEffects(res, index = "A", niche = "B")
#' @seealso [patientIntercepts()], [testNicheAbundance()] for the composition
#'   association gene by gene.
#' @export
plotPatientEffects <- function(x, index, niche = NULL, covariates = NULL, type = c("association", "pca"),
                               colour.by = NULL, n.pcs = 5L, ntop = 500L) {
  checkResults(x, "plotPatientEffects()", fit.ok = TRUE)
  type <- match.arg(type)
  fit <- if (is(x, "SpiDEResults")) x@fit else x
  condition <- if (is(x, "SpiDEResults")) x@condition else fit@condition
  checkTriplet(fit, index, niche)
  xi <- fit@index[[index]]
  checkIntercepts(xi, "plotPatientEffects()")
  if (!is.null(niche)) checkSlopesEngine(fit, "plotPatientEffects(niche = )")
  M <- .patientEffectMatrix(xi, niche, ntop)
  if (nrow(M) < 2L || ncol(M) < 3L) stop("too few genes or patients for principal components", call. = FALSE)
  pc <- stats::prcomp(t(M), center = TRUE, scale. = FALSE)
  ve <- pc$sdev^2 / sum(pc$sdev^2)
  n.pcs <- max(1L, min(as.integer(n.pcs), ncol(pc$x)))
  cv <- .patientCovariates(fit, xi, covariates)
  cvm <- cv$cov[rownames(pc$x), , drop = FALSE]
  what <- if (is.null(niche)) "the patient intercepts" else sprintf("the patient slopes on %s", niche)
  if (type == "pca") {
    if (is.null(colour.by)) colour.by <- if (length(condition) && condition %in% names(cvm)) condition else names(cvm)[1]
    if (!colour.by %in% names(cvm)) {
      stop(sprintf("'colour.by' must be one of: %s", paste(names(cvm), collapse = ", ")), call. = FALSE)
    }
    d <- data.frame(patient = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2], stringsAsFactors = FALSE)
    d[[colour.by]] <- cvm[[colour.by]]
    v <- d[[colour.by]]
    sc <- if (length(condition) && identical(colour.by, condition)) {
      .scaleCondition(attr(.conditionCoding(fit@patients, condition), "levels"))
    } else if (is.numeric(v)) {
      scale_colour_gradientn(colours = .magnitudeRamp(256L, skip = 0.1))
    } else {
      scale_colour_manual(values = spiDEColours("nuisance", nlevels(factor(v))))
    }
    return(ggplot(d, aes(.data$PC1, .data$PC2, colour = .data[[colour.by]])) +
             geom_point(size = 2.2) + sc +
             labs(x = sprintf("PC1 (%.0f%%)", 100 * ve[1]), y = sprintf("PC2 (%.0f%%)", 100 * ve[2])) +
             theme_spiDE())
  }
  a <- expand.grid(pc = seq_len(n.pcs), covariate = names(cvm), stringsAsFactors = FALSE)
  a$r2 <- mapply(function(k, v) .adjR2(pc$x[, k], cvm[[v]]), a$pc, a$covariate)
  pcl <- sprintf("PC%d\n%.0f%%", seq_len(n.pcs), 100 * ve[seq_len(n.pcs)])
  a$pc <- factor(pcl[a$pc], pcl)
  a$group <- factor(cv$group[a$covariate], c("patient", "derived"))
  a$covariate <- factor(a$covariate, rev(names(cvm)))
  ggplot(a, aes(.data$pc, .data$covariate, fill = .data$r2)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = ifelse(is.finite(.data$r2) & .data$r2 >= 0.05, sprintf("%.2f", .data$r2), ""),
                  colour = is.finite(.data$r2) & .data$r2 > 0.45), size = 3.2, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
    scale_fill_gradientn(colours = .magnitudeRamp(256L), limits = c(0, 1), na.value = .spideCols[["faint"]],
                         name = expression(adjusted ~ R^2)) +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    labs(x = sprintf("PCs of %s", what), y = NULL) +
    theme_spiDE()
}

# Genes x patients, gene-centred: the intercepts, or the slopes on one niche
# scaled by each gene's between-patient spread (as spiGSEA() scales them). A
# patient missing for more than a tenth of the genes drops out, then any gene
# not finite in every remaining patient; the ntop most expressed genes enter.
.patientEffectMatrix <- function(xi, niche, ntop) {
  if (is.null(niche)) {
    M <- xi$intercept
  } else {
    j <- match(niche, xi$niches)
    sw <- .slopeWeights(xi)
    b <- matrix(xi$beta[, , j], length(xi$genes))
    v <- matrix(sw$v[, , j], length(xi$genes))
    sig <- sqrt(sw$tau[, j] + apply(v, 1, stats::median, na.rm = TRUE))
    M <- b / sig
    dimnames(M) <- list(xi$genes, xi$patients)
  }
  M <- M[, colMeans(is.finite(M)) >= 0.9, drop = FALSE]
  ok <- rownames(M)[rowSums(!is.finite(M)) == 0]
  keep <- utils::head(intersect(names(sort(xi$mean_expr, decreasing = TRUE)), ok), ntop)
  M <- M[keep, , drop = FALSE]
  M - rowMeans(M)
}

# The covariates of plotPatientEffects(): every patient-level colData column
# that varies and is not an identifier, plus the derived ones.
.patientCovariates <- function(fit, xi, covariates = NULL) {
  pt <- fit@patients[match(xi$patients, fit@patients$patient), , drop = FALSE]
  cand <- setdiff(colnames(pt), c("patient", "ncells"))
  keep <- vapply(cand, function(cn) {
    v <- pt[[cn]]
    u <- length(unique(v[!is.na(v)]))
    u > 1L && (is.numeric(v) || u < nrow(pt))
  }, logical(1))
  pc <- pt[, cand[keep], drop = FALSE]
  pc[] <- lapply(pc, function(v) if (is.numeric(v)) v else factor(v))
  der <- data.frame(`cells (log10)` = log10(as.numeric(xi$ncells[xi$patients])),
                    `depth (mean log)` = as.numeric(xi$loglib_mean[xi$patients]), check.names = FALSE)
  nm <- as.data.frame(xi$niche_mean[xi$patients, , drop = FALSE])
  colnames(nm) <- paste(colnames(nm), "density")
  cv <- cbind(pc, der, nm)
  rownames(cv) <- xi$patients
  group <- stats::setNames(c(rep("patient", ncol(pc)), rep("derived", ncol(der) + ncol(nm))), colnames(cv))
  if (!is.null(covariates)) {
    miss <- setdiff(covariates, colnames(cv))
    if (length(miss)) {
      stop(sprintf("covariate(s) not available: %s (available: %s)", paste(miss, collapse = ", "),
                   paste(colnames(cv), collapse = ", ")), call. = FALSE)
    }
    cv <- cv[, covariates, drop = FALSE]
    group <- group[covariates]
  }
  list(cov = cv, group = group)
}

# Adjusted R^2 of y on one covariate, floored at 0; NA when it cannot be fitted.
.adjR2 <- function(y, v) {
  ok <- is.finite(y) & !is.na(v)
  y <- y[ok]
  v <- if (is.factor(v)) droplevels(v[ok]) else v[ok]
  if (length(unique(v)) < 2L || length(y) < 4L) return(NA_real_)
  s <- tryCatch(suppressWarnings(summary(stats::lm(y ~ v))$adj.r.squared), error = function(e) NA_real_)
  if (!is.finite(s)) NA_real_ else max(0, s)
}
