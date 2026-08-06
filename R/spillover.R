# Spillover diagnostic for niche associations.
#
# Transcript misassignment ("spillover") makes a cell's counts include a
# fraction of its neighbours'. To first order the contamination a cell of index
# type c receives for gene g is
#
#     kappa * sum_n d_n(cell) * m_ng
#
# with d_n the local density of type n and m_ng type n's mean for gene g. That
# is the same shape as the CellType_c:niche_n covariate, so under a pure null
# the induced coefficient is
#
#     beta_cn,g  proportional to  kappa * m_ng / M_cg
#
# -- the RATIO of the niche type's expression to the index type's, per gene.
#
# Why this matters differently in the two modes. In condition mode the two-way
# CellType:niche term is a free nuisance parameter and the tested effect is the
# three-way condition interaction; since leakage is condition-independent, the
# contamination loads onto the nuisance term and cancels out of the test. In
# NICHE mode the two-way term IS the tested term, so contamination is
# confounded with the biology at the point of estimation.
#
# Simulation backing these numbers (kappa = 0.3, 6 samples x 150 cells,
# 40 genes, leak radius 20 um), Spearman-equivalent Pearson correlations of the
# tested coefficient with the log marker ratio:
#
#                                     kappa = 0     kappa = 0.3
#   niche mode, tested (two-way)        -0.05          +0.82
#   condition mode, tested (three-way)  +0.05          +0.03
#   condition mode, two-way nuisance       -            +0.64
#
# So the signature is strong, specific to spillover, and confined to niche mode
# exactly as the algebra predicts.
#
# Two candidate CORRECTIONS were measured against this artefact, both of which
# work by estimating a niche effect common to all index cell types and removing
# it. Scored by the same signature (median |Spearman(coef, log marker ratio)|;
# the kappa = 0 noise floor is ~0.14) and by signal-to-noise against null genes
# for two kinds of planted biology -- one specific to a single index cell type,
# one shared by all of them:
#
#   D1  cell-means beta_cn                       (what spiDE tests today)
#   D2  deviation  delta_cn = beta_cn - mean_c'  (a contrast INSIDE D1's span;
#                                                 the common term is estimated
#                                                 from index cells only)
#   D3  additive   D1 + a bare niche main effect (a strictly LARGER model: the
#                                                 main effect has support on
#                                                 type-n cells, which D1 does
#                                                 not, so it reintroduces a
#                                                 self-niche slope)
#
#                              D1      D2      D3
#   signature, kappa=0.3      0.617   0.280   0.448     (floor ~0.14)
#   A-specific SNR, kappa=0    7.06    2.87    7.93
#   shared SNR,     kappa=0    5.74    0.19    8.02
#
# So neither dominates. D2 removes ~71% of the spillover signature above the
# floor but destroys cell-type-SHARED niche biology entirely (5.74 -> 0.19 with
# no spillover present at all) and attenuates cell-type-specific biology by the
# factor (k-1)/k, k being the number of index types for that niche -- worst at
# k = 2. D3 is nearly free statistically but removes only ~35% of the
# signature, because its common term is estimated using type-n cells too, and
# those barely change under same-type contamination.
#
# The destruction of shared biology is not a tuning problem: a niche response
# common to every index cell type is mathematically indistinguishable from
# spillover under the assumption that makes either correction work, namely that
# leakage affects all index cell types alike.
#
# This is deliberately a DIAGNOSTIC and not an adjustment. Regressing the
# contamination component out of beta is not identifiable from these data:
# genuine biology in which cells come to resemble their neighbours
# (convergent programs, trogocytosis, contact-induced ligands) produces the
# same cross-gene signature, and differs only in mechanism. An earlier idea --
# using the bandwidth profile as the instrument, on the theory that
# misassignment acts at a cell diameter while signalling acts over tens of
# microns -- was tested and does NOT work: the standardised artefact does not
# decay with sigma but grows with it, because the small-sigma kernel density
# estimate is noisy and attenuates the coefficient more than the short range of
# the contamination concentrates it. Absent an instrument, subtracting the
# predicted component would delete real biology in precisely the case of
# interest, so this function reports and does not correct.

#' Score niche associations for transcript-spillover contamination
#'
#' Reports, for each (index cell type, niche cell type, bandwidth), the rank
#' correlation across genes between the fitted niche coefficient and the log
#' ratio of the niche cell type's expression to the index cell type's. A high
#' score means the associations for that pair are ordered by how much of a
#' marker each gene is for the neighbouring cell type — the signature of
#' transcript misassignment rather than biology.
#'
#' The comparison uses the fit's own cell-means `CellType` coefficients, so it
#' needs no second pass over the counts. Spearman correlation is used because
#' the relationship is monotone but saturating: contamination cannot raise a
#' gene above the neighbour's own expression level.
#'
#' **Interpretation.** The score is a flag, not a correction, and it has no
#' universal threshold — read it against the `"twoway"` rows and against your
#' own negative controls. On simulated data with a 30% leak rate the tested
#' rows scored ≈ 0.8 in niche mode against ≈ 0 with no leakage, while condition
#' mode's tested rows stayed ≈ 0 in both cases (the contamination is absorbed
#' by the untested two-way term, which scored ≈ 0.6).
#'
#' A high score does not prove artefact: genuine biology in which cells come to
#' resemble their neighbours produces the same ordering. It means the result
#' cannot be distinguished from spillover on this evidence alone, and wants an
#' orthogonal check — a nucleus-only or alternative segmentation, negative
#' control probes, or an independent assay.
#'
#' @param object a [SpiDEResults] from [fitSpiDE()] or [spiDE()].
#' @param type which coefficients to score. `"tested"` (default) scores the
#'   columns the mode actually tests — the two-way `CellType:niche` terms in a
#'   condition-free fit, the three-way `CellType:condition:niche` terms
#'   otherwise. `"twoway"` scores the two-way terms whatever the mode, which in
#'   condition mode is the untested nuisance block that absorbs spillover.
#'   `"both"` returns both, distinguished by the `term` column.
#' @param ... ignored.
#'
#' @return a data.frame with one row per (bandwidth, index, niche, term) and
#'   columns `bandwidth`, `ct_index`, `ct_niche`, `term`, `score` (Spearman
#'   correlation) and `n_genes`.
#'
#' @examples
#' data(toySpiDE)
#' spe <- buildNiches(toySpiDE, sigma = 20)
#' res <- fitSpiDE(spe, condition = NULL, sigma = 20, verbose = FALSE)
#' spilloverScore(res)
#'
#' @rdname spilloverScore
#' @importFrom stats cor
#' @export
setMethod(
  "spilloverScore", "SpiDEResults",
  function(object, type = c("tested", "twoway", "both"), ...) {
    type <- match.arg(type)
    fl <- fits(object)
    if (!length(fl)) {
      stop("`object` holds no fits")
    }

    wanted <- switch(type,
                     tested = "tested",
                     twoway = "twoway",
                     both   = c("tested", "twoway"))

    out <- lapply(names(fl), function(nm) {
      f <- fl[[nm]]
      md <- .fitMode(f)
      ct <- as.character(f@covtype)
      cm <- f@coefmap
      alpha <- f@alpha

      rows <- lapply(wanted, function(tm) {
        sel <- if (tm == "tested") {
          which(.nicheTestCols(ct, md))
        } else {
          which(ct == "Niche")
        }
        # In niche mode the two sets are the same block; returning it twice
        # under different labels would be misleading, so "twoway" is dropped.
        if (tm == "twoway" && identical(md, "niche") && "tested" %in% wanted) {
          return(NULL)
        }
        if (!length(sel)) {
          return(NULL)
        }
        do.call(rbind, lapply(sel, function(k) {
          ci <- cm$index[k]
          ni <- cm$niche[k]
          cc <- paste0("CellType", ci)
          cn <- paste0("CellType", ni)
          # A merged niche column has no CellType counterpart, so no ratio can
          # be formed for it; report NA rather than guessing a member.
          if (!all(c(cc, cn) %in% colnames(alpha))) {
            sc <- NA_real_
          } else {
            # cell-means coding: alpha[, CellType<k>] is the log expected count
            # for cell type k, so the difference is the log marker ratio
            lr <- alpha[, cn] - alpha[, cc]
            sc <- suppressWarnings(
              stats::cor(alpha[, cm$covariate[k]], lr, method = "spearman",
                         use = "complete.obs"))
          }
          data.frame(bandwidth = f@sigma, ct_index = ci, ct_niche = ni,
                     term = tm, score = sc, n_genes = nrow(alpha),
                     stringsAsFactors = FALSE)
        }))
      })
      do.call(rbind, rows)
    })

    out <- do.call(rbind, out)
    if (is.null(out)) {
      return(data.frame(bandwidth = numeric(0), ct_index = character(0),
                        ct_niche = character(0), term = character(0),
                        score = numeric(0), n_genes = numeric(0),
                        stringsAsFactors = FALSE))
    }
    out <- out[order(-abs(out$score), out$bandwidth), , drop = FALSE]
    rownames(out) <- NULL
    out
  }
)
