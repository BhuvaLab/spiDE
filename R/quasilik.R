# Quasi-likelihood dispersion, edgeR v4 style.
#
# NOT WIRED INTO INFERENCE, DELIBERATELY. It was built to be, and then measured
# first: on the real cohort it centres the scale but STEEPENS the expression
# gradient it was meant to remove (1.179 -> 1.263), because the ratio it applies
# over the Pearson dispersion falls from 1.178 at the lowest abundance decile to
# 1.049 at the highest. See research/fdr-ordering/FINDINGS.md, "The QL dispersion
# does not fix it", and the spec's REFUTED status. This file is kept because it
# is correct, oracle-tested against edgeR, and is the instrument that result was
# measured with -- not because it is on a path to being used.
#
# WHY THIS EXISTS AT ALL. spiDE scales each gene's Wald covariance by that
# gene's working Pearson dispersion (R/inference.R). Gate 0
# (research/fdr-ordering/FINDINGS.md) measured what that leaves behind: the
# per-gene null t is over-dispersed by up to 2.3x, the excess sits on the
# dispersion/abundance axis, and the correction under-corrects. The naive
# alternative is worse, not better -- raw deviance / (n - p) has median 0.267
# against Pearson's 0.837 on the real cohort, and since the scale MULTIPLIES the
# covariance that would inflate t by 1.94x instead of 1.09x.
#
# The reason is that E[unit deviance] is nowhere near 1 at small counts. edgeR
# v4 fixes exactly this, per observation, and the fix is NOT the Chebyshev
# tables in its src/ql_weights.c -- those are a fast path. The definition is the
# phi >= 4.001 branch: obtain the first two moments of the unit deviance under
# the fitted NB by direct summation over the pmf, then match them to a scaled
# chi-square. If d ~ c * chisq_nu then E[d] = c*nu and Var[d] = 2 c^2 nu, so
#
#     nu = 2 E[d]^2 / Var[d]        the observation's effective df   (w1)
#     1/c = 2 E[d] / Var[d]         rescales its unit deviance       (w0)
#
# and s2 = sum(d_i * w0_i) / sum((1 - h_i) * w1_i). Implementing the definition
# rather than the approximation keeps this pure R -- a Bioconductor package may
# not call edgeR's internal C -- and `tests/testthat/test-ql.R` holds it to
# edgeR::glmQLFit as an oracle.

#' First two moments of the NB unit deviance, by direct summation
#'
#' Reproduces the definition edgeR's \code{compute_weight()} approximates with
#' Chebyshev tables. The summation window is centred on \code{mu} and sized from
#' the NB sd rather than fixed at edgeR's 50 terms from zero: edgeR can truncate
#' there because that branch is only reached when \code{phi >= 4}, where the pmf
#' piles up near zero. spiDE reaches the same code with any dispersion.
#'
#' Elements are grouped by the width they need so a few high-mu cells do not
#' impose their window on everything, and each group accumulates in two passes
#' so memory stays O(length(mu)) rather than O(width * length(mu)).
#'
#' @param mu vector of fitted means.
#' @param phi NB dispersion: scalar, or the same length as \code{mu}.
#' @param eps tail probability left outside the summation window at each end.
#' @param maxterms hard cap on the number of pmf terms summed per element.
#' @return list with \code{w0} (deviance rescaling) and \code{w1} (effective df).
#' @noRd
.nbDevianceMoments <- function(mu, phi, eps = 1e-10, maxterms = 20000L) {
  mu <- as.numeric(mu)
  phi <- rep_len(as.numeric(phi), length(mu))
  w0 <- w1 <- numeric(length(mu))
  ok <- which(is.finite(mu) & mu > 1e-32 & is.finite(phi) & phi > 0)
  if (!length(ok)) return(list(w0 = w0, w1 = w1))

  # Window from the pmf's own quantiles, not from mu +/- k*sd: at spiDE's
  # operating point (mu ~ 0.1, phi ~ 5, so size ~ 0.2) the NB is far too
  # heavy-tailed for an sd-based window -- at mu = 1, phi = 20 a 12-sd window
  # gets w1 wrong by 41%.
  sz <- 1 / phi[ok]
  lo <- stats::qnbinom(eps, size = sz, mu = mu[ok])
  hi <- stats::qnbinom(1 - eps, size = sz, mu = mu[ok])
  wid <- pmin(maxterms, pmax(10, hi - lo))
  grp <- split(seq_along(ok), pmin(30L, ceiling(log2(wid))))

  for (ix in grp) {
    m <- mu[ok][ix]; p <- phi[ok][ix]; size <- 1 / p
    l <- lo[ix]; K <- max(wid[ix])
    unit <- function(i) {
      a <- ifelse(i > 0, i * log(pmax(i, 1) / m), 0)
      2 * (a - (i + size) * log((i + size) / (m + size)))
    }
    ed <- vd <- numeric(length(ix))
    for (j in seq_len(K + 1L) - 1L) {           # pass 1: E[d]
      i <- l + j
      ed <- ed + stats::dnbinom(i, size = size, mu = m) * unit(i)
    }
    for (j in seq_len(K + 1L) - 1L) {           # pass 2: Var[d]
      i <- l + j
      vd <- vd + stats::dnbinom(i, size = size, mu = m) * (unit(i) - ed)^2
    }
    good <- is.finite(ed) & is.finite(vd) & vd > 0
    w0[ok[ix][good]] <- (2 * ed / vd)[good]
    w1[ok[ix][good]] <- (2 * ed * ed / vd)[good]
  }
  list(w0 = w0, w1 = w1)
}

#' NB unit deviance
#' @noRd
.nbUnitDeviance <- function(y, mu, phi) {
  size <- 1 / phi
  a <- ifelse(y > 0, y * log(pmax(y, .Machine$double.xmin) / mu), 0)
  2 * (a - (y + size) * log((y + size) / (mu + size)))
}

#' Per-gene quasi-likelihood dispersion
#'
#' @param y counts, genes x cells.
#' @param mu fitted means, same shape.
#' @param phi NB dispersion: scalar or per gene.
#' @param design the design matrix (cells x p), used for the leverages.
#' @param prior the average quasi-dispersion edgeR divides through by; 1 leaves
#'   the parameterisation alone.
#' @param moments \code{"grid"} evaluates the deviance moments on a log-mu grid
#'   per gene and interpolates onto the cells -- the moments depend only on
#'   (mu, phi) and phi is constant within a gene, so this is an interpolation of
#'   a smooth function, not a sample. It takes the real cohort from ~5e10 pmf
#'   evaluations to ~1e9. \code{"cell"} evaluates at every cell.
#' @param ngrid grid points per gene for \code{moments = "grid"}.
#' @param leverage \code{"trace"} spreads the p degrees of freedom evenly, which
#'   is exact to O(p/n) and is what spiDE's n >> p designs want;
#'   \code{"exact"} forms per-observation hat values (O(n p^2) per gene) and is
#'   what reproduces edgeR on its own small-n designs.
#' @return list of per-gene \code{deviance}, \code{df} and \code{s2}.
#' @noRd
.qlDispersion <- function(y, mu, phi, design, prior = 1,
                          leverage = c("trace", "exact"),
                          moments = c("grid", "cell"), ngrid = 256L) {
  leverage <- match.arg(leverage); moments <- match.arg(moments)
  y <- as.matrix(y); mu <- as.matrix(mu)
  ng <- nrow(y); n <- ncol(y); p <- ncol(design)
  phi <- rep_len(as.numeric(phi), ng)
  dev <- df <- numeric(ng)

  for (g in seq_len(ng)) {
    mg <- mu[g, ]; yg <- y[g, ]
    phi_eff <- phi[g] / prior                       # the deviance's dispersion
    h <- if (leverage == "exact") {
      zw <- sqrt(mg / (1 + mg * phi_eff))
      qrx <- qr(design * zw)
      rowSums(qr.Q(qrx)^2)
    } else {
      rep(p / n, n)
    }
    # moments are taken at mu/prior with the UNSCALED phi, matching edgeR's
    # compute_weight(u, phi, prior) call convention.
    ms <- mg / prior
    m <- if (moments == "cell" || length(unique(ms)) <= ngrid) {
      .nbDevianceMoments(ms, phi[g])
    } else {
      lg <- log(pmax(ms, 1e-32))
      kn <- seq(min(lg), max(lg), length.out = ngrid)
      mk <- .nbDevianceMoments(exp(kn), phi[g])
      list(w0 = stats::approx(kn, mk$w0, lg, rule = 2)$y,
           w1 = stats::approx(kn, mk$w1, lg, rule = 2)$y)
    }
    udp <- .nbUnitDeviance(yg, mg, phi_eff)
    hdp <- 1 - h
    drop <- !is.finite(hdp) | hdp < 1e-4
    udp[drop] <- 0; hdp[drop] <- 0
    dev[g] <- sum(udp * m$w0, na.rm = TRUE)
    df[g] <- sum(hdp * m$w1, na.rm = TRUE)
  }
  s2 <- ifelse(df < 1e-4, 0, dev / df)
  list(deviance = dev, df = df, s2 = s2)
}
