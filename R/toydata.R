# Internal synthetic SpatialExperiment used by examples and tests. Not exported.
# Plants a recoverable neighbourhood-dependent signal: in the index cell type
# "A", gene "G1" is up-regulated in the "Responder" condition in proportion to
# the local density of the niche cell type "B" (B cells are concentrated at high
# x, so the B-niche density increases with x). This lets end-to-end tests check
# that spiDE recovers a known ResponseNiche effect.
#
# The count model is a negative-binomial model shared by every generator
# (`.simGeneParams()` + `.simCounts()`), adapting the Splat model of splatter
# package (Zappia, Phipson & Oshlack, Genome Biology 2017; Bioconductor package
# 'splatter'): gene means are drawn from a Gamma so most counts are 0/1/2/3, the
# NB overdispersion follows a downward mean-variance (BCV) trend, and per-cell
# library sizes add sequencing-depth variation -- the properties of real spatial
# data.

#' Per-gene latent parameters for the realistic NB count model
#'
#' Adapts the Splat simulation model of the \pkg{splatter} package (Zappia,
#' Phipson & Oshlack, 2017). Draws a per-gene baseline abundance from a Gamma
#' distribution (splatter's \code{Splat} mean model) and derives a per-gene NB
#' dispersion from a
#' mean-variance trend: the biological coefficient of variation decreases with
#' abundance (\code{bcv = bcv.common + bcv.disp / sqrt(gmean)}), so the NB
#' \code{size = 1 / bcv^2} increases with abundance -- highly expressed genes
#' are less overdispersed, exactly the trend seen in real spatial data. A small
#' fraction of genes are expression outliers (a heavy high-abundance tail).
#' Genes named in \code{boost} are given a high fixed abundance so a planted
#' effect on them stays estimable in a small fixture.
#'
#' @param gene_names character vector of gene identifiers.
#' @param mean.shape,mean.rate Gamma shape/rate for the per-gene baseline mean.
#'   Defaults give a median per-cell mean well below 1 (counts mostly 0/1/2/3).
#' @param bcv.common,bcv.disp intercept and 1/sqrt(mean) coefficient of the BCV
#'   trend that sets the per-gene NB dispersion.
#' @param out.prob,out.facLoc,out.facScale fraction of expression-outlier genes
#'   and the log-normal (meanlog/sdlog) multiplier applied to their abundance.
#' @param mean.min a lower floor on the per-gene abundance. Zero for a fully
#'   realistic draw; a small positive value avoids degenerate all-zero genes in
#'   tiny fixtures (which break edgeR's dispersion step). Real pipelines
#'   filter such genes; the study simulator does the same downstream.
#' @param boost gene names assigned the fixed abundance \code{boost.gmean}.
#' @param boost.gmean the fixed baseline abundance for boosted genes.
#' @param total.count optional target expected library size (total counts per
#'   cell). When set, the per-gene means are rescaled to relative expression
#'   proportions that sum to \code{total.count} (splatter's library-size step),
#'   so the sequencing depth is controlled independently of the gene count and
#'   of the expression outliers. Leave \code{NULL} (default) for the tiny toy
#'   fixtures, whose few genes must stay in the low-count regime.
#' @param gmean optional pre-specified per-gene mean vector. When supplied the
#'   Gamma draw, outliers, \code{mean.min}, \code{boost} and \code{total.count}
#'   steps are all skipped and this abundance is used directly; only the BCV
#'   trend (\code{size}) is derived from it. Lets a caller impose its own
#'   expression structure (e.g. cell-type markers) while keeping the
#'   dispersion trend single-sourced.
#' @return a list with per-gene \code{gmean}, \code{bcv} and NB \code{size}.
#' @references Zappia L, Phipson B, Oshlack A (2017). "Splatter: simulation of
#'   single-cell RNA sequencing data." \emph{Genome Biology} 18:174.
#' @importFrom stats rgamma rlnorm runif
#' @noRd
.simGeneParams <- function(gene_names, mean.shape = 0.35, mean.rate = 0.9,
                           bcv.common = 0.3, bcv.disp = 1, out.prob = 0.05,
                           out.facLoc = 2, out.facScale = 0.5, mean.min = 0,
                           boost = character(0), boost.gmean = 6,
                           total.count = NULL, gmean = NULL) {
  ng <- length(gene_names)
  if (is.null(gmean)) {
    gmean <- pmax(rgamma(ng, shape = mean.shape, rate = mean.rate), 1e-3)
    # expression outliers: a few highly expressed genes (the heavy right tail)
    is_out <- runif(ng) < out.prob
    if (any(is_out)) {
      gmean[is_out] <- gmean[is_out] *
        rlnorm(sum(is_out), meanlog = out.facLoc, sdlog = out.facScale)
    }
    gmean <- pmax(gmean, mean.min)
    names(gmean) <- gene_names
    # planted/boosted genes get a high, fixed abundance so effects stay estimable
    if (length(boost)) gmean[boost] <- boost.gmean
    # scale relative expression to a target library size (controls sequencing
    # depth independently of gene count / outliers, as splatter does)
    if (!is.null(total.count)) gmean <- gmean / sum(gmean) * total.count
  } else {
    gmean <- pmax(gmean, 1e-3)
    names(gmean) <- gene_names
  }
  # trended overdispersion: BCV falls with abundance, so NB size rises with it
  bcv <- bcv.common + bcv.disp / sqrt(gmean)
  list(gmean = gmean, bcv = bcv, size = 1 / bcv^2)
}

#' Draw NB counts from the realistic mean-variance model
#'
#' Assembles the per-cell mean as
#' \code{mu = lib.size * gmean * exp(delta_ct) * exp(log_effect)} and draws
#' \code{rnbinom(mu, size)} with the per-gene trended \code{size} from
#' \code{.simGeneParams()}. \code{delta} are per-(gene, cell-type) log-fold
#' offsets (gene-specific, cell-type-specific means); \code{log_effect} carries
#' any planted interaction or random-effect terms supplied by the caller;
#' \code{lib.size} adds per-cell sequencing-depth variation.
#'
#' @param gene_params the list returned by \code{.simGeneParams()}.
#' @param cell_type character vector of length ncol giving each cell's type.
#' @param log_effect optional genes x cells matrix of additive log-mean effects.
#' @param lib.size optional per-cell multiplicative library size (default 1).
#' @param sd.ct SD of the per-(gene, cell-type) log-fold offsets.
#' @return a genes x cells integer counts matrix.
#' @importFrom stats rnbinom rnorm
#' @noRd
.simCounts <- function(gene_params, cell_type, log_effect = NULL,
                       lib.size = NULL, sd.ct = 0.5) {
  gmean <- gene_params$gmean
  size <- gene_params$size
  ng <- length(gmean)
  gene_names <- names(gmean)
  cts <- sort(unique(cell_type))
  # per-(gene, cell-type) log-fold offsets: cell types differ, gene by gene
  delta <- matrix(rnorm(ng * length(cts), 0, sd.ct), ng, length(cts),
                  dimnames = list(gene_names, cts))
  lmu <- log(gmean) + delta[, cell_type, drop = FALSE]
  if (!is.null(log_effect)) lmu <- lmu + log_effect
  if (!is.null(lib.size)) lmu <- sweep(lmu, 2, log(lib.size), `+`)
  mu <- exp(lmu)
  matrix(rnbinom(length(mu), mu = as.vector(mu), size = size),
         nrow = ng, dimnames = list(gene_names, NULL))
}

# Set the RNG seed for the duration of the caller only.
#
# These generators take a `seed` argument precisely so their output is
# reproducible, but a bare set.seed() in package code silently reseeds the
# USER's global RNG stream -- so calling a toy-data helper would perturb every
# random draw a user made afterwards. Restoring the previous .Random.seed on
# exit keeps the reproducibility and drops the side effect.
#
# BiocCheck flags the set.seed() below by a plain text search. It is the one
# call that must stay: this helper exists precisely to make seeding safe, and
# it is reached only when the caller passed an explicit `seed`.
#
# Callers must not call on.exit() without add = TRUE afterwards: that would
# discard the restore handler registered here and silently reinstate the very
# side effect this exists to remove.
.localSeed <- function(seed, env = parent.frame()) {
  if (is.null(seed)) return(invisible(NULL))
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  .spiDE_prev_seed <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  # Registered in the CALLER's frame, so the restore happens when the generator
  # returns rather than when this helper does.
  do.call(on.exit, list(quote(
    if (is.null(.spiDE_prev_seed)) {
      suppressWarnings(rm(".Random.seed", envir = globalenv()))
    } else {
      assign(".Random.seed", .spiDE_prev_seed, envir = globalenv())
    }), add = TRUE), envir = env)
  assign(".spiDE_prev_seed", .spiDE_prev_seed, envir = env)
  set.seed(seed)
}

#' Library-size variation that survives within-sample normalisation
#'
#' Counts entering spiDE have already been library-size corrected \emph{within}
#' a sample, so there is no per-cell sequencing-depth noise left to simulate.
#' Simulating it anyway (an i.i.d. per-cell log-normal, as these fixtures once
#' did) makes the fitted model \strong{misspecified}, which contaminates any
#' calibration conclusion drawn from them.
#'
#' What genuinely remains is two gene-independent scalings: a \strong{sample
#' level} depth bias, which a per-sample random intercept absorbs, and
#' \strong{cell-type level} differences in total RNA that are real biology
#' rather than technical (a tumour cell carries far more transcript than a
#' lymphocyte), which the \code{CellType} main effects absorb. Because both are
#' absorbed exactly, the fitted NB model is correctly specified and a
#' miscalibration is attributable to the inference rather than to the fixture.
#'
#' @param sample_id,cell_type per-cell labels.
#' @param sd.sample,sd.celltype log-scale SDs of the two components.
#' @return a per-cell multiplicative library size.
#' @noRd
.simLibSize <- function(sample_id, cell_type, sd.sample = 0.15,
                        sd.celltype = 0.35) {
  smp <- as.character(sample_id)
  ct <- as.character(cell_type)
  s_lib <- stats::setNames(
    stats::rnorm(length(unique(smp)), 0, sd.sample), sort(unique(smp)))
  ct_lib <- stats::setNames(
    stats::rnorm(length(unique(ct)), 0, sd.celltype), sort(unique(ct)))
  exp(s_lib[smp] + ct_lib[ct])
}

#' Build a small synthetic SpatialExperiment for examples and tests
#'
#' @param n_samples number of samples (half Responder, half Non-responder).
#' @param n_per number of cells per sample.
#' @param n_genes number of genes.
#' @param field the side length of the (square) spatial field.
#' @param beta effect size of the planted B-niche x Responder interaction on G1.
#'   This is the MAXIMUM log-fold-change across the field, since the effect is
#'   planted as \code{beta * (x / field)} -- so \code{beta = 2} is a 7-fold
#'   gradient and \code{beta = 4} a 55-fold one. The response to it is
#'   \strong{non-monotonic}: G1's own dynamic range inflates its estimated NB
#'   dispersion, which inflates its standard error, so past \code{beta ~ 2} the
#'   induced noise outgrows the signal and the statistic collapses (measured
#'   t: 2.22, 2.87, 5.28, 2.49, 1.42, 0.25 at beta 1, 1.5, 2, 2.5, 3, 4). The
#'   default is the peak; do not raise it expecting a stronger effect.
#' @param seed random seed.
#' @return a SpatialExperiment with counts, cell_type, sample_id, condition, Age,
#'   Area, and spatial coordinates.
#' @importFrom stats rlnorm rnorm runif ave
#' @noRd
.toySPE <- function(n_samples = 6, n_per = 80, n_genes = 20, field = 500,
                    sd.lib.sample = 0.15, sd.lib.celltype = 0.35,
                    beta = 2, seed = 1) {
  .localSeed(seed)
  gene_names <- sprintf("G%d", seq_len(n_genes))

  sample_ids <- sprintf("S%d", seq_len(n_samples))
  # sample-level condition: balanced two levels
  cond_levels <- rep(c("Responder", "Non-responder"), length.out = n_samples)
  names(cond_levels) <- sample_ids

  cells <- lapply(sample_ids, function(sid) {
    # coordinates on a micron-scale field so bandwidths of 10-70 are meaningful
    x <- runif(n_per, 0, field)
    y <- runif(n_per, 0, field)
    # B cells concentrated in the right 40% of the field; A and C elsewhere
    ct <- ifelse(
      x > 0.6 * field & runif(n_per) < 0.7, "B",
      sample(c("A", "C"), n_per, replace = TRUE)
    )
    data.frame(
      sample_id = sid,
      condition = cond_levels[[sid]],
      x = x, y = y,
      cell_type = ct,
      Age = rnorm(1), # sample-level nuisance, constant within sample
      stringsAsFactors = FALSE
    )
  })
  cd <- do.call(rbind, cells)
  cd$Age <- ave(cd$Age, cd$sample_id) # ensure constant within sample
  n <- nrow(cd)
  cd$cell_id <- sprintf("cell%d", seq_len(n))

  # Realistic low-count model (~60% zeros, counts mostly 0/1/2/3) with a mild
  # dispersion trend; G1 is boosted and the effect is strong so the planted
  # signal stays recoverable (|t| > 5) in this small fixture, and mean.min keeps
  # the 20-gene edgeR dispersion step free of degenerate near-zero genes.
  gene_params <- .simGeneParams(gene_names, bcv.disp = 0.6, mean.min = 0.5,
                                boost = "G1", boost.gmean = 8)

  # planted signal: G1 in A cells scales with x (B-niche proxy) in Responders
  is_A <- cd$cell_type == "A"
  is_resp <- cd$condition == "Responder"
  log_effect <- matrix(0, n_genes, n, dimnames = list(gene_names, cd$cell_id))
  log_effect["G1", ] <- beta * is_A * is_resp * (cd$x / field)

  lib.size <- .simLibSize(cd$sample_id, cd$cell_type,
                          sd.sample = sd.lib.sample,
                          sd.celltype = sd.lib.celltype)

  counts <- .simCounts(gene_params, cd$cell_type, log_effect = log_effect,
                       lib.size = lib.size)
  colnames(counts) <- cd$cell_id

  # Area is a genuine per-cell measurement, independent of the counts. It was
  # previously an alias for the library size, hence near-collinear with the
  # cell-type indicators (which sum to 1 in every row).
  cd$Area <- rlnorm(n, meanlog = 0, sdlog = 0.25)
  cd$nCount <- colSums(counts)

  SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(cd),
    spatialCoords = as.matrix(cd[, c("x", "y")])
  )
}

#' Build a small synthetic SpatialExperiment with NO condition
#'
#' Like \code{.toySPE()} but plants a condition-free neighbourhood effect: in
#' index cell type "A", gene "G1" is up-regulated in proportion to the local
#' density of niche cell type "B" (B cells concentrate at high \code{x}), for
#' every cell rather than only in responders. There is deliberately no
#' \code{condition} column at all, so a fixture built here can only be analysed
#' in the condition-free mode -- which is what the niche-only tests must prove,
#' rather than merely leaving a present condition unread.
#'
#' @inheritParams .toySPE
#' @param beta effect size of the planted B-niche association on G1 in A cells.
#' @return a SpatialExperiment with counts, cell_type, sample_id, Age, Area,
#'   nCount and spatial coordinates -- and no condition.
#' @importFrom stats rlnorm rnorm runif ave
#' @noRd
.toyNiche <- function(n_samples = 6, n_per = 80, n_genes = 20, field = 500,
                      sd.lib.sample = 0.15, sd.lib.celltype = 0.35,
                      beta = 2.5, seed = 1) {
  .localSeed(seed)
  gene_names <- sprintf("G%d", seq_len(n_genes))
  sample_ids <- sprintf("S%d", seq_len(n_samples))

  cells <- lapply(sample_ids, function(sid) {
    x <- runif(n_per, 0, field)
    y <- runif(n_per, 0, field)
    # B cells concentrated in the right 40% of the field; A and C elsewhere
    ct <- ifelse(x > 0.6 * field & runif(n_per) < 0.7, "B",
                 sample(c("A", "C"), n_per, replace = TRUE))
    data.frame(sample_id = sid, x = x, y = y, cell_type = ct,
               Age = rnorm(1), stringsAsFactors = FALSE)
  })
  cd <- do.call(rbind, cells)
  cd$Age <- ave(cd$Age, cd$sample_id) # ensure constant within sample
  n <- nrow(cd)
  cd$cell_id <- sprintf("cell%d", seq_len(n))

  gene_params <- .simGeneParams(gene_names, bcv.disp = 0.6, mean.min = 0.5,
                                boost = "G1", boost.gmean = 8)

  # planted signal: G1 in A cells scales with x (the B-niche proxy). No
  # condition gates it -- this is a pure CellType:niche effect.
  is_A <- cd$cell_type == "A"
  log_effect <- matrix(0, n_genes, n, dimnames = list(gene_names, cd$cell_id))
  log_effect["G1", ] <- beta * is_A * (cd$x / field)

  lib.size <- .simLibSize(cd$sample_id, cd$cell_type,
                          sd.sample = sd.lib.sample,
                          sd.celltype = sd.lib.celltype)
  counts <- .simCounts(gene_params, cd$cell_type, log_effect = log_effect,
                       lib.size = lib.size)
  colnames(counts) <- cd$cell_id

  # Area is a genuine per-cell measurement, independent of the counts. It was
  # previously an alias for the library size, hence near-collinear with the
  # cell-type indicators (which sum to 1 in every row).
  cd$Area <- rlnorm(n, meanlog = 0, sdlog = 0.25)
  cd$nCount <- colSums(counts)

  SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(cd),
    spatialCoords = as.matrix(cd[, c("x", "y")])
  )
}

#' Synthetic data with patient-level clustering (for the mixed-effects tests)
#'
#' Like \code{.toySPE()} but plants a per-(gene, sample) random intercept shared
#' by every cell of a sample, and NO response effect. Because the response is a
#' sample-level label, treating cells as independent makes the response tests
#' anti-conservative (cell-level pseudo-replication); the random-intercept fit
#' should test the response effect against between-sample variability and recover
#' calibration. \code{sd_patient} is the planted between-sample SD (variance
#' \code{sd_patient^2} is what \code{random = "intercept"} should recover). Genes
#' here are moderately expressed (not the very-low regime of \code{.toySPE()}) so
#' the between-sample variance component is estimable in a small fixture.
#'
#' @inheritParams .toySPE
#' @param sd_patient the planted between-sample (patient) intercept SD.
#' @return a SpatialExperiment with a patient-clustered null signal.
#' @importFrom stats rlnorm rnorm runif
#' @noRd
.toyClustered <- function(n_samples = 8, n_per = 80, n_genes = 30, field = 500,
                          sd.lib.sample = 0.15, sd.lib.celltype = 0.35,
                          sd_patient = 0.7, seed = 1) {
  .localSeed(seed)
  gene_names <- sprintf("G%d", seq_len(n_genes))
  sample_ids <- sprintf("S%d", seq_len(n_samples))
  cond_levels <- rep(c("Responder", "Non-responder"), length.out = n_samples)
  names(cond_levels) <- sample_ids

  cells <- lapply(sample_ids, function(sid) {
    x <- runif(n_per, 0, field)
    y <- runif(n_per, 0, field)
    ct <- ifelse(x > 0.6 * field & runif(n_per) < 0.7, "B",
                 sample(c("A", "C"), n_per, replace = TRUE))
    data.frame(sample_id = sid, condition = cond_levels[[sid]],
               x = x, y = y, cell_type = ct, stringsAsFactors = FALSE)
  })
  cd <- do.call(rbind, cells)
  n <- nrow(cd)
  cd$cell_id <- sprintf("cell%d", seq_len(n))

  # Moderate, tightly-spread abundances with a mild dispersion trend: the shared
  # between-sample variance component (tau2) is estimated across all genes, so a
  # few highly-overdispersed low-mean genes would inflate it. This fixture
  # exercises the mixed-effects machinery, not the low-count realism.
  gene_params <- .simGeneParams(gene_names, mean.shape = 8, mean.rate = 1.5,
                                bcv.disp = 0.3, out.prob = 0, mean.min = 1)
  # per-(gene, sample) random intercept applied to every cell of the sample
  u <- matrix(rnorm(n_genes * n_samples, 0, sd_patient), nrow = n_genes,
              dimnames = list(gene_names, sample_ids))
  log_effect <- u[, cd$sample_id, drop = FALSE]

  lib.size <- .simLibSize(cd$sample_id, cd$cell_type,
                          sd.sample = sd.lib.sample,
                          sd.celltype = sd.lib.celltype)
  counts <- .simCounts(gene_params, cd$cell_type, log_effect = log_effect,
                       lib.size = lib.size)
  colnames(counts) <- cd$cell_id

  # Area is a genuine per-cell measurement, independent of the counts. It was
  # previously an alias for the library size, hence near-collinear with the
  # cell-type indicators (which sum to 1 in every row).
  cd$Area <- rlnorm(n, meanlog = 0, sdlog = 0.25)
  cd$nCount <- colSums(counts)

  SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(cd),
    spatialCoords = as.matrix(cd[, c("x", "y")])
  )
}
