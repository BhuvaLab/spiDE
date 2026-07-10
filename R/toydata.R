# Internal synthetic SpatialExperiment used by examples and tests. Not exported.
# Plants a recoverable neighbourhood-dependent signal: in the index cell type
# "A", gene "G1" is up-regulated in the "Responder" condition in proportion to
# the local density of the niche cell type "B" (B cells are concentrated at high
# x, so the B-niche density increases with x). This lets end-to-end tests check
# that spiDE recovers a known ResponseNiche effect.

#' Build a small synthetic SpatialExperiment for examples and tests
#'
#' @param n_samples number of samples (half Responder, half Non-responder).
#' @param n_per number of cells per sample.
#' @param n_genes number of genes.
#' @param field the side length of the (square) spatial field.
#' @param beta effect size of the planted B-niche x Responder interaction on G1.
#' @param seed random seed.
#' @return a SpatialExperiment with counts, cell_type, sample_id, condition, Age,
#'   Area, and spatial coordinates.
#' @importFrom stats rnbinom runif
#' @noRd
.toySPE <- function(n_samples = 6, n_per = 80, n_genes = 20, field = 500,
                    beta = 1.5, seed = 1) {
  set.seed(seed)
  cell_types <- c("A", "B", "C")
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

  # per-cell-type baseline log-means (genes x cells)
  ct_base <- matrix(rnorm(n_genes * length(cell_types), mean = 1.5, sd = 0.6),
                    nrow = n_genes, dimnames = list(gene_names, cell_types))
  lmu <- ct_base[, cd$cell_type, drop = FALSE]

  # planted signal: G1 in A cells scales with x (B-niche proxy) in Responders
  is_A <- cd$cell_type == "A"
  is_resp <- cd$condition == "Responder"
  lmu["G1", ] <- lmu["G1", ] + beta * is_A * is_resp * (cd$x / field)

  # Area proportional to a per-sample scale (used by computeSizeFactors)
  Area <- exp(rnorm(n, sd = 0.2))

  mu <- exp(lmu)
  counts <- matrix(
    rnbinom(length(mu), mu = as.vector(mu), size = 5),
    nrow = n_genes, dimnames = list(gene_names, cd$cell_id)
  )

  cd$Area <- Area
  cd$nCount <- colSums(counts)

  SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(cd),
    spatialCoords = as.matrix(cd[, c("x", "y")])
  )
}
