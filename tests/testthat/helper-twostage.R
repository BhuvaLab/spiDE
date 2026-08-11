# A toy SPE with a hand-built SpaNormFit in metadata: design [logLS, bio.x,
# bio.y] fit with SpaNorm::fitNB, wtype marking the split. Small and fast;
# exercises the exact contract .spanormComponents() reads.
toy_spanorm_spe <- function(n_genes = 30, sigma = 30) {
  spe <- spiDE:::.toySPE(n_genes = n_genes)
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  xy <- SpatialExperiment::spatialCoords(spe)
  logLS <- log(pmax(colSums(Y), 1))
  W <- cbind(logLS = logLS - mean(logLS),
             bio.x = as.numeric(scale(xy[, 1])),
             bio.y = as.numeric(scale(xy[, 2])))
  f <- SpaNorm::fitNB(Y, W, lambda.a = 0, verbose = FALSE)
  fit <- methods::new("SpaNormFit",
    ngenes = nrow(Y), ncells = ncol(Y), gene.model = "nb",
    df.tps = c(1L, 1L, 1L, 1L), sample.p = 1, lambda.a = c(0, 0),
    batch = NULL, W = W, alpha = f$alpha, gmean = f$gmean, psi = f$psi,
    wtype = factor(c("ls", "biology", "biology")),
    loglik = f$loglik,
    sampling = if (!is.null(f$sampling)) f$sampling
               else factor(rep("fit", ncol(Y))))
  S4Vectors::metadata(spe)$SpaNorm <- fit
  buildNiches(spe, sigma = sigma, verbose = FALSE)
}
