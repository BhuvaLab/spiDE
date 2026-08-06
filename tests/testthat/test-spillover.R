# spilloverScore(): does the diagnostic actually detect transcript
# misassignment, and does it stay quiet when there is none?
#
# The contrast IS the test. A score that is high under contamination proves
# nothing on its own -- the same ordering could come from cell-type structure
# in the data -- so every assertion here is paired against a matched kappa = 0
# run built from the same seed.

# Strong per-(gene, cell type) markers (sd.ct = 1.5) so the niche type's
# profile is distinctive; that is what spillover transfers. No biological
# niche effect at all, so any ordering by marker ratio is contamination.
.spillFixture <- function(seed = 11, n_genes = 30) {
  set.seed(seed)
  field <- 500
  n_per <- 120
  ids <- sprintf("S%d", 1:4)
  gene_names <- sprintf("G%d", seq_len(n_genes))
  cells <- lapply(ids, function(sid) {
    x <- runif(n_per, 0, field)
    y <- runif(n_per, 0, field)
    ct <- ifelse(x > 0.6 * field & runif(n_per) < 0.7, "B",
                 sample(c("A", "C"), n_per, replace = TRUE))
    data.frame(sample_id = sid, x = x, y = y, cell_type = ct,
               stringsAsFactors = FALSE)
  })
  cd <- do.call(rbind, cells)
  cd$cell_id <- sprintf("cell%d", seq_len(nrow(cd)))
  gp <- spiDE:::.simGeneParams(gene_names, bcv.disp = 0.6, mean.min = 0.5,
                               boost = gene_names[1:5], boost.gmean = 8)
  counts <- spiDE:::.simCounts(gp, cd$cell_type, sd.ct = 1.5)
  colnames(counts) <- cd$cell_id
  SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(cd),
    spatialCoords = as.matrix(cd[, c("x", "y")]))
}

test_that(".toySpill conserves total counts and leaves kappa = 0 untouched", {
  spe <- .spillFixture()
  expect_identical(
    SummarizedExperiment::assay(spiDE:::.toySpill(spe, kappa = 0), "counts"),
    SummarizedExperiment::assay(spe, "counts"))

  sp <- spiDE:::.toySpill(spe, kappa = 0.3, radius = 20)
  tot0 <- sum(SummarizedExperiment::assay(spe, "counts"))
  tot1 <- sum(SummarizedExperiment::assay(sp, "counts"))
  # rounding moves the total a little; mixing must not create or destroy mass
  expect_lt(abs(tot1 - tot0) / tot0, 0.02)
  # and it must actually change the data
  expect_false(identical(SummarizedExperiment::assay(sp, "counts"),
                         SummarizedExperiment::assay(spe, "counts")))
})

test_that("spilloverScore separates contaminated from clean niche fits", {
  spe <- .spillFixture()
  clean <- buildNiches(spe, sigma = 20)
  dirty <- buildNiches(spiDE:::.toySpill(spe, kappa = 0.3, radius = 20),
                       sigma = 20)

  r_clean <- fitSpiDE(clean, condition = NULL, sigma = 20, verbose = FALSE)
  r_dirty <- fitSpiDE(dirty, condition = NULL, sigma = 20, verbose = FALSE)

  s_clean <- spilloverScore(r_clean)
  s_dirty <- spilloverScore(r_dirty)

  expect_s3_class(s_clean, "data.frame")
  expect_true(all(c("bandwidth", "ct_index", "ct_niche", "term", "score",
                    "n_genes") %in% colnames(s_clean)))
  # 3 cell types -> 3*3 - 3 self = 6 tested niche columns
  expect_equal(nrow(s_clean), 6L)
  expect_true(all(s_clean$term == "tested"))

  # the contrast: contamination orders the coefficients by marker ratio
  expect_lt(median(s_clean$score), 0.4)
  expect_gt(median(s_dirty$score), 0.6)
  expect_gt(median(s_dirty$score) - median(s_clean$score), 0.3)
})

test_that("condition mode keeps spillover out of the tested three-way term", {
  # The mechanism: the two-way CellType:niche block is a free nuisance
  # parameter there, and leakage is condition-independent, so contamination
  # lands on it and cancels out of the three-way test.
  spe <- .spillFixture()
  cd <- SummarizedExperiment::colData(spe)
  cd$condition <- ifelse(cd$sample_id %in% c("S1", "S2"), "Responder",
                         "Non-responder")
  SummarizedExperiment::colData(spe) <- cd
  dirty <- buildNiches(spiDE:::.toySpill(spe, kappa = 0.3, radius = 20),
                       sigma = 20)

  res <- fitSpiDE(dirty, condition = "condition", sigma = 20, verbose = FALSE)
  sc <- spilloverScore(res, type = "both")

  tested <- sc$score[sc$term == "tested"]   # three-way
  twoway <- sc$score[sc$term == "twoway"]   # nuisance
  expect_true(length(tested) > 0 && length(twoway) > 0)

  # contamination sits in the nuisance block, not the tested one
  expect_gt(median(twoway), median(tested))
  expect_lt(median(tested), 0.4)
})

test_that("spilloverScore reports NA for merged niches with no CellType term", {
  # A merged niche column has no CellType counterpart, so no marker ratio can
  # be formed; the score must be NA rather than silently using a member.
  spe <- buildNiches(.toySPE(), sigma = 20)
  nm <- SingleCellExperiment::reducedDim(spe, "Niche20")
  SingleCellExperiment::reducedDim(spe, "Niche20") <-
    cbind(AC = nm[, "A"] + nm[, "C"], B = nm[, "B"])
  S4Vectors::metadata(spe)[["spiDE_niche_groups"]] <-
    list(Niche20 = list(AC = c("A", "C"), B = "B"))

  res <- fitSpiDE(spe, condition = NULL, sigma = 20, verbose = FALSE)
  sc <- spilloverScore(res)
  expect_true(any(is.na(sc$score[sc$ct_niche == "AC"])))
  expect_false(any(is.na(sc$score[sc$ct_niche == "B"])))
})
