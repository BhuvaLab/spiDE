# spiGSEA(): a gene set's per-patient slope, tested like a gene's.

spe <- .toySPE(n_samples = 16, n_per = 150, n_genes = 12, seed = 7)
res <- spiDE(spe, condition = "condition", sigma = 30, index = "A", procedure = "all",
             verbose = FALSE)
sets <- list(with_g1 = paste0("G", 1:3), without = paste0("G", 4:8), rest = paste0("G", 9:12),
             tiny = "G1")

test_that("spiGSEA returns one row per usable (set, index, niche, test)", {
  gs <- spiGSEA(res, sets, min.size = 3)
  expect_true(all(c("set", "index", "niche", "test", "size", "estimate", "se", "t", "df", "p", "q",
                    "q.global", "n_patients", "direction") %in% colnames(gs)))
  expect_false("tiny" %in% gs$set)
  expect_setequal(unique(gs$test), c("pooled", "condition"))
  nn <- length(res@fit@index$A$niches)
  expect_equal(nrow(gs), 3L * nn * 2L)
  ok <- is.finite(gs$p)
  expect_true(all(gs$p[ok] >= 0 & gs$p[ok] <= 1))
  expect_true(all(gs$q[ok] >= gs$p[ok] - 1e-12))
})

test_that("a set's patient slope is the mean of its genes' scaled slopes, less the background", {
  x <- res@fit@index$A
  j <- 1L
  vp <- sweep(x$v_model, c(2, 3), x$factor, "*")
  tau <- spiDE:::.dlTau2(x$beta, vp)[, j]
  b <- x$beta[, , j]; v <- vp[, , j]
  sig <- sqrt(tau + apply(v, 1, stats::median, na.rm = TRUE))
  u <- b / sig
  g <- match(sets$with_g1, x$genes)
  self <- spiDE:::.setSlopes(x, list(s = g), "self-contained")
  comp <- spiDE:::.setSlopes(x, list(s = g), "competitive")
  expect_equal(self$b[1, , j], unname(colMeans(u[g, ])), tolerance = 1e-12)
  expect_equal(comp$b[1, , j], unname(colMeans(u[g, ]) - colMeans(u[-g, ])), tolerance = 1e-12)
})

test_that("the set holding the planted gene stands out in its niche", {
  gs <- spiGSEA(res, sets[1:3], test = "condition", min.size = 3)
  b <- gs[gs$niche == "B", ]
  expect_lt(b$p[b$set == "with_g1"], min(b$p[b$set != "with_g1"]))
  expect_gt(b$estimate[b$set == "with_g1"], 0)
})

test_that("the pooled set test never looks at the condition", {
  pt <- res@fit@patients
  set.seed(4)
  lab <- stats::setNames(sample(as.character(pt$condition)), pt$patient)
  f <- res@fit
  f@patients$perm <- factor(lab[pt$patient])
  r2 <- testSpiDE(f, condition = "perm", procedure = "all")
  g1 <- spiGSEA(res, sets, test = "pooled", min.size = 3)
  g2 <- spiGSEA(r2, sets, test = "pooled", min.size = 3)
  expect_equal(g2$p, g1$p, tolerance = 1e-12)
})

test_that("spiGSEA refuses the sandwich engine and a condition test that was not run", {
  rs <- spiDE(spe, condition = "condition", sigma = 30, index = "A", engine = "sandwich",
              procedure = "all", verbose = FALSE)
  expect_error(spiGSEA(rs, sets), "slopes engine")
  f0 <- res@fit
  f0@condition <- character()
  pooled_only <- testSpiDE(f0)
  expect_error(spiGSEA(pooled_only, sets, test = "condition"), "condition")
  expect_error(spiGSEA(res, list(paste0("G", 1:4))), "named list")
})

test_that("a set with no usable patients drops out without taking the others with it", {
  r2 <- res
  x <- r2@fit@index$A
  jB <- match("B", x$niches)
  x$beta[1:2, , jB] <- NA
  r2@fit@index$A <- x
  s3 <- list(dead = x$genes[1:2], setA = x$genes[4:6], setB = x$genes[7:9])
  gs <- spiGSEA(r2, s3, min.size = 2, test = "pooled")
  b <- gs[gs$niche == "B", ]
  expect_true(is.na(b$p[b$set == "dead"]))
  expect_true(all(is.finite(b$p[b$set != "dead"])))
})
