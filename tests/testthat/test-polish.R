# Per-gene convergence after fitNB. SpaNorm's multi-gene IRLS shares one cell
# weight vector across genes and stops on the aggregate log-likelihood, so
# bright genes sit 1-4 SE from their own optimum with a dispersion ~1.6x too
# large. These tests pin the polish stage: it never decreases the penalised
# log-likelihood, it reaches a stationary point, it leaves an already-converged
# gene alone, and the Schur-complement solve equals the dense one.

set.seed(11)

# a small design with an indicator block, in the shape .polishFit() sees
toy_design <- function(n = 240, p_x = 6, n_grp = 8) {
  X <- cbind(1, matrix(rnorm(n * (p_x - 1)), n, p_x - 1))
  colnames(X) <- c("(Intercept)", paste0("x", seq_len(p_x - 1)))
  gidx <- rep(seq_len(n_grp), length.out = n)
  Z <- stats::model.matrix(~ 0 + factor(gidx))
  colnames(Z) <- paste0("SampleCellType", seq_len(n_grp))
  W <- cbind(X, Z)
  list(W = W, nested = grepl("^SampleCellType", colnames(W)),
       pen = c(rep(0, p_x), rep(1.7, n_grp)))
}

test_that(".newtonSolver's Schur solve equals the dense solve", {
  d <- toy_design()
  w <- runif(nrow(d$W), 0.2, 3)
  s <- rnorm(ncol(d$W))

  info <- crossprod(d$W * sqrt(w))
  diag(info) <- diag(info) + d$pen
  expect_equal(as.numeric(spiDE:::.newtonSolver(d$W, d$pen, d$nested)$solve(w, s)),
               as.numeric(solve(info, s)), tolerance = 1e-8)
})

test_that(".newtonSolver's xcov equals the X-block of the dense covariance", {
  d <- toy_design()
  w <- runif(nrow(d$W), 0.2, 3)
  info <- crossprod(d$W * sqrt(w))
  diag(info) <- diag(info) + d$pen
  expect_equal(spiDE:::.newtonSolver(d$W, d$pen, d$nested)$xcov(w),
               solve(info)[!d$nested, !d$nested], tolerance = 1e-8,
               ignore_attr = TRUE)
})

test_that(".newtonSolver falls back to the dense path with no nested block", {
  d <- toy_design()
  nested <- rep(FALSE, ncol(d$W))
  w <- runif(nrow(d$W), 0.2, 3)
  s <- rnorm(ncol(d$W))
  info <- crossprod(d$W * sqrt(w))
  diag(info) <- diag(info) + d$pen
  expect_equal(as.numeric(spiDE:::.newtonSolver(d$W, d$pen, nested)$solve(w, s)),
               as.numeric(solve(info, s)), tolerance = 1e-8)
})

test_that(".polishGene reaches a stationary point and raises the log-likelihood", {
  d <- toy_design()
  a_true <- c(1.2, 0.4, -0.3, 0.2, 0, 0.1, rnorm(8, 0, 0.3))
  mu <- exp(d$W %*% a_true)
  y <- rnbinom(length(mu), mu = as.numeric(mu), size = 1 / 0.4)
  a0 <- c(0.2, rep(0, ncol(d$W) - 1))          # a deliberately poor start
  solver <- spiDE:::.newtonSolver(d$W, d$pen, d$nested)

  ll0 <- spiDE:::.nbPenLoglik(y, as.numeric(exp(d$W %*% a0)), 0.4, a0, d$pen)
  out <- spiDE:::.polishGene(y, d$W, a0, 0.4, d$pen, solver,
                             maxit = 50L, tol = 1e-10)

  expect_gt(out$loglik, ll0)
  # penalised score at the polished point, at the polished psi
  mu1 <- as.numeric(exp(d$W %*% out$alpha))
  s1 <- as.numeric(crossprod(d$W, (y - mu1) / (1 + out$psi * mu1))) -
    d$pen * out$alpha
  expect_lt(max(abs(s1)), 1e-4 * max(abs(y)))
  expect_true(is.finite(out$psi) && out$psi > 0)
  expect_gte(out$iterations, 1L)
  expect_false(out$singular)
})

test_that(".polishGene leaves an already-converged gene alone", {
  d <- toy_design()
  a_true <- c(1.0, 0.3, -0.2, 0.1, 0, 0, rnorm(8, 0, 0.2))
  mu <- exp(d$W %*% a_true)
  y <- rnbinom(length(mu), mu = as.numeric(mu), size = 1 / 0.5)
  solver <- spiDE:::.newtonSolver(d$W, d$pen, d$nested)
  # converge once, then polish again from the converged point
  first <- spiDE:::.polishGene(y, d$W, a_true, 0.5, d$pen, solver, 50L, 1e-12)
  again <- spiDE:::.polishGene(y, d$W, first$alpha, first$psi, d$pen, solver,
                               50L, 1e-12)
  expect_equal(again$alpha, first$alpha, tolerance = 1e-3)
  expect_equal(again$psi, first$psi, tolerance = 1e-3)
})

test_that(".polishGene restarts from cell-type means when the start is degenerate", {
  d <- toy_design()
  mu <- exp(d$W %*% c(1.0, rep(0.1, 5), rnorm(8, 0, 0.2)))
  y <- rnbinom(length(mu), mu = as.numeric(mu), size = 1 / 0.4)
  a_bad <- c(-25, rep(0, ncol(d$W) - 1))       # min log mu < -10
  solver <- spiDE:::.newtonSolver(d$W, d$pen, d$nested)
  out <- spiDE:::.polishGene(y, d$W, a_bad, 0.4, d$pen, solver, 50L, 1e-10)

  expect_true(out$restarted)
  expect_gt(min(as.numeric(d$W %*% out$alpha)), -10)
  expect_true(is.finite(out$loglik))
})

test_that(".polishFit is invariant to gene blocking", {
  d <- toy_design()
  ng <- 5
  A0 <- matrix(0, ng, ncol(d$W), dimnames = list(paste0("G", seq_len(ng)),
                                                 colnames(d$W)))
  A0[, 1] <- 0.5
  Y <- t(vapply(seq_len(ng), function(g) {
    mu <- exp(d$W %*% c(1 + 0.2 * g, rep(0.15, 5), rnorm(8, 0, 0.2)))
    rnbinom(nrow(d$W), mu = as.numeric(mu), size = 1 / 0.4)
  }, numeric(nrow(d$W))))
  dimnames(Y) <- list(rownames(A0), NULL)
  psi0 <- rep(0.4, ng)
  re_group <- ifelse(d$nested, "SampleCellTypeInt", NA_character_)

  a <- spiDE:::.polishFit(Y, d$W, A0, psi0, d$pen, re_group,
                          block.size = NULL)
  b <- spiDE:::.polishFit(Y, d$W, A0, psi0, d$pen, re_group,
                          block.size = 2L)
  expect_equal(a$alpha, b$alpha)
  expect_equal(a$psi, b$psi)
  expect_equal(nrow(a$polish), ng)
  expect_setequal(colnames(a$polish),
                  c("iterations", "psi_fitnb", "restarted", "capped", "singular",
                 "psi_bound", "polished"))
  expect_equal(a$polish$psi_fitnb, psi0)
})

test_that(".polishFit flags a degenerate gene and keeps finite values", {
  d <- toy_design()
  # a gene that is zero everywhere: every working weight is ~0, so the
  # information matrix is singular at the start
  Y <- matrix(0L, 1, nrow(d$W), dimnames = list("G1", NULL))
  A0 <- matrix(-40, 1, ncol(d$W), dimnames = list("G1", colnames(d$W)))
  re_group <- ifelse(d$nested, "SampleCellTypeInt", NA_character_)
  out <- spiDE:::.polishFit(Y, d$W, A0, 0.4, d$pen, re_group)
  expect_true(out$polish$singular || out$polish$restarted)
  expect_true(all(is.finite(out$alpha)))
  expect_true(is.finite(out$psi))
})

test_that("the polish stage populates @polish and raises the per-gene log-likelihood", {
  spe <- buildNiches(spiDE:::.toySPE(n_genes = 10, n_per = 50), sigma = 30)
  f0 <- fitSpiDE(spe, "condition", sigma = 30, random = "none", verbose = FALSE)
  f1 <- polishSpiDE(f0, spe, verbose = FALSE)
  a0 <- fits(f0)[[1]]
  a1 <- fits(f1)[[1]]

  expect_null(a0@polish)
  expect_s3_class(a1@polish, "data.frame")
  expect_equal(nrow(a1@polish), a1@ngenes)
  expect_true(all(a1@polish$iterations >= 1))
  expect_true(all(is.finite(a1@psi)) && all(a1@psi > 0))

  # every gene's own penalised log-likelihood is at least as high as fitNB's
  Y <- SummarizedExperiment::assay(spe, "counts")
  pen <- rep(0, ncol(a1@W))
  ll <- function(fit, g) {
    mu <- as.numeric(exp(fit@W %*% fit@alpha[g, ]))
    spiDE:::.nbPenLoglik(Y[g, ], mu, fit@psi[g], fit@alpha[g, ], pen)
  }
  base <- vapply(seq_len(a1@ngenes), function(g) ll(a0, g), numeric(1))
  gains <- vapply(seq_len(a1@ngenes), function(g) ll(a1, g), numeric(1)) - base
  expect_true(all(gains > -1e-6 * abs(base)))
  expect_gt(median(gains), 0)
})

test_that("the fit leaves the coefficients as fitNB returned them", {
  spe <- buildNiches(spiDE:::.toySPE(n_genes = 10, n_per = 50), sigma = 30)
  a <- fitSpiDE(spe, "condition", sigma = 30, random = "none", verbose = FALSE)
  b <- fitSpiDE(spe, "condition", sigma = 30, random = "none", verbose = FALSE)
  expect_identical(fits(a)[[1]]@alpha, fits(b)[[1]]@alpha)
  expect_identical(fits(a)[[1]]@psi, fits(b)[[1]]@psi)
})

test_that("a polished fit still passes validity and testSpiDE runs on it", {
  spe <- buildNiches(spiDE:::.toySPE(n_genes = 10, n_per = 50), sigma = 30)
  f <- polishSpiDE(fitSpiDE(spe, "condition", sigma = 30, random = "none", verbose = FALSE),
                   spe, verbose = FALSE)
  expect_true(validObject(fits(f)[[1]]))
  r <- testSpiDE(f, spe = spe, fdr = 1)
  expect_true(nrow(results(r)) > 0)
  expect_true(all(is.finite(results(r)$t)))
})

test_that(".newtonSolver refuses a nested block that is not a partition", {
  d <- toy_design()
  # break the partition: give one cell membership of two groups
  W <- d$W
  W[1, which(d$nested)[2]] <- 1
  expect_error(spiDE:::.newtonSolver(W, d$pen, d$nested),
               "0/1 indicators partitioning")
})

test_that(".newtonSolver recovers group membership from a float-valued product", {
  # the group index comes from a dot product; a product of 7 arriving as
  # 6.9999999 must not be truncated to group 6
  d <- toy_design(n_grp = 8)
  W <- d$W
  # perturb the indicators just inside the partition tolerance
  W[, d$nested] <- W[, d$nested] * (1 - 1e-10)
  w <- runif(nrow(W), 0.2, 3)
  s <- rnorm(ncol(W))
  info <- crossprod(W * sqrt(w))
  diag(info) <- diag(info) + d$pen
  expect_equal(as.numeric(spiDE:::.newtonSolver(W, d$pen, d$nested)$solve(w, s)),
               as.numeric(solve(info, s)), tolerance = 1e-6)
})

test_that(".polishFit splits into one block per worker when block.size is NULL", {
  d <- toy_design()
  ng <- 6
  A0 <- matrix(0, ng, ncol(d$W), dimnames = list(paste0("G", seq_len(ng)),
                                                 colnames(d$W)))
  A0[, 1] <- 0.5
  Y <- t(vapply(seq_len(ng), function(g) {
    mu <- exp(d$W %*% c(1 + 0.2 * g, rep(0.15, 5), rnorm(8, 0, 0.2)))
    rnbinom(nrow(d$W), mu = as.numeric(mu), size = 1 / 0.4)
  }, numeric(nrow(d$W))))
  dimnames(Y) <- list(rownames(A0), NULL)
  re_group <- ifelse(d$nested, "SampleCellTypeInt", NA_character_)

  # a multi-worker BPPARAM must not collapse to a single block. Fork-based
  # parallelism, so the workers inherit the loaded namespace (a snow cluster
  # cannot see package internals under devtools::load_all()). Two workers,
  # not more: R CMD check caps BiocParallel at two (_R_CHECK_LIMIT_CORES_).
  skip_on_os("windows")
  bp <- BiocParallel::MulticoreParam(2, progressbar = FALSE)
  skip_if_not(BiocParallel::bpnworkers(bp) == 2)

  msgs <- capture_messages(
    par <- spiDE:::.polishFit(Y, d$W, A0, rep(0.4, ng), d$pen, re_group,
                              BPPARAM = bp, verbose = TRUE))
  expect_match(paste(msgs, collapse = " "), "2 blocks")

  serial <- capture_messages(
    ser <- spiDE:::.polishFit(Y, d$W, A0, rep(0.4, ng), d$pen, re_group,
                              verbose = TRUE))
  expect_match(paste(serial, collapse = " "), "1 block")

  # blocking is exact, so the split must not change the answer
  expect_equal(par$alpha, ser$alpha)
  expect_equal(par$psi, ser$psi)
})

test_that("updateObject fills a NULL-prototype slot on an object serialised without it", {
  # Regression: .fillSlots() used attr(object, s) <- value, and
  # `attr(x, "s") <- NULL` REMOVES an attribute rather than setting it. So any
  # slot whose prototype is NULL could never be filled: the object stayed
  # invalid, show() errored, and updateObject() -- the documented repair path --
  # failed on the very objects it exists to repair. Fits can be hours of
  # cluster time, so this must keep working.
  p <- system.file("extdata", "testfits", "fit_toyspe_none.rds", package = "spiDE")
  skip_if(!nzchar(p), "fixture not installed")
  old <- readRDS(p)
  skip_if("polish" %in% names(attributes(old)),
          "fixture was regenerated with the polish slot present")

  u <- updateObject(old)
  expect_true(validObject(u))
  expect_null(u@polish)
  expect_silent(capture.output(show(u)))
  # every slot that WAS present is untouched
  expect_identical(u@alpha, old@alpha)
  expect_identical(u@psi, old@psi)
  expect_identical(u@W, old@W)
})

test_that(".fillSlots can fill every NULL-prototype slot of a stripped object", {
  proto <- methods::new("SpiDEFit")
  nulls <- Filter(function(s) is.null(methods::slot(proto, s)),
                  methods::slotNames("SpiDEFit"))
  expect_true(length(nulls) > 0)
  stripped <- proto
  for (s in nulls) attributes(stripped)[[s]] <- NULL
  filled <- spiDE:::.fillSlots(stripped, "SpiDEFit")
  expect_true(all(nulls %in% names(attributes(filled))))
  expect_true(all(vapply(nulls, function(s) is.null(methods::slot(filled, s)),
                         logical(1))))
})

test_that("non-integer counts are refused before the fit, not after it", {
  # dnbinom() is -Inf off the integers, so the polish would reject every Newton
  # step, leave alpha at fitNB's value and return the dispersion optimiser's
  # UPPER BOUND for every gene -- measured psi 999.96 -- while reporting success.
  spe <- buildNiches(.toySPE(), sigma = 20)
  f <- fitSpiDE(spe, "condition", sigma = 20, random = "none", verbose = FALSE)
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  SummarizedExperiment::assay(spe, "counts") <- Y + 0.5
  expect_error(polishSpiDE(f, spe, verbose = FALSE), "integer counts")
  # and the message must name the way out
  err <- tryCatch(polishSpiDE(f, spe, verbose = FALSE), error = conditionMessage)
  expect_match(err, "polish")
})

test_that("a dispersion optimum on its search bound keeps fitNB's estimate", {
  # the profile rule is the one with a search; the default keeps fitNB's value
  set.seed(4)
  n <- 300
  W <- cbind(1, scale(rnorm(n)))
  sv <- spiDE:::.newtonSolver(W, c(0, 0))
  # an all-zero gene has no information about overdispersion: the profile
  # likelihood is monotone and optimize() returns its ceiling (~976), which is
  # not an estimate. On the fixed-effects path psi scales the SE directly.
  r <- suppressWarnings(
    spiDE:::.polishGene(rep(0, n), W, c(0, 0), 0.7, c(0, 0), sv, psi.method = "profile"))
  expect_true(r$psi_bound)
  expect_equal(r$psi, 0.7)
  # a well-identified gene is unaffected
  y <- rnbinom(n, mu = exp(W %*% c(1.5, 0.3)), size = 1 / 0.4)
  r2 <- spiDE:::.polishGene(y, W, c(0, 0), 1, c(0, 0), sv, psi.method = "profile")
  expect_false(r2$psi_bound)
  expect_true(r2$polished)
  expect_gt(r2$psi, 0.1)
})

test_that("an unusable Newton keeps fitNB's fit rather than the sane start", {
  # Returning the sane start would zero every tested coefficient, which
  # inference reports as t = 0 with a finite SE: a confident null for a gene
  # that was never converged. A singular solver must fall back to the input.
  set.seed(5)
  n <- 120
  W <- cbind(1, scale(rnorm(n)))
  y <- rnbinom(n, mu = 4, size = 2)
  a0 <- c(0.9, 0.25)
  bad <- list(solve = function(w, s) NULL, xcov = function(w) NULL)
  r <- spiDE:::.polishGene(y, W, a0, 0.6, c(0, 0), bad)
  expect_identical(r$alpha, a0)
  expect_equal(r$psi, 0.6)
  expect_true(r$singular)
  expect_false(r$polished)
})

test_that("a polished fixed-effects fit is scaled by the Pearson dispersion", {
  # psi is no longer edgeR's moderated value, so using it as the SE scale
  # rescaled every gene by the sqrt of a noisy estimate: measured sd(t)
  # 2.11 -> 2.45 on this fixture before the fix.
  spe <- buildNiches(.toySPE(), sigma = 30)
  a <- spiDE(spe, "condition", sigma = 30, random = "none", polish = FALSE,
             fdr = 1, verbose = FALSE)
  b <- spiDE(spe, "condition", sigma = 30, random = "none", polish = TRUE,
             fdr = 1, verbose = FALSE)
  expect_lt(sd(fits(b)[[1]]@t_stat), sd(fits(a)[[1]]@t_stat))
})

test_that("the per-gene diagnostics are keyed by gene, not by position", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  f <- fits(polishSpiDE(fitSpiDE(spe, "condition", sigma = 20, random = "none",
                     verbose = FALSE), spe, verbose = FALSE))[[1]]
  expect_identical(rownames(f@polish), rownames(f@alpha))
  expect_false(anyNA(f@polish["G1", ]))
})

test_that("a wrongly sized lambda.a is refused with a message naming re.celltype", {
  d <- toy_design()
  A0 <- matrix(0, 2, ncol(d$W), dimnames = list(c("a", "b"), colnames(d$W)))
  Y <- matrix(rpois(2 * nrow(d$W), 3), 2, dimnames = list(c("a", "b"), NULL))
  expect_error(
    spiDE:::.polishFit(Y, d$W, A0, c(1, 1), rep(0, ncol(d$W) - 1L)),
    "re.celltype")
})

test_that("an underflowed fitted mean does not become a NaN statistic", {
  # An unwinsorised linear predictor can send exp() to exactly 0 for a few
  # cells. The Pearson working dispersion, which a polished fit now uses on
  # BOTH paths, is then (y - 0)^2 / 0 = NaN, and that propagated to the SE, the
  # statistic and every gene-set test downstream: 18 of 120 statistics were NA
  # on this fixture, and spiGSEA() died with "NAs are not allowed in
  # subscripted assignments".
  spe <- buildNiches(.toyNiche(), sigma = 30)
  r <- spiDE(spe, condition = NULL, sigma = 30, random = "none", fdr = 1,
             verbose = FALSE)
  f <- fits(r)[[1]]
  expect_false(anyNA(f@t_stat))
  expect_false(anyNA(f@se))
  # and the gene-set path that surfaced it runs
  gs <- list(set1 = rownames(f@alpha)[1:5])
  expect_no_error(spiGSEA(r, spe = spe, genesets = gs, type = "niche"))
})

test_that("the absorbed covariance equals the dense one on the TESTED columns", {
  # The Schur complement gives the X-block of the full penalised covariance,
  # which is exactly what the tested columns need -- so inference can avoid
  # building a gram over the ~660 indicator columns. This is the equivalence
  # that substitution rests on, at the sub-block inference actually reads.
  d <- toy_design()
  set.seed(11)
  w <- runif(nrow(d$W), 0.2, 3)
  sv <- spiDE:::.newtonSolver(d$W, d$pen, d$nested)

  info <- crossprod(d$W * sqrt(w))
  diag(info) <- diag(info) + d$pen
  dense <- solve(info)

  xi <- which(!d$nested)
  absorbed <- sv$xcov(w)
  expect_equal(absorbed, dense[xi, xi], tolerance = 1e-8)

  # and the quantities inference reads off it: the diagonal, and a quadratic
  # form in a contrast supported on the tested columns
  k <- seq_len(min(4L, length(xi)))
  expect_equal(diag(absorbed)[k], diag(dense)[xi[k]], tolerance = 1e-8)
  wv <- numeric(length(xi)); wv[k] <- c(0.4, 0.3, 0.2, 0.1)[seq_along(k)]
  wfull <- numeric(ncol(d$W)); wfull[xi] <- wv
  expect_equal(as.numeric(crossprod(wv, absorbed %*% wv)),
               as.numeric(crossprod(wfull, dense %*% wfull)), tolerance = 1e-8)
})

# The re-polish after a variance-component step starts from a converged fit at
# a nearby penalty. Measured on the 0.99.19 cohort runs, a re-polish pass that
# treated that start like fitNB's cost half a cold pass (110-168 min against
# 245-395 on the panel): three information matrices and two profile-psi
# searches per gene, and at bandwidth 10 about 100 genes thrown back to the
# sane start on every pass because their converged fit has a fitted log-mean
# below -10. A warm pass is a few damped Newton steps at the held dispersion.

test_that(".polishFit's warm re-polish at a nearby penalty matches the cold one", {
  d <- toy_design()
  ng <- 5
  A0 <- matrix(0, ng, ncol(d$W), dimnames = list(paste0("G", seq_len(ng)),
                                                 colnames(d$W)))
  A0[, 1] <- 0.5
  Y <- t(vapply(seq_len(ng), function(g) {
    mu <- exp(d$W %*% c(1 + 0.2 * g, rep(0.15, 5), rnorm(8, 0, 0.2)))
    rnbinom(nrow(d$W), mu = as.numeric(mu), size = 1 / 0.4)
  }, numeric(nrow(d$W))))
  dimnames(Y) <- list(rownames(A0), NULL)
  re_group <- ifelse(d$nested, "SampleCellTypeInt", NA_character_)
  first <- spiDE:::.polishFit(Y, d$W, A0, rep(0.4, ng), d$pen, re_group)

  # one Schall step's worth of change in the nested penalty
  pen2 <- d$pen
  pen2[d$nested] <- 2 * pen2[d$nested]
  cold <- spiDE:::.polishFit(Y, d$W, first$alpha, first$psi, pen2, re_group)
  warm <- spiDE:::.polishFit(Y, d$W, first$alpha, first$psi, pen2, re_group,
                             warm = TRUE)

  expect_equal(warm$alpha, cold$alpha, tolerance = 1e-3)
  # the dispersion is held: its profile optimum moves at second order in the
  # penalty change
  expect_equal(warm$psi, first$psi)
  expect_lt(max(abs(warm$loglik - cold$loglik) / abs(cold$loglik)), 1e-4)
  expect_true(all(warm$polish$iterations <= 5L))
  expect_true(all(warm$polish$polished))
})

test_that("a warm re-polish keeps a converged fit whose fitted log-mean is below -10", {
  # a gene absent from one cell type: its cell-type intercept runs off toward
  # -Inf under a zero ridge and the converged fit has log mu < -10 there. That
  # is what the cold path's restart check catches in fitNB's output, and it
  # must not throw a converged fit away on a re-polish.
  n <- 240
  ct <- rep(c("A", "B"), length.out = n)
  gidx <- rep(seq_len(8), length.out = n)
  Z <- stats::model.matrix(~ 0 + factor(gidx))
  colnames(Z) <- paste0("SampleCellType", seq_len(8))
  W <- cbind(CellTypeA = as.numeric(ct == "A"), CellTypeB = as.numeric(ct == "B"),
             x1 = rnorm(n), Z)
  nested <- grepl("^SampleCellType", colnames(W))
  pen <- c(0, 0, 0, rep(1.7, 8))
  ct_cols <- c(TRUE, TRUE, rep(FALSE, 9))
  y <- rnbinom(n, mu = exp(1.5), size = 2)
  y[ct == "B"] <- 0L
  solver <- spiDE:::.newtonSolver(W, pen, nested)

  first <- spiDE:::.polishGene(y, W, rep(0, ncol(W)), 0.5, pen, solver,
                               ct_cols = ct_cols)
  expect_lt(min(as.numeric(W %*% first$alpha)), -10)
  # the cold path restarts it -- the trap the warm path must avoid
  cold <- spiDE:::.polishGene(y, W, first$alpha, first$psi, pen, solver,
                              ct_cols = ct_cols)
  expect_false(cold$restarted)   # a converged zero-count intercept is not a degenerate start
  warm <- spiDE:::.polishGene(y, W, first$alpha, first$psi, pen, solver,
                              ct_cols = ct_cols, warm = TRUE)
  expect_false(warm$restarted)
  expect_true(warm$polished)
  # the identified coefficients stay put; the unidentified intercept stays
  # where the converged fit left it rather than at the sane start (-6.9)
  expect_equal(warm$alpha[-2], first$alpha[-2], tolerance = 1e-4)
  expect_lt(warm$alpha[2], -10)
  expect_equal(warm$psi, first$psi)
})

test_that("forked polish workers run BLAS single-threaded and leave the parent alone", {
  skip_if_not_installed("RhpcBLASctl")
  skip_on_os("windows")
  bp <- BiocParallel::MulticoreParam(2, progressbar = FALSE)
  skip_if_not(BiocParallel::bpnworkers(bp) == 2)
  prev <- RhpcBLASctl::blas_get_num_procs()
  withr::defer(RhpcBLASctl::blas_set_num_threads(prev))
  RhpcBLASctl::blas_set_num_threads(2L)
  skip_if_not(RhpcBLASctl::blas_get_num_procs() == 2L)

  got <- unlist(BiocParallel::bplapply(1:2, function(i) {
    spiDE:::.workerBLAS()
    RhpcBLASctl::blas_get_num_procs()
  }, BPPARAM = bp))
  expect_equal(got, c(1L, 1L))
  expect_equal(RhpcBLASctl::blas_get_num_procs(), 2L)

  # and .polishFit applies it whenever it has more than one worker
  d <- toy_design()
  ng <- 4
  A0 <- matrix(0, ng, ncol(d$W), dimnames = list(paste0("G", seq_len(ng)),
                                                 colnames(d$W)))
  Y <- t(vapply(seq_len(ng), function(g) {
    mu <- exp(d$W %*% c(1 + 0.2 * g, rep(0.15, 5), rnorm(8, 0, 0.2)))
    rnbinom(nrow(d$W), mu = as.numeric(mu), size = 1 / 0.4)
  }, numeric(nrow(d$W))))
  dimnames(Y) <- list(rownames(A0), NULL)
  re_group <- ifelse(d$nested, "SampleCellTypeInt", NA_character_)
  msgs <- capture_messages(
    spiDE:::.polishFit(Y, d$W, A0, rep(0.4, ng), d$pen, re_group,
                       BPPARAM = bp, verbose = TRUE))
  expect_match(paste(msgs, collapse = " "), "one BLAS thread per worker")
  expect_equal(RhpcBLASctl::blas_get_num_procs(), 2L)
})
