# The two switches measured against the shipped profile-ML polish + Pearson
# scale: a moderated dispersion kept at the converged mean, and the
# quasi-likelihood scale on top. The QL machinery's oracle tests are in
# test-quasilik.R.

test_that("polish.psi = 'moderated' converges the mean and keeps fitNB's dispersion", {
  spe <- spiDE:::.toySPE()
  spe <- buildNiches(spe, sigma = 30, verbose = FALSE)
  fm <- fits(fitSpiDE(spe, "condition", sigma = 30, polish.psi = "moderated", verbose = FALSE))[[1]]
  expect_equal(unname(fm@psi), unname(fm@polish$psi_fitnb))
  expect_true(all(fm@polish$polished | fm@polish$iterations == 0L))
  fp <- fits(fitSpiDE(spe, "condition", sigma = 30, verbose = FALSE))[[1]]
  expect_false(isTRUE(all.equal(unname(fp@psi), unname(fm@psi))))
})

test_that("dispersion = 'ql' scores every gene with a finite statistic close to the Pearson one", {
  spe <- spiDE:::.toySPE()
  spe <- buildNiches(spe, sigma = 30, verbose = FALSE)
  fit <- fitSpiDE(spe, "condition", sigma = 30, verbose = FALSE)
  rp <- fits(testSpiDE(fit, spe = spe, fdr = 1))[[1]]
  rq <- fits(testSpiDE(fit, spe = spe, fdr = 1, dispersion = "ql"))[[1]]
  expect_true(all(is.finite(rq@t_stat)))
  ratio <- rq@t_stat / rp@t_stat
  expect_true(all(abs(log(ratio[is.finite(ratio) & rp@t_stat != 0])) < log(3)))
  # the QL scale needs a Pearson-scaled fit: a fixed-effects, unconverged fit
  # has none, and the guard must say so rather than silently fall back
  f0 <- fitSpiDE(spe, "condition", sigma = 30, random = "none", converge = FALSE,
                 verbose = FALSE)
  expect_error(testSpiDE(f0, spe = spe, dispersion = "ql"), "mixed or converged")
})
