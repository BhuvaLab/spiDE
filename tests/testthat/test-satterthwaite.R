test_that(".ptByCol handles NULL / scalar / vector / matrix df", {
  tm <- matrix(c(0, 1, -1, 2), nrow = 2)         # 2 genes x 2 cols
  # NULL -> normal
  expect_equal(spiDE:::.ptByCol(tm, NULL), stats::pnorm(tm))
  # scalar -> same df everywhere
  expect_equal(spiDE:::.ptByCol(tm, 5), stats::pt(tm, df = 5))
  # per-column vector -> col 1 uses df=5, col 2 uses df=50
  ref <- cbind(stats::pt(tm[, 1], df = 5), stats::pt(tm[, 2], df = 50))
  expect_equal(spiDE:::.ptByCol(tm, c(5, 50)), ref)
  # single-gene vector input with per-column df
  expect_equal(spiDE:::.ptByCol(c(0, 2), c(5, 50)),
               c(stats::pt(0, 5), stats::pt(2, 50)))
  # lower.tail forwarded
  expect_equal(spiDE:::.ptByCol(tm, 5, lower.tail = FALSE),
               stats::pt(tm, df = 5, lower.tail = FALSE))
})
