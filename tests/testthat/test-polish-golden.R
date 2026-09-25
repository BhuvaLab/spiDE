# Safety net for moving the polish machinery to SpaNorm
# (design/plans/2026-09-25-polish-to-spanorm.md): pins the production
# fit -> polish -> test outputs at tolerance 0 against a golden snapshot
# captured with the pre-move code. The `meta` element records the SpaNorm/
# spiDE versions and seed the golden was captured with and is not compared;
# after the move, the golden must still reproduce exactly under the SpaNorm
# version recorded there.
test_that("polish + inference reproduce the pre-move golden outputs exactly", {
  skip_on_cran()
  g <- readRDS(test_path("_golden", "golden_polish.rds"))
  now <- golden_polish_fits()
  for (nm in setdiff(names(g), "meta")) {
    expect_identical(names(now[[nm]]$fits), names(g[[nm]]$fits), info = nm)
    for (bw in names(g[[nm]]$fits)) {
      for (s in names(g[[nm]]$fits[[bw]])) {
        expect_equal(now[[nm]]$fits[[bw]][[s]], g[[nm]]$fits[[bw]][[s]],
                     tolerance = 0, info = paste(nm, bw, s))
      }
    }
    expect_equal(now[[nm]]$results, g[[nm]]$results, tolerance = 0, info = nm)
  }
})
