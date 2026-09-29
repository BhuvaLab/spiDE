# Builds data/toySpiDE.rda, the small example SpatialExperiment shipped with
# spiDE. Run with: source("data-raw/make_toySpiDE.R")
devtools::load_all()

# a seeded synthetic SpatialExperiment with a planted neighbourhood signal:
# gene G1 is up-regulated in the index cell type "A" in Responders in
# proportion to the local density of the niche cell type "B". Sixteen patients
# (eight per condition): spiDE's tests take patients as their units, and the
# condition-specific test needs at least three patients per condition with
# usable slopes (0.99.30; 6 patients before).
toySpiDE <- .toySPE(n_samples = 16, n_per = 150, n_genes = 20, seed = 7)

usethis::use_data(toySpiDE, overwrite = TRUE)
