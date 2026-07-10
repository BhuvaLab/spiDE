# Builds data/toySpiDE.rda, the small example SpatialExperiment shipped with
# spiDE. Run with: source("data-raw/make_toySpiDE.R")
devtools::load_all()

# a seeded synthetic SpatialExperiment with a planted neighbourhood signal:
# gene G1 is up-regulated in the index cell type "A" in Responders in
# proportion to the local density of the niche cell type "B".
toySpiDE <- .toySPE(n_samples = 6, n_per = 80, n_genes = 20, seed = 1)

usethis::use_data(toySpiDE, overwrite = TRUE)
