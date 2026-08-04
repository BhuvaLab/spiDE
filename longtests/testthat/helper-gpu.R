# longtests/ runs as its own testthat directory, so it does not see
# tests/testthat/helper-gpu.R. The GPU tolerance is duplicated here rather than
# sourced across directories, because the long-tests tarball is built from this
# directory alone and a cross-directory source() would work locally and fail on
# the builder.

gpu_tol <- function() {
  if (SpaNorm::getBackendDevice() == "mps") 1e-4 else 1e-8
}
