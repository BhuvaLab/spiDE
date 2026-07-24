# GPU test helpers, mirroring SpaNorm's tests/testthat/helper-gpu.R conventions.
# spiDE's GPU inference path is built on SpaNorm's exported device/tensor layer,
# so the availability probe and tolerances defer to SpaNorm's own resolution.

skip_if_no_gpu <- function() {
  testthat::skip_if_not(SpaNorm::checkGPU(), "torch GPU/MPS not available")
}

# accelerator tolerance: near machine precision on float64 devices (cuda/cpu),
# looser on MPS which can only run float32
gpu_tol <- function() {
  if (SpaNorm::getBackendDevice() == "mps") 1e-4 else 1e-8
}
