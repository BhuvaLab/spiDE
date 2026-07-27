# spiDE 0.99.4

* Changed `fitSpiDE()`'s default `re.prop` from `0.2` to `1` (no cell
  subsampling in the random-effect variance-component loop). A replicate
  study on real data (`vignettes/spiDE-mixed-benchmark.Rmd`) found that
  `re.prop < 1` doesn't just add noise to the fitted `tau2` — the noise
  itself stays roughly flat from `re.prop = 0.2` to `0.8` (never shrinking
  below the genuine between-patient signal it's confounded with), and the
  **mean** `tau2` is systematically biased downward at every `re.prop < 1`
  tested, an attenuation that more replicates cannot average away. GPU
  backends make `re.prop = 1` affordable in practice (see the benchmark
  vignette), so it is now the default; lowering it remains possible but is
  rarely advised (see `?fitSpiDE`).
* Bumped the `SpaNorm` requirement to `>= 1.7.7`, which fixes
  `getGPUMemoryBudget()` reporting a MIG-partitioned GPU's whole physical
  card instead of the process's assigned slice — spiDE's GPU inference path
  (`R/inference.R`) calls this directly, so on MIG hardware the bug could
  feed a many-fold-too-large budget into `.covBatchSize()` and cause an
  out-of-memory failure instead of the intended blocked, bounded-memory fit.
* Fixed a memory blowup that made `combine = "cauchy"` (the default) unusable
  with random effects on realistically-sized data, on **both** backends. The
  batched Wald covariance introduced in 0.99.3 precomputed a Khatri-Rao cross
  term of the design (`ncells x p^2`) once per bandwidth. Being built outside
  the gene loop, its size could not be bounded by `block.size`, and it scales
  quadratically in the design width: a 602-column random-intercept design over
  21,843 cells needs 63 GB, and a 4,906-column random-slope design 4.2 TB.
  The per-gene Gram matrices are now built with a batched matmul over a
  bounded sub-batch of genes instead — memory is linear in `p` and capped
  automatically (`options(spiDE.cov.mem.budget = <bytes>)` on the CPU path,
  the GPU budget otherwise), independently of `block.size`. Results are
  unchanged (bit-identical on the CPU path) and invariant to the sub-batch
  size.

# spiDE 0.99.3

* `testSpiDE()` / `spiDE()` gain a `backend` argument (`"auto"`, `"cpu"`,
  `"gpu"`) that GPU-accelerates the inference stage: the per-gene Wald
  covariance and negative-binomial working weights are batched across each
  gene-block and computed on the accelerator via SpaNorm's tensor engine,
  rather than one gene at a time. On the GPU backend the block size is
  auto-selected to keep peak device memory within the detected budget
  (override via the new `gpu.mem.budget`), and a serial `BPPARAM` is used
  automatically to avoid multiple processes contending for one device. The
  batched Wald covariance also speeds up the CPU path. GPU results match the
  CPU path to a small tolerance (single precision on Metal/MPS); the CPU path
  stays exactly reproducible and block-size invariant. Requires
  `SpaNorm (>= 1.7.6)`.

# spiDE 0.99.0

* Initial development version.
* `buildNiches()`, `mergeNiches()`, `computeSizeFactors()` for constructing
  niche covariates and per-cell-type size factors.
* `fitSpiDE()` fits a per-gene negative binomial GLM over the neighbourhood
  interaction design using the SpaNorm `fitNB()` engine.
* `testSpiDE()` combines Wald statistics across correlated covariates (Brown's
  method) and bandwidths (Cauchy combination), with hierarchical FDR control.
* `spiDE()` convenience wrapper for the full workflow.
