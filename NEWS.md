# spiDE 0.99.7

* **Behaviour change:** `fitSpiDE()` / `spiDE()` now default to
  `df.method = "satterthwaite"` (previously `"between"`). Consequently a
  mixed-effects fit's `SpiDEFit@df` is a **named per-tested-coefficient vector**
  by default, where it used to be the scalar `S - 2`; code that assumed a scalar
  (e.g. `stats::pt(t, fit@df)` across genes) must now index the column it is
  testing, or pass `df.method = "between"` explicitly to restore the old
  behaviour. Nothing else about the fit changes — only the reference df used at
  inference time.

  The change follows the completed benchmark study (`research/`), which measured
  both arms on identically seeded data. `"between"` is severely over-conservative
  when samples are few — null type-I error `~0.001` at `S = 4` against a nominal
  `0.05`, with correspondingly near-zero power — and only approaches nominal at
  the largest sample sizes studied. `"satterthwaite"` holds type-I error within
  `0.042`–`0.065` across the whole sampled range and gains `~0.10` mean TPR
  (paired, FDR 0.05), with observed FDP no worse. The trade is a mild liberal
  drift at larger `S` (worst measured `~0.065` against nominal `0.05`); a
  Kenward–Roger bias correction of the variance-parameter covariance is the
  indicated next step for closing it. Use `df.method = "between"` where strict
  conservatism matters more than power, or for back-compatibility.

# spiDE 0.99.6

* `fitSpiDE()` / `spiDE()` gain a `df.method` argument for the mixed-effects fit
  (`random != "none"`). The default `"between"` keeps the single between-patient
  `S - 2` reference df (back-compatible). The new `"satterthwaite"` computes a
  **per-coefficient** Satterthwaite reference df, so the response-niche
  interactions — which carry within-patient information and were over-conservative
  under a flat `S - 2` — get their larger effective df, while the response main
  effect stays at `~ S - 2` by construction. `SpiDEFit@df` is accordingly a scalar
  under `"between"` or a named per-tested-coefficient vector under
  `"satterthwaite"`. The df is invariant to the per-gene dispersion, so it is
  computed once per bandwidth and adds nothing to the per-gene inference loop. See
  the *spiDE model* and *Speeding up the mixed-effects fit* vignettes.

# spiDE 0.99.5

* `mergeNiches()` now records the merged niche group membership on the object,
  and the design matrix built by `fitSpiDE()` / `nicheDesign()` uses it to drop
  every covariate whose index cell type is a member of the merged niche it is
  tested against (previously only the exact index-equals-niche self interaction
  was dropped). This matches the neighbourhood self-density exclusion of the
  original analysis scripts.

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
