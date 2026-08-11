# spiDE 0.99.14

## New Features

* `twoStageSpiDE()`: a second, structurally different estimator of
  condition-dependent niche effects that removes cell-level
  pseudo-replication *by construction* rather than by correction. Stage 1
  estimates each patient's niche slopes with a joint weighted fit per
  (sample, index cell type) -- anchored on a stored `SpaNorm` fit by default
  (`stage1 = "spanorm"`), with `"ols"` and `"nb"` alternatives -- pools
  cores within patients by precision, and stage 2 contrasts the patient
  slopes with a `limma` moderated t-test weighted by `1/(v + tau2)`
  (DerSimonian-Laird `tau2`, pooled over genes). Returns the familiar tidy
  results schema plus three diagnostics (`r2`, `inclusion` with an
  informative-dropout warning, `tau2`). In the paired benchmark it is the
  best-calibrated method in the study, trading power for that calibration;
  restrict `index`/`niche`/genes to a pre-specified hypothesis. Documented
  in full, with equations, in the model vignette.

## Changes

* The vignette suite is trimmed to the three user-facing documents
  (quickstart, model, calibration): the eight built vignettes totalled
  10.85 MB of HTML, over Bioconductor's 10 MB tarball cap on their own. The
  five validation studies and the benchmark tables they read moved to the
  research repository and are published at
  https://bhuvalab.github.io/spiDE-research/ -- pkgdown's "Statistical
  validation" menu links there.
* The shared palette's two-stage colour is now CVD-safe (`#C51B7D` with a
  triangle marker as a secondary channel; the previous `#E7298A` was
  indistinguishable from the intercept teal under deuteranopia).

# spiDE 0.99.12

## New Features

* `updateObject()` methods for `SpiDEFit` and `SpiDEResults`, filling any slots
  that did not exist when an object was serialised.

  Slots added to a class do not appear in objects pickled before them. Reading
  such an object still works, but anything that triggers validity fails --
  including `initialize()`, which is the documented way to re-combine
  single-bandwidth fits across bandwidths:

  ```
  invalid class "SpiDEResults" object: slots in class definition but not in
  object: "results.celltype", "results.patient"
  ```

  This bit a real analysis: fits stored before the `CellType:condition` result
  layers existed could no longer be combined, and those objects represent hours
  of cluster time. `updateObject()` fills absent slots from the class prototype
  and leaves every present slot untouched; it also descends into the contained
  `SpiDEFit` objects. Verified on a four-bandwidth production result: validity
  goes from failing to passing with all fits and result rows preserved.

# spiDE 0.99.11

## Behaviour change

* `spiGSEA()` now defaults to `test = "competitive"` (was `"self-contained"`).
  The simulation benchmark (`research/`, scenario `gsea`) measured both on
  ground truth and the self-contained form **does not control error on
  correlated data**: with nothing planted, at realistic inter-gene correlation
  and a nominal FDR of 0.05, it called 20.6 of 208 sets per replicate, every
  one of them false -- a realised FDP of 1.00 -- against 0.05 sets for the
  competitive test.

  The cause is not the inter-gene correlation term. Setting `rho = 0` roughly
  quadruples the damage, so that term is doing real work, but the
  self-contained test stays badly anti-conservative with the correct `rho`. It
  assumes the gene-level statistics being averaged have unit spread, which
  holds under the null and fails under signal (measured spread 1.0 null, 1.8 at
  the largest effect tested). The competitive form divides by the observed
  spread and is immune.

  Two further results from the same benchmark. Under a global shift affecting
  most genes, the self-contained test called 96.3% of random sets against 0.7%
  for competitive. And because the per-gene statistics track expression, the
  self-contained test's false calls concentrate in abundant sets: called sets
  sat at the 86th expression percentile against the 51st for random sets --
  closely reproducing what is seen on real data.

  `test = "self-contained"` is retained for comparison with the flat-script
  pipeline it replaces, and is documented as unsuitable for inference.

  Power for the competitive test, at FDR 0.05, is usable from about 25 genes
  per set at moderate effects, and saturates by 50 genes at larger ones.

## New Features

* A gene-set benchmark in `research/`: scenarios `gsea` and `gseacal`, with
  gene-gene correlation induced at the rate observed in real spatial data and
  calibrated against the realised residual correlation rather than assumed.

# spiDE 0.99.10

## Improvements

* `testSpiDE()` now stores that correlation on each `SpiDEFit` (new `rho`
  slot), computed as a by-product of the gene blocks inference already loads --
  so it costs no extra pass over the counts. Read it as `fits(res)[[i]]$rho`,
  like any other slot. It is a useful diagnostic in its own right: it says how
  much residual variation is shared across genes, and so how far a set of `m`
  genes falls short of carrying `m` genes' worth of independent evidence. `spiGSEA()` uses it by default and
  therefore needs no counts pass at all, which makes its `spe` argument
  optional: a fitted result can be shared and queried against many gene-set
  collections without the counts matrix travelling with it.

  Fits serialised before this slot existed still load and still work; they fall
  back to estimating the correlation from the counts, exactly as they did when
  they were written.

## Bug Fixes

* `spiGSEA()` errored on a single gene set. The per-bandwidth `vapply()`
  returned a length-k vector rather than a 1 x k matrix, which the Cauchy
  combination then transposed against the weight matrix. Testing one set is
  ordinary usage.

# spiDE 0.99.9

## New Features

* `spiGSEA()` adds a gene-set layer over a fitted model. The per-gene niche
  tests are individually under-powered on sparse spatial data; averaging a
  statistic over the genes of a pathway recovers power, with the average
  inter-gene correlation carried explicitly so that co-regulation is not
  mistaken for evidence.

  Two nulls are available and they answer different questions.
  `test = "self-contained"` (default) asks whether the set's mean statistic
  differs from zero; `test = "competitive"` asks whether it differs from the
  genes outside the set, as `limma::camera` does. The default is the more
  permissive of the two: under a global shift it will call most sets, correctly
  but uninformatively.

  Works on either result layer via `type = "niche"` (the three-way
  celltype:condition:niche statistics) or `type = "celltype"` (the
  CellType:condition statistics).

  Three details differ from the flat-script gene-set code this replaces, all of
  them corrections:

  - statistics become z before averaging (as `camera` does), because the set
    statistic assumes unit variance while a t with v df has variance v/(v-2) --
    immaterial for the niche coefficients, but not for the CellType:condition
    ones where v is the between-patient S-2;
  - bandwidths combine on two-sided p-values, so a set shifting up at one scale
    and down at another no longer cancels to nothing;
  - the inter-gene correlation is estimated per bandwidth rather than once,
    since each bandwidth is a different design and leaves different residuals.

  The correlation is estimated without ever forming the gene x gene matrix
  (1.4 GB at 13,000 genes): standardised residuals are streamed one gene block
  at a time and accumulated through an identity that is exact, not an
  approximation, in memory linear in the number of cells. That pass carries the
  same `backend` / `BPPARAM` / `block.size` controls as the inference stage, and
  is verified to give identical answers across block sizes, worker counts and
  the CPU and GPU backends.

# spiDE 0.99.8

## New Features

* The niche design now carries a `CellType:condition` block, so a response that
  is **cell-type-specific but niche-independent** has a term of its own instead
  of being forced into the three-way niche interaction. Benchmarked over 2,840
  simulated design points (`research/`): the previous niche-only design produces
  **6.8-9.5x more spurious niche calls** on such genes. Calibration is unchanged
  (null type-I `0.0450` vs `0.0449`).

  This is a trade, not a free win. On truth that really is niche-only, power is
  lower (TPR `0.183` vs `0.319`), because the new block absorbs part of the
  signal the three-way term used to carry alone. `vignettes/spiDE-mixed-benchmark.Rmd`
  reports both sides.

  Note the **cell-means coding**: there is no bare `condition` main effect. The
  condition contrast is carried by the `CellType:condition` (`ResponseCellType`)
  columns, one per cell type. Code matching a single `"Response"` column finds
  nothing and should match `ResponseCellType`.

* `results()` gains `type = "celltype"` and `type = "patient"` alongside the
  default `"niche"`: cell-type-specific response calls keyed by
  `(gene, ct_index)`, and one abundance-weighted response contrast per gene.

## Improvements

* **All three result layers now combine evidence across every bandwidth** with
  the log-likelihood-weighted Cauchy combination. The cell-type and patient
  layers previously reported the last bandwidth's fit alone. On a four-bandwidth
  production fit the cell-type layer returns 107 genes combined against 93 for
  the best single bandwidth, and 0 under the old behaviour -- the widest
  bandwidth was the least informative of the four. If you compared these layers
  across analyses with different `sigma` sets, those results change.

* Cauchy combination now consumes **two-sided** p-values. Under one-sided input
  `tan((0.5 - p)pi)` diverges to `-Inf` as `p -> 1`, so a gene up in one niche
  and down in another cancelled exactly. Brown's method stays one-sided, where
  `-2log(p)` is bounded at 0 and safe.

* Consistent FDR scale across the three layers: `.dirBH()` tested each direction
  at `q/2` but returned the per-direction `q` unscaled, understating cell-type
  and patient q-values roughly two-fold relative to the niche layer.

* The toy data generators no longer reseed the caller's global RNG. They take a
  `seed` argument for reproducibility but previously called `set.seed()`
  directly, perturbing every random draw the caller made afterwards.

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
