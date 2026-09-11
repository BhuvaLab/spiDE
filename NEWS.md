# spiDE 0.99.20

## Changes

* **The polish stage costs one cold pass, not four.** On the 0.99.19 cohort
  runs (769 genes, four workers) the stage took 620-700 min: a cold pass of
  245-395 min and three re-polish passes of 110-168 min each, and four of
  fifteen tasks hit their 12 h limit inside the fourth pass. Three causes,
  each measured, each fixed:
  * **Forked polish workers now run BLAS single-threaded.** A forked worker
    inherits the parent's OpenBLAS thread count, so four workers on eight
    cores ran 32 threads. At the cohort's design shape, four workers with
    four threads each on four cores take 2.3-2.6 s per Newton step against
    0.24-0.28 s with one thread each (9x), while a single worker gains only
    1.5x from four threads. Both blocked stages, the polish and the inference (whose per-gene
    gram is the same BLAS work), now dispatch through a wrapper that sets
    one BLAS thread inside each worker when RhpcBLASctl (now in Suggests)
    is installed; the parent process is left alone.
  * **The re-polish after a variance-component step is warm.** It starts
    from a converged fit at a nearby penalty, so it is a few damped Newton
    steps per gene at the held dispersion: no profile-psi search, and no
    log-mean restart check, which had thrown about 100 of 769 genes back to
    the sane start on every re-polish pass at bandwidth 10 (3-5 in the cold
    pass). It agrees with a cold re-polish to 1e-3 on the coefficients; on
    the clustered fixture the variance components agree to four digits and
    the t-statistics to 0.01, and re-profiling the held dispersion at the
    end of the loop would move them by at most 0.014.
  * **The variance-component loop converges instead of stopping at a cap.**
    Schall's update is a linearly convergent fixed-point map, and with the
    cap of three the (sample x cell type) component was still 6-12% above its
    extrapolated limit in thirteen of fifteen cohort runs (the per-sample one
    is within 1% after two steps). That matters: a tenfold error in a
    near-zero nested component moves individual t-statistics by up to 0.13.
    The loop is now Steffensen-accelerated (Aitken extrapolation after every
    two plain steps, per component, only while the steps contract), with
    `tau2.maxit = 10`, `tau2.tol = 1e-2` on `log(tau2)` and
    `tau2.accelerate = TRUE` (`FALSE` is the plain map), and it ends with one
    more warm pass at the reported penalty so coefficients, penalty and df
    are consistent.
* `@polish` gains `repolish.iterations`, the Newton steps spent in the
  re-polish passes, `repolish.capped` and `repolish.singular` (a warm pass
  that hit its cap or a singular system, which leaves the gene at its previous
  converged fit), and `psi_fitnb` is `fitNB`'s dispersion again: the
  0.99.19 loop handed each re-polish the previous pass's polished value, so
  the "polished / fitNB" ratio the cohort driver reports read about 1.
* `.toyClustered(sd_nested = )` plants a (sample x cell type) intercept, so a
  fixture can have an interior nested variance component.
* The polish's two dispersion rules are documented with their cost
  (`polishSpiDE(psi = "moderated")` skips the dispersion search and the two
  re-polishes it triggers) and benchmarked from one fit polished both ways:
  the harness times the polish on its own and the timing table carries one
  row per arm and size (`research/config.R`, arm `polish-rules`).

# spiDE 0.99.19

## Changes

* **The pipeline is fit -> polish -> test -> gsea, and the polish is a stage.**
  `fitSpiDE()` fits the shared model and nothing more: `converge`,
  `converge.maxit`, `converge.tol` and `polish.psi` are gone from it.
  `polishSpiDE()` is the stage the user runs, or skips, according to their
  data; `spiDE()` runs it by default (`polish = TRUE`) and `polish = FALSE`
  tests the shared fit as `fitNB` returned it.
* **The variance components are re-estimated from the converged fit.** The
  fit's Schall loop reads the shared fit's own coefficients and dispersion,
  which can sit far from every gene's optimum: on the clustered test fixture
  it reported a between-sample variance of 10 against a planted 0.49 (from
  sample intercepts three times too wide), under every release since
  2026-08-08. `polishSpiDE(tau2 = TRUE)` (the default) takes one Schall step
  on the polished coefficients with the gene-averaged weights at the
  polished mean and dispersion, re-polishes at the new penalty, iterates to
  `tau2.tol`, and refreshes the Satterthwaite reference df. The fixture's
  component lands in the planted window. Calibration had survived the old
  value because the QL scale is robust to the dispersion and an inflated
  component only weakens the ridge, but the reference df, the shrinkage and
  the Newton weights all read it.
* **The polish's dispersion rule is the profile value again**
  (`polishSpiDE(psi = "profile")`). The moderated rule, made the default in
  0.99.18 for a third of the cost, keeps whatever the shared fit left, which
  on the same fixture is fifteen times the converged value; the
  variance-component step needs a dispersion consistent with the converged
  mean. The two rules remain indistinguishable on the null and in the real
  cohort's calls, and the cost argument has weakened now that the GPU
  carries the cohort's real grid in a third of the CPU time.
* **The two-stage estimator has left the package.** `twoStageSpiDE()` and its
  stages are archived as a standalone research package, `spiDEtwostage`
  (`research/twostage/` in the research repository), which returns its own
  result object; `SpiDEResults` loses the `diagnostics` slot and the validity
  relaxation that allowed an empty `fits`. The mixed-effects estimator is the
  recommended approach: it fits all cells jointly, so thin cell types borrow
  strength, and on the real cohort it was the better-calibrated of the two.
  The two-stage benchmark stays on the research site as the record of that
  comparison.

# spiDE 0.99.18

## Changes

* **The standard errors are scaled by the quasi-likelihood dispersion.**
  `testSpiDE(dispersion = "ql")` is the default: each gene's edgeR-v4
  quasi-likelihood dispersion (its adjusted NB deviance over its effective
  residual df, `SpaNorm::qlDispersion()`) moderated across genes with
  `limma::squeezeVar()`, in place of the working Pearson dispersion, which
  remains available as `dispersion = "pearson"`. Measured on the synthetic
  benchmark (40 replicates, intercept mode, Satterthwaite df) it is the only
  configuration that holds the nominal type-I error at every sample size:
  0.050–0.055 at S = 4, 10, 16 and 30 against 0.091, 0.076, 0.071 and 0.068
  for the Pearson scale. Its recall at FDR 0.05 is that of the design with
  neither 0.99.17 switch (0.360 against 0.361 at S = 30, and 0.402 for the
  Pearson scale, whose extra recall is bought with its inflated null) and its
  realised FDP sits far below nominal (0.014 at S = 30, 0.000 at S = 4 where
  the Pearson scale reaches 0.59). On the real cohort it calls the same
  triplets as before (111 of 114 shared) at no extra cost. It runs on either
  backend: the deviance moments are one shared (log mu, log phi) table per
  gene block and the rest is elementwise over genes x cells on the device.
  A fixed-effects, unconverged fit (`random = "none", converge = FALSE`)
  has no such scale and keeps its legacy NB-dispersion scale with a message.
* **The convergence stage keeps the moderated dispersion.**
  `fitSpiDE(polish.psi = "moderated")` is the default: the coefficients are
  converged per gene under `fitNB`'s cross-gene moderated dispersion, and the
  per-gene profile-ML re-estimate is available as `polish.psi = "profile"`.
  The two are indistinguishable on the null at every sample size (to the
  third decimal) and call the same triplets on the real cohort, and the
  moderated one runs the cohort's real grid in two thirds of the wall time
  (178 against 266 minutes), because the profile step was the expensive part
  of the stage. `polishSpiDE()` takes the same argument with the same default.
* The quasi-likelihood machinery now lives in SpaNorm (>= 1.7.10) as
  `nbUnitDeviance()`, `nbDevianceMoments()` and `qlDispersion()`, with both
  backends and the oracle tests against `edgeR::glmQLFit()`; spiDE carries
  only the wiring. The 2026-08-31 refutation of the QL dispersion is
  withdrawn: it was scored on the composition bias that the nested intercept
  removes, which no per-gene scale could fix.

# spiDE 0.99.17

## New Features

* `re.celltype` (default `TRUE`) adds a nested (sample x cell type) random
  intercept to the niche design alongside the per-sample one, so that every
  tested niche slope is a within-(sample, cell type) slope. Without it the
  slopes also carry the between-sample composition effect -- a patient-level
  association between a cell type's mean niche density and its mean expression
  in that type -- which a shuffle null that permutes within (sample, cell type)
  preserves exactly. On the YTMA cohort that confound was the cause of the
  triplet-level FDR failure: with the nested block the shuffle null is
  calibrated in every expression band (`sd(null t)` 0.96-0.99, no expression
  gradient, no `|t| > 4.89` exceedances against 95 before) where before it ran
  to 1.42 for the brightest genes. `re.celltype = FALSE` reproduces
  pre-correction fits. `random = "slope"` is not a substitute: per-sample
  slopes on the `CellType:niche` bases leave the group means untouched.

* `converge` (default `TRUE`) converges each gene to its own penalised
  negative-binomial optimum after `SpaNorm::fitNB` returns, and re-estimates
  its dispersion by profile maximum likelihood at the converged mean.
  `fitNB`'s multi-gene IRLS shares one cell weight vector across genes and
  stops on the aggregate log-likelihood, which leaves bright,
  cell-type-restricted genes one to four standard errors short of their own
  optimum with a dispersion about 1.6 times too large. The stage is per-gene,
  so it is blocked and parallelised over `block.size` / `BPPARAM` as inference
  is, with the nested indicator block absorbed by a Schur complement so the
  per-iteration cost does not grow with the number of groups. New
  `SpiDEFit@polish` records the per-gene diagnostics. Note this replaces
  edgeR's cross-gene moderated dispersion with a per-gene one, and that it
  *raises* the null variance of bright genes -- their previous standard errors
  were inflated by a dispersion estimated off the optimum.

  It also **sharpens real signal substantially**, which was not the reason for
  building it. On the toy fixture's planted niche effect the test statistic
  goes from 1.68 (FDR 0.09, not called) to 10.19 (FDR 2e-15) in condition mode,
  because the unconverged fit put the planted gene's dispersion at 3.09 where
  its own optimum is 0.28 and the inflated standard error buried the effect. In
  niche mode a SPURIOUS competing association that outranked the true one
  (|t| 7.82 against 5.63) disappears, and the true one goes to 14.27. Two unit
  tests had encoded those artefacts and were updated.

## Changes

* `nicheDesign()` gains `re.celltype`, defaulting to `FALSE` -- the same
  deliberate asymmetry as `random`, since a design returned with
  penalty-identified columns is rank-deficient.
* `.toySPE()` gains `composition`, which plants a between-sample composition
  confound with zero within-sample niche slope.

* `polishSpiDE()`: the per-gene convergence stage as a **post-hoc
  adjustment** on an existing fit — a fit made with `converge = FALSE`, one
  serialised by an older version, or one to re-converge with other settings —
  without refitting. It shares one implementation and one penalty rule with
  `fitSpiDE(converge = TRUE)`, so the two routes give the same fit (tested to
  1e-8), and it clears every inference slot derived from the old
  coefficients so `testSpiDE()` recomputes them.
* `compositionTest()`: the between-sample composition association that the
  nested intercept absorbs, tested on its own terms at the patient level —
  pseudobulk per (sample, index type), the sample's mean niche density around
  those cells, a `limma` moderated *t* across samples, with `"niche"` (pooled)
  and `"condition:niche"` (the patient-level counterpart of the three-way
  term) reported per (gene, index, niche). It is real signal, and it is not
  neighbourhood-dependent DE; the two now have two tests.
* `.toySPE(composition = )` now plants a confound that reaches the tested
  slope: each sample's A cells are shifted toward or away from the B-rich
  region, so samples differ in *mean* niche density the way tissue does, and
  G2's baseline in Responders' A cells shifts with that mean. The earlier
  version tied B-cell prevalence to condition, which moved the between-sample
  mean too little to matter next to the within-sample spread.

## Documentation

* The vignettes are reorganised into four, written top-down rather than in
  the order the work happened: the quickstart (calls and result layers, no
  justification), *The spiDE model* (the model stated once, organised around
  the Frisch–Waugh–Lovell theorem, with the per-gene convergence — damped
  Newton, profile-ML dispersion, the Schur-complement absorption of the
  nested block — and the combination and FDR steps explained for the first
  time), *The two-stage estimator* (split out), and *Calibration: reading
  lambda, and what the benchmarks measured*, now the only place in the
  package that quotes benchmark numbers. Research history has been removed
  from all four.

## Bug Fixes

* **Non-integer counts are refused when `converge = TRUE`.** `dnbinom()` is
  `-Inf` off the integers, so the convergence stage rejected every step, left
  the coefficients at `fitNB`'s values and returned the dispersion optimiser's
  *upper bound* for every gene (measured 999.96) while its own diagnostics
  reported success. `checkCounts()` now refuses such an assay before the fit.
* **The convergence stage no longer densifies the whole counts matrix.** Its
  automatic block size only applied with more than one worker, so under the
  default `SerialParam()` a single block realised every gene at once (8.3 GB on
  a real cohort), breaking the package's never-densify invariant.
* **A polished fit is scaled by the Pearson working dispersion on both paths.**
  The fixed-effects path scaled its standard errors by `psi`, which is safe only
  while `psi` is edgeR's moderated value; after convergence it is not. Measured
  on the toy fixture, converging took `sd(t)` from 2.11 to 2.45 before this fix
  and to 1.92 after it.
* **Fit and inference now use one mean function.** The convergence stage
  maximises the unclamped likelihood, so inference evaluates it unclamped rather
  than at `calculateMu()`'s winsorised mean; the planted toy effect's statistic
  goes 3.98 to 5.78. With `converge = TRUE`, `winsor` therefore sets only the
  Newton starting point.
* **The Satterthwaite fallback no longer aborts the fit.** Its documented
  degrade-to-`"between"` path returned a bare number from `.fitNBmixed()`, so
  the caller died with `$ operator is invalid for atomic vectors`.
* A dispersion optimum sitting on its search bound, and a gene whose Newton
  cannot run, both now keep `fitNB`'s estimate and are flagged in `@polish`
  (`psi_bound`, `polished`) rather than storing a bound or a zeroed fit.
* `@polish` is keyed by gene rather than by position, and a wrongly sized
  `lambda.a` is refused with a message naming `re.celltype`.
* **An underflowed fitted mean no longer becomes a `NaN` statistic.** An
  unwinsorised linear predictor can send `exp()` to exactly 0 for a few cells,
  and the Pearson working dispersion is then `0/0`; it propagated to the
  standard error, the statistic and `spiGSEA()`, which failed outright. The
  fitted mean is floored at `exp(-30)` in the convergence stage and at
  inference, at the same value.
* **A covariate with missing or non-finite values is refused by name.**
  `model.matrix()` drops those rows, so the design no longer matched the
  random-effect block and the run died in `cbind()` with `number of rows of
  matrices must match`. `log()` of a QC column that is zero for some cells is
  the usual way in.
* **Inference absorbs the nested indicator block.** The (sample x cell type)
  columns are a 0/1 partition of the cells, so the per-gene covariance comes
  from a Schur complement over the dense columns alone rather than a gram over
  the full design. Measured at the real cohort shape (77,454 cells, 345 dense +
  660 nested columns): **2.16 s -> 0.35 s per gene**, or 8.0 h -> 1.3 h for a
  13,348-gene transcriptome, agreeing with the dense path to 7e-21. The saving
  is size-dependent -- at a smaller shape the two are within 20% -- so it is
  kept unconditional only because both are cheap there.
* Examples now use a single bandwidth, and the fast fixed-effects path where
  their subject is an accessor rather than the mixed correction, taking
  `R CMD check`'s example time from ~23 minutes to ~10.

* `updateObject()` could not repair an object serialised before a slot whose
  prototype is `NULL` was added (`polish`, and latently `re_group`, `tau2`,
  `penalty`, `df`). `.fillSlots()` used `attr(object, s) <- value`, and
  `attr(x, "s") <- NULL` *removes* an attribute rather than setting it, so the
  slot stayed absent: the object remained invalid, `show()` errored, and
  `updateObject()` -- the documented repair path -- failed on exactly the
  objects it exists to repair. It now assigns through `methods::slot<-`.
* The per-gene convergence stage ran as a single block regardless of
  `BPPARAM`, so a caller requesting several workers got no parallelism. Absent
  an explicit `block.size` it now splits one block per worker; blocking is
  exact, so the result is unchanged.

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
