# Vignette reorganisation — plan for review

Date: 2026-09-07. Status: **awaiting approval**; nothing has been changed yet.

## 1. Diagnosis of the current three vignettes

Total 2,123 lines: quickstart 299, model 1,532, calibration 292.

**The model vignette is in discovery order, not logical order.** Its overview
promises four parts (model, pseudo-replication, mixed-effects correction,
two-stage), but the body has grown to roughly fifteen sections, in the order
the work happened:

1. "The full spiDE model" states a model with a per-sample library-size term,
   a global response effect `eta_g` and *no random effects and no nested block*.
2. "Mixed-effects correction" then presents "The model" a second time, with a
   sample random intercept and optional random slopes — still no nested
   (sample × cell type) block, which is the default since 0.99.17.
3. The nested block is introduced ~600 lines later, in "Why the niche slope
   must be a within-(sample, cell type) slope", as a third model statement.
4. "Why the design carries `CellType:condition`" sits between the two, added
   this week, and is written as a benchmark decomposition rather than as part
   of the model.
5. Research history that users do not need: the "nested fixed model built,
   measured and rejected" subsection with two tables from a refuted
   construction; the imbalance sweep; "What the study found" for the df
   arms; the two-stage section's "how it measured", withdrawn numbers and the
   superseded-sweep narrative.
6. The two-stage estimator takes 470 lines (a third of the file) for a method
   the same vignette says not to use by default.
7. Never explained anywhere in the model vignette: the within-gene combination
   of correlated niche p-values (Cauchy / Brown), the cross-bandwidth
   combination, and the hierarchical FDR cascade — the three steps that turn
   Wald statistics into the results table. The quickstart's introduction has
   one sentence on them.

**Redundancy across the three files.** The composition confound is explained
in the quickstart twice ("The full workflow in one call" and "Mixed-effects
inference"), plus a third section ("The between-sample composition
association"), in the model vignette, and in the calibration vignette twice
("What drives lambda", "When lambda exceeds 1"). The 0.99.17 benchmark numbers
(0.071 vs 0.058, 0.36 vs 0.17, 0.7× SE) appear in all three vignettes. The
two-stage estimator is described in the quickstart, the model vignette and
the calibration vignette's cross-references.

**The quickstart is no longer a quickstart.** A user's first page now carries
justification paragraphs, benchmark numbers and the estimand argument before
the second code chunk.

**The calibration vignette is coherent** but has acquired design arguments
("do not reach for a smaller model either", the composition driver) that
belong with the model.

## 2. Proposed set: four vignettes

Split the two-stage estimator out; give the evidence one home; make the
model vignette a single top-down statement of the framework.

### A. `spiDE.Rmd` — Quickstart (target ~180 lines, from 299)

Task-oriented only. Order:

1. What spiDE tests, in three sentences, and the three stages.
2. Example data.
3. One call, `results()`, and **reading the result layers**: `"niche"`,
   `"celltype"`, `"patient"`, plus `compositionTest()` as the fourth question
   (patient-level composition) — each defined by *what quantity it estimates*
   in one line, no numbers, no history.
4. Step by step (build / fit / test), with `polishSpiDE()` as the post-hoc
   route.
5. Niche-only mode (`condition = NULL`) — kept, trimmed to what it does.
6. Scaling: blocks, `BPPARAM`, GPU.
7. "The defaults, in one box": `random = "intercept"`, `re.celltype = TRUE`,
   `converge = TRUE`, `df.method = "satterthwaite"`, `combine = "cauchy"` —
   one line each saying what it does, with a single pointer to the model
   vignette. No justification prose.

Removed from here: all composition-confound paragraphs, all benchmark
numbers, the two-stage discussion (one sentence + pointer remains), the
"reading the fit's diagnostics" subsection (moves to the model vignette's
slot table).

### B. `spiDE-model.Rmd` — The model (target ~750 lines, from 1,532)

One narrative, top-down, the model stated **once** in full. The organising
idea, stated early and used throughout, is the **Frisch–Waugh–Lovell
theorem**: every nuisance block in the design partials one component out of
the tested slope, and the design is read as a sequence of such blocks.

1. **The question and the estimand.** Within an index cell type, how does the
   condition's effect on a gene change along the local density of a niche
   cell type? Name the three quantities the design separates: the flat
   cell-type response (`CellType:condition`), the condition-independent niche
   slope (`CellType:niche`), and the tested three-way slope. One paragraph on
   what a "niche-dependent" call means and does not mean.
2. **The effective niche** (current text, unchanged in substance).
3. **The model, stated once.** One equation with every block: within-sample
   fixed terms, the tested response terms, the sample random intercept
   `u_s`, the nested `(sample × cell type)` intercept `c_{sk}`, optional
   random niche slopes. Coefficient table. Cell-means coding. Variance
   function. Then the FWL section:

   *Frisch–Waugh–Lovell, and what each block does to the tested slope.* State
   the theorem: in a (weighted) least-squares fit of `y` on `[X1 | X2]`, the
   coefficient on `X2` equals the coefficient from regressing `y` on `X2`
   after both have been residualised on `X1`. Then apply it three times.
   (i) `X1` = the `(sample × cell type)` indicators: residualising is
   centring within group, so the tested slope is estimated from
   within-(sample, cell type) deviations of niche density and expression, and
   the between-sample covariance — patients whose index cells sit in denser
   niches also expressing the gene differently in that type — is removed
   exactly. Without the block it loads onto the slope and is reported with a
   cell-level standard error. (ii) `X1` = `CellType:condition`: the flat
   response is partialled out, so the slope is the gradient and not the
   mean; without it the slope is a regression through the origin whose
   standard error is too small. (iii) IRLS: each Newton/IRLS step is a
   weighted least-squares problem, so FWL applies step by step and hence at
   the converged fit. Then the one caveat: the block is ridge-penalised, so
   the centring is *partial* — the group means are shrunk toward the pooled
   mean by `tau2` — and exact only as `tau2 → ∞`; on the real cohort the
   shrinkage reached the exact-centring null within measurement (one sentence,
   pointer to the calibration vignette for the numbers).
4. **Why cells are not replicates.** The pseudo-replication demonstration:
   the 300-gene null simulation, fixed vs mixed rejection rate, the QQ figure.
   Kept, shortened by removing the "why not sample as a fixed effect" toy and
   its collinearity walk-through to one paragraph (the `NA` coefficient
   demonstration stays as one chunk; the rest is prose). **Cut:** the nesting
   aside, the imbalance table, and the "nested fixed model built, measured and
   rejected" subsection with both tables. They move to nowhere in the package;
   the research site already carries `between-sample-stratum.html`, and
   CLAUDE.md records the outcome.
5. **Random effects as ridge.** The identity, the Bayesian and shrinkage
   readings, the BLUP — condensed to about half its current length. Keep the
   `glmmPQL` cross-check, the convergence trace figure and the shrinkage-arrows
   figure: they are the vignette's best teaching material. Drop the Gaussian
   density check figure (it makes a point already made by the trace).
6. **Fitting.** (a) The outer Schall loop for `tau2` (moment update, effective
   df, `re.maxit`, subsampling note) — condensed. (b) **Per-gene convergence**,
   new and explicit:
   - *Why:* the shared-weight IRLS is not at each gene's own optimum for
     bright, cell-type-restricted genes (one sentence with the magnitude).
   - *Damped Newton on the gene's own penalised log-likelihood:* the objective,
     the score and information, step-halving, the sane-start restart.
   - *Profile-ML ψ, defined:* with the coefficients held at their converged
     value `α̂` (so `μ̂` is fixed), the NB log-likelihood becomes a function of
     ψ alone, `ℓ_g(ψ | μ̂)`; it is maximised by a bounded one-dimensional
     search over `log ψ` (`stats::optimize` on `[1e-3, 1e3]`). "Profile"
     because the other parameters are held at their optimum rather than
     integrated or moderated; the coefficients are then re-polished at the
     new ψ and the pair of steps repeated once. Contrast with `fitNB`'s
     edgeR-moderated dispersion (shrunk toward a cross-gene trend): the
     profile value is the gene's own, which is why a gene carrying no
     information about overdispersion returns a bound and keeps the
     moderated value (`psi_bound`), and why the standard error uses the
     Pearson working dispersion rather than ψ.
   - *The Schur-complement absorption, defined:* the design is `[X | Z]` with
     `Z` the nested indicators (one column per non-empty (sample, cell type)
     group, each cell in exactly one group). The Newton information is the
     block matrix `[[A, B], [Bᵀ, C]]` with `A = Xᵀ diag(w) X + Λ_x`,
     `B = Xᵀ diag(w) Z`, `C = Zᵀ diag(w) Z + Λ_z`. Because each cell belongs
     to one group, `C` is **diagonal** (the group's summed weight plus its
     penalty), so the step for the dense columns solves the `p × p` system
     with the Schur complement `S = A − B C⁻¹ Bᵀ`, and the group steps follow
     by back-substitution. The cost per Newton step is therefore one dense
     gram plus group sums, independent of the number of groups (660 on the
     cohort). Then the connection: `S` is exactly the information of the
     dense columns *after* partialling out the group indicators —
     Frisch–Waugh again, computed by algebra instead of by centring. The
     same absorption is used by inference for the per-gene covariance.
7. **Inference.** The penalised covariance, the Pearson dispersion, the
   Satterthwaite df (kept, with the toy df demonstration); then the missing
   pieces, written for the first time: the within-gene combination of a
   gene's correlated niche p-values (Cauchy two-sided by default, Brown
   one-sided; why sidedness differs), the cross-bandwidth Cauchy combination
   weighted by relative log-likelihood, and the three-level hierarchical BH
   cascade that produces `results()`. Niche-only mode as a subsection: what
   the design drops and why, what is tested.
8. **Reading a fit.** The slot table (`tau2`, `df`, `penalty`, `polish`), the
   diagnostics chunk from the quickstart, and the usage caveats (sample-level
   covariates, few samples, `re.maxit = 10` for slopes).
9. **The patient-level question.** `compositionTest()`: the pseudobulk, the
   mean niche density, the limma model, the two terms; one paragraph on
   reading it beside the niche results (the cohort example in one sentence).
10. **Choosing the estimator.** Short: `fitSpiDE(random = "intercept")` for
    discovery; `twoStageSpiDE()` for pre-specified hypotheses; pointer to C.

Removed: all benchmark numbers except a single sentence per default with a
pointer to D; "What the study found"; the design-decomposition narrative
(its *conclusion* is now the FWL paragraph (ii), its *numbers* go to D).

### C. `spiDE-twostage.Rmd` — The two-stage estimator (target ~330 lines, from ~470)

The current section as its own vignette: why a second estimator, the
pipeline figure, stage 1 (offset, joint fit, the other responses), pooling,
stage 2 (`tau2`, the moderated contrast), diagnostics, worked example, when
to use it. **Cut:** "how it measured", the withdrawn permutation pair, the
superseded-sweep story, the two "measurement notes". One paragraph of
operating characteristics with a pointer to D and the research site.

### D. `spiDE-calibration.Rmd` — Calibration and evidence (target ~330 lines, from 292)

Keep the λ material as is (definition, the below-1 demonstration, the cost
table, per-cell-type diagnosis, do-not-rescale, restrict what you test).
**Remove** the design arguments that duplicate B: the composition-driver
paragraph collapses to two sentences with a pointer; "do not reach for a
smaller model either" is deleted (its content is FWL (ii) in B); the
"third possibility" under λ > 1 collapses to a pointer.

**Add** one section, "What the benchmarks measured", which becomes the *only*
place in the package that quotes numbers, each as a short paragraph with the
research-site pointer:

- `random = "none"` vs `"intercept"` null type-I (0.71 vs 0.04).
- `df.method`: between vs satterthwaite (0.001 at S = 4; 0.042–0.065).
- the 0.99.17 switches on the synthetic null and what the simulator does not
  contain; the S-dependence of the two mechanisms; the two candidate cures.
- the niche-only design's recall as a different estimand (0.36 vs 0.17;
  0.7× SE; 41% more triplets on the cohort with the term).
- the two-stage arms (ols/nb/spanorm) in one paragraph.
- the real-cohort shuffle null: the confound, the fix, the calibrated null.

Retitled: "Calibration: reading lambda, and what the benchmarks measured".

## 3. Cross-cutting rules

- A number appears in **one** vignette (D). B and A say "measured; see the
  calibration vignette".
- History stays out of all four. What was refuted lives on the research site
  and in CLAUDE.md; a vignette says what the method *is*.
- Every symbol defined once, at first use, in B; A, C and D reuse B's
  notation and link to it.
- The FWL theorem is stated once (B §3) and referred to by name afterwards;
  "absorbed" and "partialled out" are used as synonyms for it explicitly.

## 4. Execution and verification

1. Write B first (the reference), then C (moved text), then A and D.
2. Knit all four standalone; then `R CMD check` (vignette rebuild) and
   `BiocCheck`; the built HTML total must stay well under the 10 MB tarball
   cap (today's three total 3.2 MB).
3. Update the cross-references in `man/` (`@seealso`), CLAUDE.md's "Where the
   evidence lives", `NEWS.md`, and `_pkgdown.yml` if it lists vignettes.
4. Commit per vignette, then push.

Estimated size: about 1,600 lines across four files (from 2,123 across
three), with the model vignette halved and the evidence in one place.

## 5. Decisions for the author

1. Four vignettes (split two-stage out) versus three (two-stage stays inside
   the model vignette as an appendix-style final section). Recommendation:
   four.
2. Keep the two pedagogical figures (glmmPQL trace, shrinkage arrows) and the
   pseudo-replication QQ, drop the Gaussian-density figure. Recommendation:
   as stated.
3. Whether "What the benchmarks measured" should carry a compact table per
   item rather than paragraphs. Recommendation: paragraphs with one small
   table for the switch arms.
