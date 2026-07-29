# Reply to the statistician's review of the spiDE mixed-effects model

Thank you — these are exactly the right questions, and running them changed how we frame the demonstration. Short version: **your mechanism is real and we can reproduce it, but it is only part of the story; the anti-conservatism does *not* disappear once the group means are equalised, because the real problem is the *error stratum*, not the group-mean difference. A properly-tested nested fixed-effects model fixes it exactly as you say — and that fix is, for the main effect, equivalent to the random-effects model. We prefer the RE model for scalability and, decisively, for the niche-*slope* tests that are the actual target.**

All numbers below are from `analysis/statistician-review/re-vs-fixed.R` (the vignette's own `sim_clustered` design; `n_samples = 12`, `n_per = 90`, `n_genes = 200`, `sd_patient = 0.7`, 5 seeds; Response-effect null type-I at α = 0.05, so the honest target is ≈ 0.05).

## 1. How the individual effects are generated

`sim_clustered()` draws a per-(gene, sample) intercept `u[g,s] ~ N(0, 0.7^2)`, **i.i.d. and independent of condition**, adds it on the log-mean scale (`log μ = base[celltype] + u[g,s]`), and generates NB counts (`size = 5`). Condition is a fixed patient-level label (alternating across samples); cells per sample are equal. So for each gene the responder and non-responder group means of `u` differ only by chance, \~`N(0, 4·τ²/S)`. This is precisely the setup you diagnosed.

## 2. Your zero-mean-within-group test — you are right about the mechanism

We re-ran with the per-sample effects centred to **zero mean within each group** (your `c(-0.5, 0.5, -0.3, 0.3)` generalised), and separately with non-zero but **equal** group means:

| per-sample effects `u`     | fixed-effects type-I | mixed type-I | mean\|Δ\| |
|----------------------------|----------------------|--------------|-----------|
| i.i.d. (current)           | **0.835**            | 0.020        | 0.32      |
| zero-mean within group     | 0.614                | 0.000        | 0         |
| equal non-zero group means | 0.644                | 0.000        | 0         |
| responder shifted by Δ=0.6 | **0.993**            | 0.015        | 0.60      |

and confirmed the mechanism directly: under the i.i.d. design the per-gene fixed z-statistic correlates with the realised per-gene group-mean difference `Δ_g` (cor ≈ 0.54). So yes — a chance difference in the group means of `u` is picked up by `resp`, and inflating it (row 4) drives type-I to \~1. That part is exactly as you describe.

## 3. But equalising the group means does *not* remove the bias

This is the part that surprised us. Even with the group means of `u` **exactly equal** (rows 2–3, mean\|Δ\| = 0), the fixed-effects test is still badly anti-conservative (\~0.61–0.64, versus the nominal 0.05). Two reasons, and the second is the important one:

-   **The log link.** Centring the effects on the log-mean scale does not equalise the group means on the *count* scale: `E[exp(u)]` depends on the within-group spread of the realised `u`, which differs between two finite groups. So a count-scale group difference survives even when the log-scale means match — and it *grows*, relative to its cell-level SE, as cells per sample increase (the signature of pseudo-replication, not of a point-estimate bias).
-   **The error stratum (the crux).** The fixed cell-level test judges `resp` against `~C` cells of residual, when the honest reference is the `S` patients. This is independent of the group mean. We show it cleanly in the next section.

So the i.i.d. design is not an unfair artifact: real patient effects are *not* constrained to zero-mean-within-group, and even if they were, a count-scale between-patient signal remains that must be tested at the between-patient stratum.

## 4. Your nested (contr.sum) fixed model — the key experiment

You are right that the naive "sample as an independent fixed effect" is collinear with `resp`, and that nesting individuals within group via `contr.sum` makes `resp` estimable. We built exactly that (`A=(-1,0), B=(1,0), C=(0,-1), D=(0,1)` within each group) and tested `resp` three ways, on the same i.i.d. data (Gaussian on `log1p` counts, so the arithmetic is transparent):

| test of `resp` on the nested model                    | type-I   |
|-------------------------------------------------------|----------|
| Wald, **cell-level** residual SE                      | **0.84** |
| F against the **between-individual** stratum (S−2 df) | 0.013    |
| random intercept (`lmer`, Satterthwaite)              | 0.033    |

This is the whole argument in one table. **The nesting alone does not fix anything** — tested at the cell residual it is as anti-conservative as before (0.84). What fixes it is testing `resp` against the **between-individual error stratum** (0.013). And the random-effects model gives the *same* calibrated answer (0.033). So:

-   You are correct that a *properly tested* nested fixed-effects model is calibrated and avoids the bias.
-   The random effect is not doing anything mysterious — its mean-zero prior is, as you say, the shrinkage counterpart of the sum-to-zero-within-group constraint; the operative change in both is moving `resp` to the between-patient stratum.

## 5. So do we need the random-effects model?

For the **main effect** η, no — a nested fixed model tested against the individual stratum is equivalent (row 2 ≈ row 3 above). And the niche interactions **β do not need random slopes either**: on our real data the between-patient slope variance τ²_slope ≈ 0 (the niche→expression slope does not vary appreciably between patients), so a random slope is shrunk away and adds nothing. (An earlier draft of this reply claimed RE was "decisively" needed for the slopes — that was wrong, and your question is what made us check.)

What the niche interactions actually need is the correct **reference df**, not a random slope. With τ²_slope ≈ 0 the niche covariate varies *within* patients, so β is a within-patient contrast carrying far more than S independent observations (the high end of an effective-df range spanning ≈ 4 to ≈ 2400 across the tested columns). Referencing it against S−2 is therefore badly **over-conservative** — the actual real-data symptom — and the fix is a per-coefficient (Satterthwaite) df that gives β its large within-patient df (≈ 775 in the toy fit, versus ≈ S−2 = 8 for the between-patient main effect). This is a degrees-of-freedom question, orthogonal to the random-vs-fixed debate.

So the honest case for the random-effects model, given τ²_slope ≈ 0, is:

1.  **Scalability and shrinkage.** One τ² is pooled across all ~13k genes (a closed-form Schall/MoM update per iteration), each patient shrunk toward the shared distribution — more stable for sparse genes than free per-patient parameters.
2.  **Convenience for NB / covariates / imbalance.** The mixed model handles the NB-GLM and nuisance covariates without hand-rolling error strata. In fairness, we *tested* the imbalance concern across a **1× to 10×** cell-count sweep and it did **not** materialise: at every fold both your nested-between test and the random-intercept model stay calibrated (≈ 0.04, essentially unchanged from balanced), while only the plain cell-level fixed test degrades (0.73 → 0.91) as imbalance grows. With a substantial τ² the unequal within-sample precisions are a minor perturbation, so imbalance is *not* a differentiator even at 10-fold — your nested-fixed approach is robust to it too. The mixed model's edge here is just that it does the right thing automatically rather than requiring the correct stratum to be assembled by hand.
3.  **The per-coefficient df falls out of the same framework.** The Satterthwaite df that relieves the β over-conservatism is a by-product of the mixed model, not a separate bolt-on.

Random *slopes* are deliberately **not** on this list: with τ²_slope ≈ 0 they are inert, so we do not recommend them for the default analysis. (A stress test that *plants* between-patient slope variation does leave even the random-slope model partly anti-conservative — a covariance-level effect — but that regime is exactly the one τ²_slope ≈ 0 says real data is not in, which is why we are not pursuing a Kenward–Roger covariance correction at this stage.)

## 6. Method-of-moments vs MLE for `τ²`

You are right that Schall/MoM is less efficient than (RE)ML. Two mitigating points, one empirical:

-   **Pooling across genes largely compensates.** The pooled MoM estimate on this design lands at `τ̂² = 0.454` (true 0.49). Per-gene REML (matched NB `glmmPQL`, θ = 5) on 20 genes gives mean 0.395, median 0.351, **sd 0.21** — i.e. the per-gene MLE is *noisier and slightly biased low* (PQL bias) than the pooled MoM, at \~100× the per-gene cost. At 13k genes the per-gene-optimisation MLE is the expensive path; the shared MoM update is closed-form.
-   We would happily revisit an MLE/REML variance component (your derivation) if a setting arises where `τ²`'s absolute accuracy, rather than the response-niche t-statistics it feeds, is what matters — but on the evidence here the pooled MoM is both cheaper and, per gene, no less accurate.

## Summary

-   Your mechanism is real: a chance group-mean difference in the per-sample effects is picked up by `resp` (§2, §5-decomposition).
-   But equalising the group means does **not** remove the anti-conservatism (§3): the real issue is the error stratum, not the group mean.
-   Your nested `contr.sum` model, tested against the between-individual stratum, **is** calibrated and equivalent to the random-intercept model for the main effect (§4) — the nesting alone (cell SE) is not enough.
-   We keep the random-**intercept** model for scalability/shrinkage/unbalanced-NB handling, and rely on the per-coefficient (Satterthwaite) df — **not** random slopes (τ²_slope ≈ 0, inert) — to relieve the over-conservative niche-interaction tests (§5).
-   MoM is less efficient than MLE, but pooling across genes makes it both cheaper and, per gene, no less accurate here (§6).

Happy to share the script (`analysis/statistician-review/re-vs-fixed.R`) and dig further into any of these.