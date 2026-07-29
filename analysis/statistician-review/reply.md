# Reply to the statistician's review of the spiDE mixed-effects model

Thank you — these are exactly the right questions, and running them changed how we
frame the demonstration. Short version: **your mechanism is real and we can
reproduce it, but it is only part of the story; the anti-conservatism does *not*
disappear once the group means are equalised, because the real problem is the
*error stratum*, not the group-mean difference. A properly-tested nested
fixed-effects model fixes it exactly as you say — and that fix is, for the main
effect, equivalent to the random-effects model. We prefer the RE model for
scalability and, decisively, for the niche-*slope* tests that are the actual
target.**

All numbers below are from `analysis/statistician-review/re-vs-fixed.R` (the
vignette's own `sim_clustered` design; `n_samples = 12`, `n_per = 90`,
`n_genes = 200`, `sd_patient = 0.7`, 5 seeds; Response-effect null type-I at
α = 0.05, so the honest target is ≈ 0.05).

## 1. How the individual effects are generated

`sim_clustered()` draws a per-(gene, sample) intercept `u[g,s] ~ N(0, 0.7^2)`,
**i.i.d. and independent of condition**, adds it on the log-mean scale
(`log μ = base[celltype] + u[g,s]`), and generates NB counts (`size = 5`).
Condition is a fixed patient-level label (alternating across samples); cells per
sample are equal. So for each gene the responder and non-responder group means of
`u` differ only by chance, ~`N(0, 4·τ²/S)`. This is precisely the setup you
diagnosed.

## 2. Your zero-mean-within-group test — you are right about the mechanism

We re-ran with the per-sample effects centred to **zero mean within each group**
(your `c(-0.5, 0.5, -0.3, 0.3)` generalised), and separately with non-zero but
**equal** group means:

| per-sample effects `u` | fixed-effects type-I | mixed type-I | mean\|Δ\| |
|---|---|---|---|
| i.i.d. (current) | **0.835** | 0.020 | 0.32 |
| zero-mean within group | 0.614 | 0.000 | 0 |
| equal non-zero group means | 0.644 | 0.000 | 0 |
| responder shifted by Δ=0.6 | **0.993** | 0.015 | 0.60 |

and confirmed the mechanism directly: under the i.i.d. design the per-gene fixed
z-statistic correlates with the realised per-gene group-mean difference `Δ_g`
(cor ≈ 0.54). So yes — a chance difference in the group means of `u` is picked up
by `resp`, and inflating it (row 4) drives type-I to ~1. That part is exactly as
you describe.

## 3. But equalising the group means does *not* remove the bias

This is the part that surprised us. Even with the group means of `u` **exactly
equal** (rows 2–3, mean\|Δ\| = 0), the fixed-effects test is still badly
anti-conservative (~0.61–0.64, versus the nominal 0.05). Two reasons, and the
second is the important one:

- **The log link.** Centring the effects on the log-mean scale does not equalise
  the group means on the *count* scale: `E[exp(u)]` depends on the within-group
  spread of the realised `u`, which differs between two finite groups. So a
  count-scale group difference survives even when the log-scale means match — and
  it *grows*, relative to its cell-level SE, as cells per sample increase (the
  signature of pseudo-replication, not of a point-estimate bias).
- **The error stratum (the crux).** The fixed cell-level test judges `resp`
  against `~C` cells of residual, when the honest reference is the `S` patients.
  This is independent of the group mean. We show it cleanly in the next section.

So the i.i.d. design is not an unfair artifact: real patient effects are *not*
constrained to zero-mean-within-group, and even if they were, a count-scale
between-patient signal remains that must be tested at the between-patient
stratum.

## 4. Your nested (contr.sum) fixed model — the key experiment

You are right that the naive "sample as an independent fixed effect" is collinear
with `resp`, and that nesting individuals within group via `contr.sum` makes
`resp` estimable. We built exactly that (`A=(-1,0), B=(1,0), C=(0,-1), D=(0,1)`
within each group) and tested `resp` three ways, on the same i.i.d. data
(Gaussian on `log1p` counts, so the arithmetic is transparent):

| test of `resp` on the nested model | type-I |
|---|---|
| Wald, **cell-level** residual SE | **0.84** |
| F against the **between-individual** stratum (S−2 df) | 0.013 |
| random intercept (`lmer`, Satterthwaite) | 0.033 |

This is the whole argument in one table. **The nesting alone does not fix
anything** — tested at the cell residual it is as anti-conservative as before
(0.84). What fixes it is testing `resp` against the **between-individual error
stratum** (0.013). And the random-effects model gives the *same* calibrated
answer (0.033). So:

- You are correct that a *properly tested* nested fixed-effects model is
  calibrated and avoids the bias.
- The random effect is not doing anything mysterious — its mean-zero prior is, as
  you say, the shrinkage counterpart of the sum-to-zero-within-group constraint;
  the operative change in both is moving `resp` to the between-patient stratum.

## 5. So do we need the random-effects model?

For the **main effect**, no — a nested fixed model with an F-test against the
individual stratum is equivalent (row 2 ≈ row 3 above). We use random effects for
three practical reasons, the last decisive:

1. **Scalability and shrinkage.** One `τ²` is pooled across all ~13k genes
   (a single Schall/MoM update per iteration), and each patient's estimate is
   shrunk toward the shared distribution — more stable for the many sparse genes
   than a free per-patient parameter.
2. **Unbalanced / covariate-adjusted / NB designs.** The clean ANOVA error-strata
   arithmetic degrades with unequal cells, nuisance covariates, and the NB-GLM;
   the mixed-model machinery handles all of these uniformly.
3. **The niche-slope tests `β` — the actual target.** spiDE's scientific quantity
   is the response×celltype×niche interaction, i.e. how the *slope* of expression
   on local niche density differs by condition. Under a null with genuine
   between-sample slope variation, the fixed test is again grossly
   anti-conservative, and the fixed-nested analogue would require
   individual×niche interaction contrasts — a large, ill-conditioned block that is
   collinear with the effect being tested. The random-*slope* model is the only
   practical route:

   > Under a null with genuine between-sample niche-slope variation (true
   > β = 0), the fixed test is grossly anti-conservative (type-I ≈ 0.45–0.59).
   > The random-*slope* model roughly halves it (≈ 0.24–0.29) but does **not**
   > fully calibrate it, and — unlike the main effect — this does not improve
   > materially from S = 12 to S = 24. The per-coefficient (Satterthwaite) df
   > makes it slightly *worse* (≈ 0.31–0.36), which localises the residual to the
   > **covariance** (the reduced-design SE), not the reference df.

   We are being candid about this: full calibration of the slope tests is the
   genuinely hard part, and on this small-S / many-cells toy design the
   random-slope model reduces but does not eliminate the anti-conservatism — a
   known small-sample property of the reduced-design Wald covariance, and the
   reason a Kenward–Roger covariance correction is held in reserve. What the
   experiment *does* establish is the ordering: the fixed cell-level test is
   unusable for the slopes, the fixed-nested analogue (individual×niche contrasts)
   is intractable, and the random-slope model is the only practical route and the
   least anti-conservative of the three. Calibration of the slope tests on
   realistic-scale data is being assessed separately.

## 6. Method-of-moments vs MLE for `τ²`

You are right that Schall/MoM is less efficient than (RE)ML. Two mitigating
points, one empirical:

- **Pooling across genes largely compensates.** The pooled MoM estimate on this
  design lands at `τ̂² = 0.454` (true 0.49). Per-gene REML (matched NB `glmmPQL`,
  θ = 5) on 20 genes gives mean 0.395, median 0.351, **sd 0.21** — i.e. the
  per-gene MLE is *noisier and slightly biased low* (PQL bias) than the pooled
  MoM, at ~100× the per-gene cost. At 13k genes the per-gene-optimisation MLE is
  the expensive path; the shared MoM update is closed-form.
- We would happily revisit an MLE/REML variance component (your derivation) if a
  setting arises where `τ²`'s absolute accuracy, rather than the response-niche
  t-statistics it feeds, is what matters — but on the evidence here the pooled MoM
  is both cheaper and, per gene, no less accurate.

## Summary

- Your mechanism is real: a chance group-mean difference in the per-sample effects
  is picked up by `resp` (§2, §5-decomposition).
- But equalising the group means does **not** remove the anti-conservatism (§3):
  the real issue is the error stratum, not the group mean.
- Your nested `contr.sum` model, tested against the between-individual stratum,
  **is** calibrated and equivalent to the random-intercept model for the main
  effect (§4) — the nesting alone (cell SE) is not enough.
- We keep the random-effects model for scalability/shrinkage/NB and, decisively,
  for the niche-*slope* tests where the fixed-nested analogue is intractable (§5).
- MoM is less efficient than MLE, but pooling across genes makes it both cheaper
  and, per gene, no less accurate here (§6).

Happy to share the script (`analysis/statistician-review/re-vs-fixed.R`) and dig
further into any of these.
