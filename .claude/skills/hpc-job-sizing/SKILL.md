---
name: hpc-job-sizing
description: Determine the best SLURM job structure (array vs. monolithic, CPU vs. GPU, partition/QOS, resource sizing) for a computational job on Bunya HPC. Use when asked to run something "on the cluster" / "via sbatch", to size a batch job's cores/mem/walltime, to decide whether a job should be an array, or whether GPU makes sense for a workload. Encodes lessons learned the hard way on this project's mixed-effects diagnostics workflow (research/diagnostics/step03/06/08).
---

# Sizing an HPC job on Bunya

This skill is a decision procedure, not a one-shot script: work through it in
order for any new job, and update it (via a short PR-style note in this file)
if you learn something that contradicts it. It exists because guessing
resource sizes on this cluster has repeatedly cost hours of wasted compute --
every rule below traces back to a real failure.

**Golden rule: never guess sizing at a new scale. Calibrate on one representative
unit of work first, look at real `sacct` numbers, then size the batch.** A
`sbatch --test-only` dry run only validates syntax/directives -- it tells you
nothing about whether the job will actually fit in the memory/time you asked
for.

## Step 1 -- characterize the job

Answer these before touching `sbatch`:

- **Is it embarrassingly parallel?** Many independent units of work (genes,
  samples, simulated datasets, held-out folds) that don't need to share state
  -> favors an **array** job. A single fit that can't be decomposed -> a
  monolithic job (optionally with internal `BiocParallel`/OpenMP parallelism).
- **Is there a fault-prone stage and an expensive-but-reliable stage in the
  same script?** E.g. an hours-long model fit followed by a
  memory-hungry inference/test stage that's more likely to OOM or crash.
  If so, **checkpoint after the reliable stage and before the fragile one**,
  so a crash in the fragile stage doesn't force redoing the expensive part.
  This is not optional polish -- it's what separates a 5-minute retry from
  losing 4+ hours of compute (see `research/diagnostics/step03_prep_and_fit_intercept.R`
  and `step08_diag5_slope_needed.R`, both of which learned this after an OOM
  ate a completed fit).
- **What's the actual data scale?** Row/column counts of the real matrices
  involved, not a stand-in dataset. Scaling from a smaller calibration point
  to a bigger one is *not* safe to assume linear -- treat it as a rough lower
  bound, not a real estimate (see Step 2).

## Step 2 -- calibrate before batching

1. `sbatch --test-only <script>.sbatch` first, always -- catches directive
   typos, invalid partition/QOS/gres combos, and bad `--array` syntax for
   free, in about 1 second, before spending a queue slot.
2. If the job is an array (or otherwise has many similar units), **submit
   ONE task first** -- ideally the most expensive/representative one (e.g.
   the single full-data fit-only run, or the largest held-out subset) --
   sized generously but reasoned, and confirmed with the user (see Step 4).
3. Read the real numbers from `sacct -j <jobid> --format=JobID,State,ExitCode,Elapsed,MaxRSS -n`
   once it finishes. `MaxRSS` on the `.0` (actual srun/task) step is the
   number that matters, not the `.batch` step (that's just the wrapper
   shell, a few MB). Only *then* size mem/time for the full batch, with
   comfortable margin above the observed peak (not the requested amount --
   the requested amount tells you nothing about the real need).
4. If a later stage of the *same* script needs more resources than the
   calibrated stage (e.g. a blocked multi-worker inference step after a
   single-threaded fit), that's a **separate, larger calibration problem** --
   don't assume the fit's footprint bounds the whole job (this bit us in
   step03 and step08: both stages OOM'd independently, at very different
   memory scales, within the same script).

## Step 3 -- diagnose queue delays before resizing blindly

If a job sits `PD` for a long time, check *why* before assuming your request
is the problem:

```bash
squeue -j <jobid> -o "%.10i %.9P %.8u %.2t %.10M %.6D %R"   # look at the REASON column
sacctmgr show qos <qos_name> format=Name,MaxTRESPerUser,GrpTRES,MaxTRESPA
squeue -A <account>                                          # is it just YOUR jobs?
```

`QOSGrpCpuLimit` / `Priority` often mean **cluster-wide** contention on a
shared QOS pool (`GrpTRES`), not something your own core/mem request can fix
-- reducing your request may not shorten the wait at all if the whole QOS
pool is saturated by other users. Distinguish "my account is capped" from
"the cluster is busy right now" before resizing reactively.

## Step 4 -- always confirm with the user before submitting

Per this project's standing HPC rules: never guess account/partition/QOS/
resource sizing, and never run a real fit on the login/interactive node
(only syntax checks belong there). Concretely:

- State the exact `--cpus-per-task`/`--mem`/`--time` you intend to use and
  *why* (calibration data point, or explicit reuse of a previously-proven
  config for a structurally identical job) before calling `sbatch` for real.
- If a job fails and you're resubmitting with different resources, say what
  changed and why (OOM -> more mem; timeout -> more time or more parallelism)
  -- don't silently retry with a bigger number.
- Cancelling a job that's still `PD` (never started) costs nothing -- check
  `squeue` before assuming a resize means lost work. Cancelling a `R`unning
  job does lose whatever that job hadn't checkpointed yet.

## Step 5 -- array vs. monolithic, concretely

Prefer an array over one job doing internal `BiocParallel::MulticoreParam`-style
forking when the units of work are each individually expensive (minutes to
hours) and independent:

| | Array (one task per unit) | Monolithic (internal fork) |
|---|---|---|
| Fault isolation | one bad unit costs one resubmit | one worker crash can kill the whole `bplapply` (seen firsthand: an OOM'd worker surfaced as BiocParallel's cryptic `wrong args for environment subassignment`, taking down every other in-flight result with it) |
| Memory accounting | clean per-task `MaxRSS` from `sacct` | conflated across all forked workers -- can't isolate one worker's real footprint |
| Scheduling | many small requests backfill into gaps easily | one job needing many cores/a large memory block at once schedules less easily on a contended QOS |
| Overhead | more submissions, needs an aggregation step afterward | simpler to write for genuinely small units of work |

Key naming lesson: **key per-task outputs by a stable identity (e.g. sample
name), not by array index**, if the set of units could ever change size or
selection. An earlier version of this project's LOO jackknife capped the
array at a random subset of samples and keyed outputs by position -- when
the scope later expanded to the full sample set, the index-to-sample mapping
silently shifted and would have corrupted already-computed results if not
caught first.

## Step 6 -- BLAS / CPU-bind gotchas specific to this cluster (Bunya)

Both bit us for real; apply both every time an R script does
`BiocParallel::MulticoreParam` or similar forking:

1. **BLAS oversubscription.** A multithreaded BLAS (OpenBLAS/MKL) competing
   with N forked workers oversubscribes the node by ~N^2 and can crash a
   worker mid-job -- surfaces as the same cryptic
   `wrong args for environment subassignment` error from BiocParallel's
   result-collection code, not an obviously BLAS-related message. Fix, both
   layers every time:
   - In the R script: `RhpcBLASctl::blas_set_num_threads(1)` and
     `RhpcBLASctl::omp_set_num_threads(1)`, before any fitting call.
   - In the sbatch script: `export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
     MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1` before `srun`.
2. **`srun` CPU-bind conflicts when submitting from within an active
   interactive job.** `srun --export=ALL` can re-export the *submitting*
   interactive job's `SLURM_CPU_BIND*` env vars into the new batch job,
   which then tries to apply a bind mask sized for a different node's core
   layout -- fails fast with `CPU binding outside of job step allocation`.
   Fix: `unset SLURM_CPU_BIND SLURM_CPU_BIND_LIST SLURM_CPU_BIND_TYPE
   SLURM_CPU_BIND_VERBOSE` and add `--cpu-bind=none` to the `srun` call in
   every sbatch script submitted from an interactive session.

## Step 7 -- CPU vs. GPU

Don't reach for GPU just because a workload is slow or because smaller GPUs
have shorter queues -- check these first, in order:

1. **Does the tool actually have a working GPU code path, right now, in this
   environment?** Check for the accelerated backend's actual dependency
   (e.g. `requireNamespace("torch", quietly = TRUE)`), not just a
   `backend = "gpu"` argument existing in a function signature -- an
   argument accepting `"gpu"` doesn't mean the dependency is installed.
2. **Is the bottleneck actually compute time, not memory?** If every failure
   so far has been an OOM (not a walltime timeout on an otherwise-healthy
   job), GPU doesn't address that -- GPU memory is a *harder*, smaller,
   per-device ceiling than CPU system RAM, not a bigger one.
3. **Does the workload's real memory footprint fit the specific GPU tier
   you'd request?** "Smaller/less busy" GPUs have less VRAM, not more --
   check the actual card's memory (see the table below) against your
   calibrated CPU footprint (Step 2) before assuming it'll fit. Workloads
   that stay sparse/un-densified on CPU often *increase* their memory
   footprint when ported to GPU tensors (most GPU tensor stacks favor dense
   layouts), which can make a memory-bound CPU problem worse on GPU, not
   better.
4. **Never adopt a "reinstall the GPU toolchain every job" pattern.** A
   CUDA-built `torch`/libtorch works across all NVIDIA cards in the same
   partition (they share a CUDA driver stack) without reinstalling per job;
   only a genuinely different accelerator family (e.g. AMD ROCm vs. NVIDIA
   CUDA) needs a separate persistent install, not a fresh one per
   submission. Compute nodes may not have outbound internet access at all,
   so a per-job install step is also a second, silent point of failure.
   Set up once per architecture/partition family (per the Bunya guide's own
   advice: log onto a GPU node interactively first to see what modules are
   available for that architecture), and confirm with the user before
   installing or upgrading anything.

## Bunya reference facts

(As of the last time this was checked -- reverify against
[the Bunya User Guide](https://github.com/UQ-RCC/hpc-docs/blob/main/guides/Bunya-User-Guide.md)
if anything here looks stale, especially GPU counts/types.)

- **CPU partitions**: `general` is the standard partition; QOS `normal` has
  a **cluster-wide** `GrpTRES` cpu cap shared by all users (not per-account),
  which is what an unexplained long `PD`/`QOSGrpCpuLimit` wait usually means.
- **GPU partitions**: `gpu_cuda` (NVIDIA: A100, H100, L40, L40s),
  `gpu_rocm` (AMD: Mi210, Mi300x), `gpu_sxm` (H100 SXM5, approved users
  only), `gpu_viz` (visualization: L40/L40s/A16, on-Bunya only).
- **GPU memory tiers** (check against your calibrated footprint before
  requesting): H100 80GB, A100 80GB, A100 MIG slices 10-40GB (fractional/
  shared), L40/L40s 48GB, Mi210 64GB, Mi300x 192GB.
- **Requesting a GPU**: `--gres=gpu:<type>:<count>` (e.g. `--gres=gpu:a100:1`,
  `--gres=gpu:l40:1`; MIG slices use the longer
  `--gres=gpu:nvidia_a100_80gb_pcie_1g.10gb:1` form) plus `--qos=gpu` (or
  `--qos=mig` for MIG slices). Max 4 H100s across all QOS.
- **CPU architecture split**: newer EPYC4 nodes and older EPYC3 nodes are
  *not* binary-compatible for compiled R/Python packages -- a package built
  on one won't necessarily run on the other. If a job mysteriously fails
  only on some nodes, check which CPU generation it landed on.
- **Software**: modules differ per node/architecture; the guide's own advice
  is to start an interactive session on the target node type first to see
  what's available there before assuming a module exists cluster-wide.
