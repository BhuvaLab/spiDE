---
name: run-benchmark-arm
description: Launch a spiDE benchmark arm, null-grid set or real-cohort run on Bunya with reproducibility built in - pre-registered criteria, frozen code (package snapshot or pinned library), task-private driver copies, recorded seed and provenance, atomic writes, resume-on-partial, memory sized from measurement, and every arm scored on the public cohorts. Use when launching a new engine/depth arm, a block or permutation null, a bench2 scenario, or a real-cohort run whose numbers anyone will later cite.
disable-model-invocation: true
---

# Run a benchmark arm

Every rule here exists because its absence cost a run on this project. Follow
them in order; each step is cheap next to the multi-hour job it protects.

The harnesses:
- **`research/bench2/`** is the new benchmark: a template-fitted simulator
  scoring depth handling and the two engines. It is being built.
- **`research/simplify/`** is the study behind the 0.99.30 engines, with its
  nulls and plasmode.
- **`research/R/run_task.R`** is the old harness, with its tables in
  `research/reports/benchmarks/tables/`. It runs the **mixed model** of spiDE
  <= 0.99.22 from pinned legacy libraries (`research/libs/mixed-final` for
  `spiDEmixed`, and the design libraries of `research/R/build_design_libs.R`).
  It is historical: use it only to reproduce an archived number, and call the
  old functions package-qualified (`spiDEmixed::fitSpiDE()`).

## 0. Pre-register

Write the question, the arms, the scoring and the pass/fail criteria into
the study's `README.md` **before** the results they govern. Examples are
`research/simplify/README.md` and `research/bench2/README.md` (the realism
gate). The README is then edited only to add dated decisions, never to move
a criterion after its data have been seen.

Name everything in words: engines, depth arms, scenarios and gates (the
mechanism, power-ceiling, calibration and power gates). Never use a
letter-and-number code.

## 1. Freeze the code

There are two routes, depending on how the driver loads spiDE.

**A driver that `load_all()`s a tree** (the old harness, cohort drivers):

```bash
SNAP=$(.claude/skills/run-benchmark-arm/scripts/freeze_snapshot.sh <label>)
```

Pass `SPIDE_PKG=$SNAP` to every task. SLURM tasks load the tree at *task
start*, so array tasks starting at different times pick up different code if
you edit while a run is in flight. That has corrupted 94 in-flight tasks on
this project. The snapshot records the branch, HEAD and the uncommitted diff
in `SNAPSHOT_PROVENANCE.txt`.

**A driver that loads an installed build** (simplify, bench2: their
`R/common.R` prepends `research/libs/<name>`): install spiDE and SpaNorm from
a frozen snapshot into a pinned library `research/libs/<name>`. Record the
source snapshot, the versions and the date in `research/libs/<name>.SNAPSHOT`,
as `simplify.SNAPSHOT` and `mixed-final.SNAPSHOT` do. Never reinstall a
pinned library while an array reads it.

The default user library holds an old spiDE (0.99.19, the mixed model) and
SpaNorm 1.7.10. A driver that forgets its `.libPaths()` line runs the wrong
model; the per-patient engines need SpaNorm >= 1.7.14.

## 2. Seed, and record the seed

The package's fit and test path draws no random numbers (0.99.30). The
harness around it does:
- the block null (`.blockShuffle()` in `research/simplify/R/common.R`,
  `set.seed(7000 + grid)`);
- the within-slide permutations (`.permLabels()`; the study drew seed
  20260928 for development and a fresh seed for each held-out or validation
  round: 20260930, 20261001, 20261002);
- plasmode plants;
- the bench2 simulator.

Derive each task's seed from its grid, permutation or replicate index, so
tasks stay independent. Write it into the output alongside:
- the spiDE and SpaNorm versions (`fit@params$version` records spiDE's);
- the library or snapshot path;
- the engine, `sigma`, `depth`, `covariates`, `strata`, `procedure` and
  `fdr`;
- the index types and the gene filter.

Keep validation seeds unseen until the candidates are registered: the study
drew a fresh block and permutation set for each validation round.

## 3. Size from measurement, not from fear

Use the measured per-grid timings and the procedure of the `hpc-job-sizing`
skill. Check `sstat -j <jobid> --format=MaxRSS` on the *same code path*, then
add about 3x. Over-requesting is not free insurance: the account QOS caps
mean it throttles your own concurrency. On 2026-08-13 a mixed-model cohort
run peaked at 21.7 GB against a 220 GB request, and three 24-core jobs
saturated `QOSGrpCpuLimit` at 80 CPUs.

For block nulls, remember that every grid is a refit. For permutation nulls
of the slopes engine only `testSpiDE()` reruns, because the fit never sees
the condition: one fit serves a thousand permutations. The sandwich engine
refits per permutation.

## 4. Write atomically, and resume

Write each output to `<path>.tmp$$`, then `mv` it into place, so a job killed
mid-write leaves the previous file intact rather than a truncated one. On
restart, load existing outputs and skip completed units. Carry the PRIOR
list forward rather than rebuilding it, or the resume silently drops earlier
replicates. Key outputs by a stable identity (cohort, grid, seed), never by
array position.

**Do not edit a driver an array is running.** R parses a driver script
incrementally through an 8 KB buffer, so a task inside its fit reads the
next expression from the modified file at a stale byte offset. On 2026-09-07
that killed 255 benchmark tasks *after* they had written their result (array
marked FAILED, `afterok` aggregation at `DependencyNeverSatisfied`; cancel it
and re-run the aggregation without the dependency, since the outputs are
intact). It then killed eight cohort tasks six to seven hours into their fit.

Every sbatch runs a **task-private copy** of its driver:
- the simplify and bench2 sbatch files copy `R/*.R` to `$(mktemp -d)` at task
  start;
- the old harness prefers the frozen copy `freeze_snapshot.sh` places under
  `$SPIDE_PKG/drivers/`.

Keep that pattern in any new sbatch. A task-private copy is taken at *task
start*, so an edit between the starts of two tasks still splits the array:
edit between arrays, never during one.

To rescue a live task, find its read position via
`srun --overlap --cpu-bind=none --jobid=<raw task id>`. A reader that has not
yet refilled its buffer (`pos: 8192` in `/proc/<pid>/fdinfo/3`) survives if
the file is truncated and rewritten **in place** with the bytes it expects
(`git show <rev>:path > path`). `git checkout` and `sed -i` create a new
inode and do not reach it. The `protect-running-drivers` hook blocks in-place
writes to `research/**` drivers while jobs run.

Inside the driver, fork no more BLAS threads than cores: export
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1` in the sbatch.
spiDE's `.bpGenes()` and SpaNorm's `polishNB()` set one BLAS thread per worker
themselves (when `RhpcBLASctl` is installed); a driver's own `mclapply` does
not. From inside an allocation, pass the environment explicitly
(`sbatch --export=ALL,VAR=...`), because `SBATCH_EXPORT=NONE` strips it, and
read the first task's log before calling a submission done.

## 5. Rows, not new files

A new arm is extra **rows** in the one canonical table per scenario, never a
parallel file that would carry a stale copy of the other arms:
- bench2: `research/bench2/tables/*.rds`, with the method, engine and depth
  arm as columns;
- the old harness: `research/reports/benchmarks/tables/*.rds`, the mixed
  model's, historical.

Only the table's installer writes them (for the old harness,
`research/R/install_results.R` and `research/plasmode/install_twostage.R`).
The `protect-canonical-tables` hook blocks direct writes to both sets.

## 6. Score on every cohort

Every arm is scored on the three public cohorts (GSE250346, GSE282639,
GSE289194) as well as on YTMA. A gate or a winner needs agreement across
them, or the disagreement is reported. YTMA is down-weighted (user directive
2026-09-29, weight 0.25 per YTMA cohort): a conclusion resting on YTMA alone
is unsupported until a public cohort replicates it. Where a pre-registered
rule weighted cohorts equally, report its verdict unchanged next to the
re-weighted one.

## 7. Gate anything GPU

spiDE 0.99.30 has no GPU path; SpaNorm's GPU backend does (e.g. template
fits).
- GPU work needs h100: a100 fails fp64 `digamma` inside NB dispersion
  *mid-run*, and l40s runs fp64 at about 1:64.
- Gate the real run behind a verification job with `--dependency=afterok`,
  and make sure the gate cannot pass vacuously. A gate that `source()`d a
  script ending in `quit()` "passed" in 68 seconds having tested nothing
  beyond its first layer. Set a flag per layer and assert all of them at the
  end.
- Check `squeue -h -j <id> -o %S` before promising a turnaround.

## 8. Score before quoting

Run the `calibration-check` skill on the output, including the arm's own
block and permutation nulls, before reporting any discovery count. Then
record the result with the `record-finding` skill.
