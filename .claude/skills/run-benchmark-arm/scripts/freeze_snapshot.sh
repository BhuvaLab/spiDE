#!/usr/bin/env bash
# Freeze an immutable snapshot of the spiDE package for a long HPC run.
#
#   freeze_snapshot.sh [label]   -> prints the snapshot path on stdout
#
# WHY: SLURM jobs load the package with devtools::load_all() on the LIVE working
# tree at task start. Array tasks start at different times, so an edit made
# while a run is in flight silently gives later tasks DIFFERENT code. This has
# corrupted 94 in-flight tasks on this project once, and nearly did so again on
# 2026-08-13. Point SPIDE_PKG at the snapshot instead.
set -euo pipefail
SRC="${SPIDE_SRC:-/scratch/project_mnt/S0249/R_projects/spiDE}"
LABEL="${1:-run}"
SNAP="/scratch/project_mnt/S0249/R_projects/spiDE_snapshots/$(date +%Y%m%d_%H%M%S)_${LABEL}"
mkdir -p "$SNAP"
rsync -a --exclude='.git' --exclude='research' --exclude='docs' \
      --exclude='*.tar.gz' --exclude='.Rproj.user' "$SRC"/ "$SNAP"/
{ echo "spiDE snapshot for an HPC run"
  echo "taken:  $(date -Is)"
  echo "label:  $LABEL"
  echo "source: $SRC"
  ( cd "$SRC" && echo "branch: $(git branch --show-current)" \
      && echo "HEAD:   $(git rev-parse HEAD)" \
      && echo && echo "--- git status ---" && git status --porcelain \
      && echo && echo "--- diff stat vs HEAD ---" && git diff --stat )
} > "$SNAP/SNAPSHOT_PROVENANCE.txt"
# The Rscript drivers too: R parses a script incrementally, so a task that is
# inside its multi-hour fit reads the NEXT expression from whatever the live
# file holds by then. The sbatch scripts run the copy under drivers/ when it
# exists (and copy it again to a task-private file at task start).
mkdir -p "$SNAP/drivers"
for d in research/R/run_task.R research/fdr-ordering/R/package_fixed_design.R; do
  [ -f "$SRC/$d" ] && cp "$SRC/$d" "$SNAP/drivers/"
done
chmod -R a-w "$SNAP/R" "$SNAP/drivers" 2>/dev/null || true
echo "$SNAP"
