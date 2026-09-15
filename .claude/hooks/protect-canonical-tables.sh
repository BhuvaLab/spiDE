#!/usr/bin/env bash
# PreToolUse: refuse direct writes to the canonical benchmark tables.
#
# WHY (CLAUDE.md, "Where the evidence lives"): research/reports/benchmarks/
# tables/*.rds are THE canonical tables the reports read -- "One canonical
# table per scenario: a new method arm is extra *rows*, not a parallel file
# that would carry a stale copy of the others." They are produced by the
# installers (research/R/install_results.R, research/plasmode/install_twostage.R)
# from aggregated part files, never hand-edited. A direct write silently
# decouples a published report from the run that justified it.
#
# Exit 2 = block the call and return this message to Claude.
set -uo pipefail
input=$(cat)
f=$(printf '%s' "$input" | jq -r '.tool_input.file_path // .tool_input.filePath // empty' 2>/dev/null)
[ -z "$f" ] && exit 0
case "$f" in
  *research/reports/benchmarks/tables/*.rds)
    cat >&2 <<'MSG'
BLOCKED: that is a canonical benchmark table.

CLAUDE.md: "One canonical table per scenario: a new method arm is extra *rows*,
not a parallel file that would carry a stale copy of the others."

These tables are regenerated from aggregated part files, not written directly:
  research/R/install_results.R
  research/plasmode/install_twostage.R

To add a method arm, emit rows from the benchmark runner (method / df.method /
stage1 / ls.model columns) and re-run the installer. If you truly intend to
replace the canonical table, do it in a shell command the user can see and
approve, not through a file write.
MSG
    exit 2 ;;
  *data/toySpiDE.rda)
    cat >&2 <<'MSG'
BLOCKED: data/toySpiDE.rda is the shipped example dataset and is generated,
not written directly.

CLAUDE.md: regenerate it with `source("data-raw/make_toySpiDE.R")` after
`devtools::load_all()` (it calls the internal .toySPE()). The exported
examples and the vignette depend on its planted G1/A/B effect and on the
field/n_per tuning that lets all four default bandwidths converge; a hand
edit or an ad-hoc save() breaks that silently.
MSG
    exit 2 ;;
esac
exit 0
