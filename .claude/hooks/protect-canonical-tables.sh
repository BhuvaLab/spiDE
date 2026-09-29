#!/usr/bin/env bash
# PreToolUse (Edit|Write): refuse direct writes to the canonical benchmark
# tables and to the shipped example dataset.
#
# WHY: a canonical table is THE table a report or the calibration vignette
# reads -- one table per scenario, where a new method / engine / depth arm is
# extra ROWS, never a parallel file that would carry a stale copy of the other
# arms. Each table is produced by its installer from aggregated outputs, never
# hand-edited; a direct write silently decouples a published number from the
# run that justified it. Two sets are protected:
#   research/bench2/tables/*.rds             the benchmark of spiDE >= 0.99.30,
#                                            written only by bench2's installer
#   research/reports/benchmarks/tables/*.rds the mixed model's (spiDE <= 0.99.22),
#                                            archived; research/R/install_results.R
#                                            and research/plasmode/install_twostage.R
#
# Matching is a bash `case` glob on the normalised path, not a grep regex, so it
# behaves the same under GNU grep and ugrep (Bunya's grep, where `^` inside a
# mid-pattern group never matches). In a glob `*` also crosses `/`, so a
# subdirectory of a tables/ directory is protected too.
#
# Exit 2 = block the call and return this message to Claude.
set -uo pipefail
input=$(cat)
f=$(printf '%s' "$input" | jq -r '.tool_input.file_path // .tool_input.filePath // empty' 2>/dev/null)
[ -z "$f" ] && exit 0
# resolve ./ and ../ so research/bench2/R/../tables/x.rds cannot slip past
if command -v realpath >/dev/null 2>&1; then
  n=$(realpath -m -- "$f" 2>/dev/null) && [ -n "$n" ] && f="$n"
fi
case "$f" in
  *research/bench2/tables/*.rds|*research/reports/benchmarks/tables/*.rds)
    cat >&2 <<MSG
BLOCKED: that is a canonical benchmark table.

  target: $f

One canonical table per scenario: a new method, engine or depth arm is extra
ROWS (its own value in the table's arm columns), never a parallel file that
would carry a stale copy of the other arms, and never a hand edit.

These tables are regenerated from aggregated outputs by their installer:
  research/bench2/tables/             bench2's installer (research/bench2/R/)
  research/reports/benchmarks/tables/ research/R/install_results.R,
                                      research/plasmode/install_twostage.R
                                      (the archived mixed model's; historical)

To add an arm, emit rows from the runner and re-run the installer. If you truly
intend to replace a canonical table, do it in a shell command the user can see
and approve, not through a file write.
MSG
    exit 2 ;;
  *data/toySpiDE.rda)
    cat >&2 <<'MSG'
BLOCKED: data/toySpiDE.rda is the shipped example dataset and is generated,
not written directly.

Regenerate it with `source("data-raw/make_toySpiDE.R")` (it runs
devtools::load_all() and calls the internal .toySPE()). The exported examples,
the vignettes and the tests depend on its planted G1/A/B effect and on its
sixteen patients (eight per condition, the condition test needing at least
three with usable slopes per condition); a hand edit or an ad-hoc save()
breaks that silently.
MSG
    exit 2 ;;
esac
exit 0
