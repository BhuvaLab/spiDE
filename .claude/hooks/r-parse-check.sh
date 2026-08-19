#!/usr/bin/env bash
# PostToolUse: parse-check any .R file that was just edited or written.
#
# WHY: a scripted edit to R/twostage-stage1.R once replaced far more than
# intended (151 insertions for a small substitution). It was caught only
# because parse() happened to be run by hand. R is not compiled, so a
# structurally broken file stays silent until something loads it -- which on
# this project can be an 80-minute SLURM job.
#
# Advisory, never blocking: exit 0 always. A parse failure is reported on
# stderr so it reaches the transcript immediately.
set -uo pipefail
input=$(cat)
f=$(printf '%s' "$input" | jq -r '.tool_input.file_path // .tool_input.filePath // empty' 2>/dev/null)
[ -z "$f" ] && exit 0
case "$f" in
  *.R|*.r) ;;
  *) exit 0 ;;
esac
[ -f "$f" ] || exit 0
if ! err=$(Rscript --vanilla -e "invisible(parse('${f//\'/\\\'}'))" 2>&1); then
  echo "R PARSE FAILURE in $f" >&2
  echo "$err" | head -5 >&2
  echo "The file is syntactically broken -- fix before anything loads it." >&2
fi
exit 0
