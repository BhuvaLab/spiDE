#!/usr/bin/env bash
# Stop: warn once when a roxygen block in R/ changed but NAMESPACE and man/
# did not, i.e. devtools::document() has not been run since the edit.
#
# WHY: CI (.github/workflows/check-bioc.yml) runs R CMD check, which fails on
# a stale NAMESPACE or man/*.Rd, and BiocCheck flags undocumented exports.
# The parse-check hook confirms an edited file is syntactically sound; nothing
# confirms the documentation was regenerated after it. This check is pure git,
# no R process, so it costs nothing.
#
# Behaviour: on the first Stop with drift, block once with a reason so the
# reminder reaches Claude; on the re-entry (stop_hook_active = true) always
# exit 0, so it cannot loop. Only uncommitted changes are inspected: a commit
# is the user's statement that the tree is as they want it.
set -uo pipefail
input=$(cat)
active=$(printf '%s' "$input" | jq -r '.stop_hook_active // false' 2>/dev/null)
[ "$active" = "true" ] && exit 0

cd "${CLAUDE_PROJECT_DIR:-.}" 2>/dev/null || exit 0
git rev-parse --is-inside-work-tree >/dev/null 2>&1 || exit 0

# R/ files with uncommitted changes to roxygen lines (#' ...), staged or not
changed=$(git diff HEAD --unified=0 -- 'R/*.R' 2>/dev/null \
  | awk '/^\+\+\+ b\//{f=substr($0,7)} /^[-+][[:space:]]*#'"'"'/{print f}' | sort -u)
[ -z "$changed" ] && exit 0

docs=$(git diff HEAD --name-only -- NAMESPACE man/ 2>/dev/null)
[ -n "$docs" ] && exit 0

files=$(printf '%s\n' "$changed" | tr '\n' ' ')
reason="Roxygen comments changed in ${files}but NAMESPACE and man/ are unchanged: run \`Rscript -e 'devtools::document()'\` (R CMD check in CI fails on stale Rd/NAMESPACE), or say why it is not needed, then stop."
jq -n --arg r "$reason" '{decision: "block", reason: $r}'
exit 0
