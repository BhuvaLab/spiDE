#!/usr/bin/env bash
# PreToolUse (Edit|Write|Bash): refuse in-place rewrites of research drivers
# and sbatch scripts while this user has SLURM jobs running or queued.
#
# WHY (CLAUDE.md, "A warning about the benchmark harness"): R parses a script
# incrementally through an 8 KB stdio buffer, so a task that is inside its
# multi-hour fit reads the *next* expression from whatever the file holds by
# then, at the old byte offset. On 2026-09-07 an in-place edit killed 255
# benchmark tasks after they had written their result, and then eight cohort
# tasks six to seven hours into their fit. The Edit and Write tools, and a
# shell `>` / `>>` / `cp` / `tee`, rewrite the EXISTING inode -- the dangerous
# kind. `git checkout`, `git show > newfile && mv`, and `sed -i` create a new
# inode and cannot reach a running reader.
#
# Scope: research/**/*.R, research/**/*.sbatch, research/**/*.sh, and any path
# under a drivers/ directory (the frozen copies freeze_snapshot.sh places in
# $SPIDE_PKG/drivers/). The interactive session's own job (SLURM_JOB_ID) and
# VS Code / RStudio / Jupyter server jobs do not count as "running".
#
# Escape hatch: SPIDE_ALLOW_DRIVER_EDIT=1 in the environment, for the case
# where you have checked the running jobs do not read the file.
#
# Exit 2 = block the call and return this message to Claude.
set -uo pipefail
[ "${SPIDE_ALLOW_DRIVER_EDIT:-0}" = 1 ] && exit 0
command -v squeue >/dev/null 2>&1 || exit 0

input=$(cat)
tool=$(printf '%s' "$input" | jq -r '.tool_name // empty' 2>/dev/null)

# path pattern for guarded files (extended regex, applied to a path string)
guard_re='(^|/)research/.*\.(R|r|sbatch|sh)$|(^|/)drivers/[^/[:space:]]+$'

target=""
case "$tool" in
  Edit|Write)
    f=$(printf '%s' "$input" | jq -r '.tool_input.file_path // .tool_input.filePath // empty' 2>/dev/null)
    [ -n "$f" ] && printf '%s' "$f" | grep -Eq "$guard_re" && target="$f"
    ;;
  Bash)
    cmd=$(printf '%s' "$input" | jq -r '.tool_input.command // empty' 2>/dev/null)
    [ -z "$cmd" ] && exit 0
    # in-place writers: a redirect into a guarded path, tee onto it, or cp over it
    redir_re='>>?[[:space:]]*"?(([^[:space:]"|;&]*/)?research/[^[:space:]"|;&]*\.(R|r|sbatch|sh)|[^[:space:]"|;&]*/drivers/[^[:space:]"|;&/]+)([[:space:]]|$|["|;&])'
    if printf '%s' "$cmd" | grep -Eq "$redir_re"; then
      target=$(printf '%s' "$cmd" | grep -Eo "$redir_re" | head -1 | sed -E 's/^>>?[[:space:]]*"?//; s/[[:space:]"|;&]*$//')
    else
      # no `^` inside a mid-pattern group: ugrep (the grep on Bunya) does not
      # match it there, GNU grep does; both accept the forms below
      path_re='"?(\./|[^[:space:]"|;&]*/)?research/[^[:space:]"|;&]*\.(R|r|sbatch|sh)|[^[:space:]"|;&]*/drivers/[^[:space:]"|;&/]+'
      write_re="(^|[[:space:]|;&])(tee|cp)([[:space:]]+[^[:space:]|;&]+)*[[:space:]]+($path_re)([[:space:]]|\$|[\"|;&])"
      if printf '%s' "$cmd" | grep -Eq "$write_re"; then
        target=$(printf '%s' "$cmd" | grep -Eo "$path_re" | tail -1 | sed 's/^"//')
      fi
    fi
    ;;
  *) exit 0 ;;
esac
[ -z "$target" ] && exit 0

# running or pending jobs, excluding this session and interactive servers
self="${SLURM_JOB_ID:-}"
jobs=$(squeue -u "$USER" -h -t R,PD,CG -o '%i %T %j' 2>/dev/null \
  | grep -Ev '^([0-9]+)_?[0-9]*[^ ]* [A-Z]+ (vscc|vscode|code-server|rstudio|jupyter|interactive|bash|sys/dashboard)' \
  | { if [ -n "$self" ]; then grep -Ev "^${self}(_|\.| )"; else cat; fi; })
[ -z "$jobs" ] && exit 0

n=$(printf '%s\n' "$jobs" | wc -l)
cat >&2 <<MSG
BLOCKED: in-place write to a driver while $n SLURM job(s) are running or queued.

  target: $target

$(printf '%s\n' "$jobs" | head -8 | sed 's/^/  /')

R parses a driver incrementally (8 KB buffer). A task inside its fit reads
the NEXT expression from the file as it is now, at the old byte offset --
this killed 255 benchmark tasks and eight 6-7 h cohort tasks on 2026-09-07.
Edit, Write, '>' and 'cp' rewrite the existing inode and reach the reader.

Safe alternatives:
  - wait until the array finishes (squeue -u \$USER), then edit;
  - write a NEW inode: edit a copy and 'mv' it over, or 'git checkout' /
    'git show <rev>:path > path.new && mv path.new path' -- a rename swaps
    the directory entry and the running reader keeps its old inode;
  - if the jobs demonstrably do not read this file, re-run with
    SPIDE_ALLOW_DRIVER_EDIT=1 in the environment.
MSG
exit 2
