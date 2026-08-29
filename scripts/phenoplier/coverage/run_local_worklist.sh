#!/usr/bin/env bash
# Run the LOCAL share of the ARCHS4 coverage study, ONE model at a time (each
# uses the whole box). Ordered small->large so the coverage curve fills in from
# the fast end. Idempotent (run_one_coverage.sh skips finished stages), so
# re-launching resumes.
#
# WHY 1-wide (not 2 concurrent): a workstation has no delegated cpuset cgroup,
# and phenoplier's step-3/6/7 worker subprocesses reset their own CPU affinity,
# so `taskset` does NOT hold them -- two concurrent runs then oversubscribe the
# box (observed load >100 on 24 cores, everything crawls). Running one model at
# a time gives step 6/7 the full core count and avoids oversubscription entirely.
set -uo pipefail   # not -e: one failed/locked model must not halt the worklist
cd "$(dirname "$0")"

# rs:seed, small->large. rs100/seed1 is byte-identical to the final archs4 and
# is reused (not listed). Alpine runs rs1/rs5 (the 6 smallest).
WORKLIST=(${CLAMP_LOCAL_WORKLIST:-10:1 10:2 10:3 25:1 25:2 25:3 50:1 50:2 50:3 75:1 75:2 75:3 100:2 100:3})

echo "[worklist-1wide] ${#WORKLIST[@]} models, $(date -Is)"
for p in "${WORKLIST[@]}"; do
  echo "[run] rs${p%%:*} seed${p##*:}  $(date -Is)"
  bash run_one_coverage.sh "${p%%:*}" "${p##*:}" || echo "[warn] rs${p%%:*} seed${p##*:} exited non-zero (see its log); continuing"
done
echo "[all-local-coverage-done] $(date -Is)"
