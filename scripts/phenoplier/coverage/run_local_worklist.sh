#!/usr/bin/env bash
# Run the LOCAL share of the ARCHS4 coverage study: a worklist of (rs, seed)
# models, 2 at a time, each pinned to a disjoint set of physical cores.
#
# WHY the pinning: there is no SLURM cgroup on a workstation, so phenoplier
# would size its step-6/7 process pools from the whole box (sched_getaffinity)
# in EACH run -- two concurrent runs then oversubscribe the machine. `taskset`
# confines each run to half the physical cores, so each sizes its pool to ~12
# workers and the two don't fight. Cores are AMD 24C/48T: physical core N =
# logical {N, N+24}, so slot 0 = physical 0-11 = "0-11,24-35", slot 1 = 12-23.
#
# Ordered small->large so the coverage curve fills in from the fast end.
# Idempotent (run_one_coverage.sh skips finished stages), so re-launching resumes.
set -euo pipefail
cd "$(dirname "$0")"

# rs:seed, small->large. rs100/seed1 is byte-identical to the final archs4 and
# is reused, so it is NOT in this list. Alpine takes rs1/rs5 (the 6 smallest).
WORKLIST=(${CLAMP_LOCAL_WORKLIST:-10:1 10:2 10:3 25:1 25:2 25:3 50:1 50:2 50:3 75:1 75:2 75:3 100:2 100:3})
SLOT_CORES=("${CLAMP_SLOT0_CORES:-0-11,24-35}" "${CLAMP_SLOT1_CORES:-12-23,36-47}")

echo "[worklist] ${#WORKLIST[@]} models, 2-wide, $(date -Is)"
declare -A slotpid
idx=0
while (( idx < ${#WORKLIST[@]} )) || (( ${#slotpid[@]} > 0 )); do
  for s in 0 1; do
    if [[ -z "${slotpid[$s]:-}" ]] && (( idx < ${#WORKLIST[@]} )); then
      p="${WORKLIST[$idx]}"; idx=$((idx+1))
      echo "[launch] slot $s cores ${SLOT_CORES[$s]} -> rs${p%%:*} seed${p##*:}  $(date -Is)"
      taskset -c "${SLOT_CORES[$s]}" bash run_one_coverage.sh "${p%%:*}" "${p##*:}" &
      slotpid[$s]=$!
    fi
  done
  wait -n 2>/dev/null || true
  for s in 0 1; do
    if [[ -n "${slotpid[$s]:-}" ]] && ! kill -0 "${slotpid[$s]}" 2>/dev/null; then
      wait "${slotpid[$s]}" 2>/dev/null || echo "[warn] slot $s job exited non-zero"
      unset 'slotpid[$s]'
    fi
  done
done
echo "[all-local-coverage-done] $(date -Is)"
