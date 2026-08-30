#!/usr/bin/env bash
# Run the LOCAL share of the ARCHS4 coverage study as an N-wide pool
# (CLAMP_CONCURRENCY, default 3), small->large. Idempotent + resume-safe
# (run_one_coverage.sh skips finished stages and clears stale locks), so
# re-launching continues where it left off.
#
# WHY N-wide (with the BLAS cap in run_one_coverage.sh, not taskset):
#   * The pipeline is dominated by step 3 (gene correlations) whose wall time is
#     set by the 1-2 biggest chromosomes -- single-threaded under the BLAS cap,
#     so ONE model leaves this 24-core box ~90% idle for hours (a clean rs10
#     took ~9h, mostly a single busy core).
#   * `taskset` does NOT confine phenoplier's workers (they reset their own
#     affinity), so it can't be used to partition the box.
#   * Instead cap BLAS/OpenMP to 1 (run_one_coverage.sh) so nothing explodes to
#     all cores, and overlap several models so one's step-3 tail / step 1-2 runs
#     while another is in its parallel phases. Worst case is a brief single-
#     threaded oversubscription during step-7 overlaps -- work-conserving, not
#     the load-100 thrash that uncapped BLAS caused.
set -o pipefail   # not -e: one failed/locked model must not halt the pool; not -u: empty assoc-array count trips old bash
cd "$(dirname "$0")"

CONC="${CLAMP_CONCURRENCY:-3}"
# rs:seed, small->large. rs100/seed1 is byte-identical to the final archs4 and
# is reused (not listed). Alpine runs rs1/rs5 (the 6 smallest).
WORKLIST=(${CLAMP_LOCAL_WORKLIST:-10:1 10:2 10:3 25:1 25:2 25:3 50:1 50:2 50:3 75:1 75:2 75:3 100:2 100:3})

echo "[worklist ${CONC}-wide] ${#WORKLIST[@]} models, $(date -Is)"
declare -A pids
idx=0
while (( idx < ${#WORKLIST[@]} )) || (( ${#pids[@]} > 0 )); do
  while (( ${#pids[@]} < CONC )) && (( idx < ${#WORKLIST[@]} )); do
    p="${WORKLIST[$idx]}"; idx=$((idx+1))
    echo "[launch] rs${p%%:*} seed${p##*:}  $(date -Is)"
    bash run_one_coverage.sh "${p%%:*}" "${p##*:}" &
    pids[$!]=1
  done
  wait -n 2>/dev/null || true
  for pid in "${!pids[@]}"; do kill -0 "$pid" 2>/dev/null || unset 'pids[$pid]'; done
done
echo "[all-local-coverage-done] $(date -Is)"
