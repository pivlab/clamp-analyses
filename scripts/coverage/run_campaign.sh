#!/usr/bin/env bash
#
# Start, or restart, the host half of the GO:BP coverage campaign.
#
# The campaign is long enough (days) that its drivers must outlive the shell,
# editor or agent session that started them.  Every component here is launched
# with setsid so it is reparented to init, and every component is idempotent:
# running this script while the campaign is healthy starts nothing.
#
# Snakemake itself decides what still needs doing, so nothing already finished
# is recomputed.  The high-memory ARCHS4 fits are not started here; they run on
# a Slurm cluster, and this script only reports on them.
#
# Site-specific values (SSH alias, cluster paths, environment locations) are
# read from a site environment file that is deliberately not committed.  See
# "Running ARCHS4" in README.md for the template.
#
# Usage:
#   scripts/coverage/run_campaign.sh [--site-env FILE] [--status] [--dry-run]

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

SITE_ENV="${CLAMP_SITE_ENV:-.campaign.env}"
STATUS_ONLY=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --site-env) SITE_ENV="$2"; shift 2 ;;
    --status)   STATUS_ONLY=1; shift ;;
    --dry-run)  DRY_RUN=1; shift ;;
    -h|--help)  sed -n '2,20p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

if [[ -f "$SITE_ENV" ]]; then
  # shellcheck disable=SC1090
  source "$SITE_ENV"
else
  echo "No site environment file at '$SITE_ENV'." >&2
  echo "Copy the template from README.md ('Running ARCHS4') and fill it in." >&2
  exit 2
fi

: "${CLAMP_SNAKEMAKE_BIN:?set CLAMP_SNAKEMAKE_BIN in $SITE_ENV}"
: "${CLAMP_ANALYSES_BIN:?set CLAMP_ANALYSES_BIN in $SITE_ENV}"
: "${CLAMP_REMOTE:?set CLAMP_REMOTE in $SITE_ENV}"
: "${CLAMP_REMOTE_ROOT:?set CLAMP_REMOTE_ROOT in $SITE_ENV}"
: "${CLAMP_REMOTE_PROFILE:?set CLAMP_REMOTE_PROFILE in $SITE_ENV}"

# Host memory the campaign may schedule in total, and how much of it a single
# comparator fit may claim.  Lower CLAMP_FIT_SLOTS to 1 while another workload
# holds a large share of host memory.
CLAMP_TOTAL_MEM_MB="${CLAMP_TOTAL_MEM_MB:-96000}"
CLAMP_FIT_MEM_MB="${CLAMP_FIT_MEM_MB:-64000}"
CLAMP_ORA_MEM_MB="${CLAMP_ORA_MEM_MB:-64000}"
CLAMP_FIT_SLOTS="${CLAMP_FIT_SLOTS:-2}"
CLAMP_ORA_SLOTS="${CLAMP_ORA_SLOTS:-8}"
CLAMP_POLL_SECONDS="${CLAMP_POLL_SECONDS:-300}"

export PATH="$CLAMP_SNAKEMAKE_BIN:$CLAMP_ANALYSES_BIN:$PATH"
export PYTHONUNBUFFERED=1

COVERAGE_DIR="output/03_model_biology/02_archs4/00_coverage"
MODEL_DIR="$COVERAGE_DIR/models"
ORA_DIR="$COVERAGE_DIR/ora"
mkdir -p "$COVERAGE_DIR"

SNAKE=(snakemake -s workflow/Snakefile --profile workflow/profiles/local --nolock)

# Total published ORA directories the campaign must produce.  Derived from the
# config rather than hardcoded, so adding a fraction, seed or database does not
# silently leave the report waiting on a stale number.
expected_ora_count() {
  python3 - <<'PY'
import yaml
config = yaml.safe_load(open("workflow/config/archs4.yaml"))["archs4"]
coverage = config["coverage"]
levels = len(coverage["levels"])
seeds = len(coverage["seeds"])
databases = len(config["ora"]["databases"])
comparators = len(coverage["comparators"])
archs4 = levels * seeds * databases * 2                 # CLAMPfull and CLAMPbase
comparator_full = comparators * seeds * databases
comparator_base = comparators * databases               # one shared base per dataset
print(archs4 + comparator_full + comparator_base)
PY
}

published_ora_count() {
  find "$ORA_DIR" -type f -name complete 2>/dev/null | wc -l
}

# The comparator models the ORA stage has to wait for, also config-derived.
comparator_validated_files() {
  python3 - "$MODEL_DIR" <<'PY'
import sys, yaml
model_root = sys.argv[1]
coverage = yaml.safe_load(open("workflow/config/archs4.yaml"))["archs4"]["coverage"]
print(" ".join(
    f"{model_root}/{dataset}/rs100/seed{seed}/validated.json"
    for dataset in coverage["comparators"]
    for seed in coverage["seeds"]
))
PY
}

validated_models() {
  find "$MODEL_DIR" -name validated.json 2>/dev/null | wc -l
}

# A component is identified by a stable substring of its command line.
component_pid() {
  pgrep -u "$USER" -f "$1" | head -1
}

launch() {
  local name="$1" pattern="$2" logfile="$3" script="$4"
  local existing
  existing="$(component_pid "$pattern" || true)"
  if [[ -n "$existing" ]]; then
    printf '  %-22s already running (pid %s)\n' "$name" "$existing"
    return 0
  fi
  if [[ "$DRY_RUN" == 1 ]]; then
    printf '  %-22s would start\n' "$name"
    return 0
  fi
  setsid bash -c "$script" </dev/null >>"$logfile" 2>&1 &
  disown || true
  sleep 1
  printf '  %-22s started (log: %s)\n' "$name" "$logfile"
}

EXPECTED_ORA="$(expected_ora_count)"
COMPARATOR_VALIDATED="$(comparator_validated_files)"

echo "Coverage campaign status"
echo "  validated models        $(validated_models)"
echo "  published ORA dirs      $(published_ora_count) / $EXPECTED_ORA"
if ssh -o BatchMode=yes -o ConnectTimeout=10 "$CLAMP_REMOTE" true 2>/dev/null; then
  echo "  cluster jobs            $(ssh -o BatchMode=yes "$CLAMP_REMOTE" 'squeue -u $USER -h | wc -l' 2>/dev/null | tr -d '[:space:]')"
  echo "  cluster drivers         $(ssh -o BatchMode=yes "$CLAMP_REMOTE" 'tmux ls 2>/dev/null | wc -l' 2>/dev/null | tr -d '[:space:]') tmux session(s)"
else
  echo "  cluster                 unreachable (host half continues regardless)"
fi

if [[ "$STATUS_ONLY" == 1 ]]; then
  exit 0
fi

echo "Host components"

# 1. Comparator model fits (GTEx and recount2, three seeds each).
launch "comparator fits" \
  "coverage_bp_comparator_models" \
  "$COVERAGE_DIR/local-driver.log" \
  "${SNAKE[*]} coverage_bp_comparator_models \
     --resources mem_mb=$CLAMP_TOTAL_MEM_MB fit_mem_mb=$CLAMP_FIT_MEM_MB \
     ora_mem_mb=$CLAMP_ORA_MEM_MB fit_slots=$CLAMP_FIT_SLOTS \
     ora_slots=$CLAMP_ORA_SLOTS large_fit_slot=0"

# 2. Comparator ORAs, once every comparator model has been validated.
launch "comparator ORA" \
  "coverage_bp_comparator_ora" \
  "$COVERAGE_DIR/comparator-ora-driver.log" \
  "while :; do
     ready=1
     for marker in $COMPARATOR_VALIDATED; do
       test -f \"\$marker\" || ready=0
     done
     test \"\$ready\" -eq 1 && break
     sleep $CLAMP_POLL_SECONDS
   done
   ${SNAKE[*]} coverage_bp_comparator_ora --cores 8 \
     --resources mem_mb=$CLAMP_ORA_MEM_MB fit_mem_mb=0 ora_mem_mb=$CLAMP_ORA_MEM_MB \
     fit_slots=0 ora_slots=$CLAMP_ORA_SLOTS large_fit_slot=0"

# 3. Dynamic ORA router: cluster when it has spare capacity, host otherwise.
launch "ORA router" \
  "scripts/coverage/route_ora.py" \
  "$COVERAGE_DIR/router.log" \
  "python scripts/coverage/route_ora.py \
     --remote $CLAMP_REMOTE --remote-root $CLAMP_REMOTE_ROOT \
     --remote-profile $CLAMP_REMOTE_PROFILE \
     ${CLAMP_REMOTE_PATH_DIRS:-} \
     --local-total-mem-mb $CLAMP_TOTAL_MEM_MB \
     --poll-seconds $CLAMP_POLL_SECONDS"

# 4. Aggregation and the coverage notebook, once every ORA has been published.
launch "aggregation + report" \
  "coverage_report_archs4" \
  "$COVERAGE_DIR/report-driver.log" \
  "while :; do
     completed=\$(find $ORA_DIR -type f -name complete 2>/dev/null | wc -l)
     test \"\$completed\" -ge $EXPECTED_ORA && break
     sleep $CLAMP_POLL_SECONDS
   done
   ${SNAKE[*]} aggregate_bp_coverage coverage_report_archs4 \
     --forcerun aggregate_bp_coverage coverage_report_archs4 \
     --cores 2 --resources mem_mb=16000 fit_mem_mb=0 ora_mem_mb=0 \
     fit_slots=0 ora_slots=0 large_fit_slot=0"

echo
echo "Re-run this script any time; it restarts only what is not running."
echo "Cluster side: ssh \$CLAMP_REMOTE 'tmux ls; squeue -u \$USER'"
