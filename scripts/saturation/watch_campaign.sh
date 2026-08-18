#!/usr/bin/env bash
#
# Watchdog for the saturation campaign.
#
# A Snakemake driver on a cluster reports failures only in its own log, which
# nobody reads until something looks wrong.  This polls the real signals --
# failed Slurm jobs, R errors in fit logs, published models -- and writes a
# one-line status plus an explicit ALERT when jobs start failing, so a broken
# run is obvious within one poll instead of hours later.
#
# Usage: scripts/saturation/watch_campaign.sh <repo_root> [poll_seconds]

set -uo pipefail

ROOT="${1:?repo root required}"
POLL="${2:-300}"
SAT="$ROOT/output/01_model_building/02_archs4/03_saturation_study"
STATUS="$ROOT/output/03_model_biology/02_archs4/01_saturation/STATUS.txt"
mkdir -p "$(dirname "$STATUS")"
SINCE="$(date '+%Y-%m-%dT%H:%M:%S')"
echo "[$SINCE] watchdog started; counting failures from here" >> "$STATUS"

while :; do
    now=$(date '+%Y-%m-%d %H:%M:%S')
    models=$(find "$SAT/models" -name CLAMPfull_bp -type d 2>/dev/null | wc -l)
    full_ora=$(find "$SAT/ora" -path '*CLAMPfull*' -name complete 2>/dev/null | wc -l)
    base_ora=$(find "$SAT/ora" -path '*CLAMPbase*' -name complete 2>/dev/null | wc -l)
    queued=$(squeue -u "$USER" -h 2>/dev/null | wc -l)

    # Failures since this watchdog started, not a rolling day. Counting the
    # whole day mixed 22 already-fixed failures in with 1 new one and buried
    # the signal, which is how a live OOM went unnoticed.
    # Count job rows only. sacct also emits a row per step (.batch, .0,
    # .extern), which reported a single failed job as three.
    sacct_jobs() {
        sacct -u "$USER" -S "$SINCE" --format="$1" --noheader 2>/dev/null \
            | grep -vE '^[[:space:]]*[0-9_]+\.'
    }
    failed=$(sacct_jobs JobID,State | grep -cE 'FAILED|OUT_OF_MEMORY|TIMEOUT' || true)
    oom=$(sacct_jobs JobID,State | grep -c 'OUT_OF_MEMORY' || true)

    # R-level errors sitting in fit logs, which Slurm reports only as exit 1.
    r_errors=$(grep -rilE '^Error|invalid filename' "$SAT" \
                    --include=fit.log --include=clampbase.log 2>/dev/null | wc -l)

    # Our own share of the shared node. Snakemake's resource ledger is per
    # driver, so two coexisting drivers (e.g. a tmux restart that left the old
    # driver's jobs running) each stay under the cap while the sum blows past
    # it. That happened once and went unnoticed; only the total is meaningful.
    my_cpu=$(squeue -u "$USER" -h -t RUNNING -o '%C' 2>/dev/null | paste -sd+ | bc 2>/dev/null || echo 0)
    my_mem=$(squeue -u "$USER" -h -t RUNNING -o '%m' 2>/dev/null \
             | sed 's/M$//;s/G$/*1024/' | paste -sd+ | bc 2>/dev/null || echo 0)
    my_cpu=${my_cpu:-0}; my_mem=${my_mem:-0}

    line="[$now] models $models/107  fullORA $full_ora/321  baseORA $base_ora/321  queued $queued  failed $failed (oom $oom)  fitLogErrors $r_errors  cpu ${my_cpu}/126  mem ${my_mem%%.*}/675000"
    echo "$line" >> "$STATUS"

    if (( $(echo "$my_cpu > 126" | bc -l) )) || (( $(echo "$my_mem > 675000" | bc -l) )); then
        {
            echo "  ALERT: over the 75% node cap (cpu ${my_cpu}/126, mem ${my_mem%%.*}/675000)."
            echo "  Likely a second Snakemake driver's jobs still running. Check submit times:"
            squeue -u "$USER" -h -o '    %.8i %.9T %.9m %.20V' 2>/dev/null | sort -k5 | head -8
        } >> "$STATUS"
    fi

    if [[ "$failed" -gt 0 || "$r_errors" -gt 0 ]]; then
        {
            echo "  ALERT: $failed failed jobs since start ($oom OOM), $r_errors fit logs with R errors."
            sacct -u "$USER" -S "$SINCE" --format=JobID%12,JobName%20,State%16,ReqMem,MaxRSS --noheader 2>/dev/null \
                | grep -E 'FAILED|OUT_OF_MEMORY|TIMEOUT' | head -8 | sed 's/^/    /' 
            echo "  Recent fit-log errors:"
            grep -rhoE '^Error.*' "$SAT" \
                --include=fit.log --include=clampbase.log 2>/dev/null | sort -u | head -5 | sed 's/^/    /'
        } >> "$STATUS"
    fi

    # Stop once every model and every ORA has landed.
    if [[ "$models" -ge 107 && "$full_ora" -ge 321 && "$base_ora" -ge 321 ]]; then
        echo "[$now] COMPLETE" >> "$STATUS"
        exit 0
    fi
    sleep "$POLL"
done
