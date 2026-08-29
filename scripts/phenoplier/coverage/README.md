# ARCHS4 coverage study — GLS run (multi-machine)

Runs `phenoplier shortcut gls` for the ARCHS4 **coverage study** models — how the
LV–trait signal grows as more of the compendium is used (coverage fraction rs).

Models (on pico): `…/02_coverage_study/models/archs4/rs{1,5,10,25,50,75,100}/seed{1,2,3}/CLAMPfull_bp/CLAMPfull_bp.rds`
— 7 fractions × 3 seeds = **21**. `rs100/seed1` is byte-identical to the final
archs4 (already run) → **20 new**. (Distinct from the saturation study in
`03_saturation_study/`, which varies K.)

## Why multi-machine

pico was busy (and slow) on the saturation run, so the coverage study was run in
parallel off-pico, on **local** (primary) + **alpine** (backfill). Both envs were
brought to phenoplier-cli **v0.5.2**. Every run uses **`--trait-filter biomedical`**
(2,366 of 4,049 traits), same as the finals/saturation.

| where | models | how |
|---|---|---|
| **local** (24C/48T, 251G) | 14: rs10×3, rs25×3, rs50×3, rs75×3, rs100×2 | `run_local_worklist.sh` — **1-wide**, small→large |
| **alpine** (acpu) | 6 smallest: rs1×3, rs5×3 | `alpine_coverage.sbatch` — `--array=0-5%6` backfill |
| (reused) | rs100/seed1 | the finals archs4 result |

### Local: `run_local_worklist.sh` + `run_one_coverage.sh`
Runs models **one at a time**, each using the whole box (`--n-jobs 24`; step 7 gets
~24 workers). `run_one_coverage.sh` streams each model from pico, runs GLS, then
`store build`; it is idempotent (skips finished stages), so the worklist re-launches
to resume.

**Why 1-wide, not 2 concurrent:** a workstation has no *delegated* cpuset cgroup, and
phenoplier's step-3/6/7 worker subprocesses **reset their own CPU affinity**, so
`taskset` does not confine them. Two concurrent runs therefore each size their pools
from the full box and oversubscribe it (observed load >100 on 24 cores → everything
crawls). One model at a time avoids this and, since step 6/7 are the bottleneck and
scale with cores, gives comparable throughput without the thrash.

```bash
# local (env + workspace at /media/data already set up):
setsid nohup bash scripts/phenoplier/coverage/run_local_worklist.sh \
    > /media/data/clamp_coverage/logs/worklist_master.log 2>&1 < /dev/null &
```

### Alpine: `alpine_coverage.sbatch`
The 6 small models as a `%6` array; SLURM's cgroup (`-c 16`) caps the pools, so no
`taskset` there. Models are streamed pico→local→alpine first (alpine can't reach pico).

```bash
sbatch --array=0-5%6 scripts/phenoplier/coverage/alpine_coverage.sbatch
```

## Outputs
- GLS summaries: `<workspace>/projects/cov_rs<f>_seed<s>_CLAMPfull_bp/results/gls/phenoplier/`
- HDF5 stores: `<coverage-dir>/stores/cov_rs<f>_seed<s>_CLAMPfull_bp.h5`
- Inspect: `nbs/03_model_biology/00_phenoplier/02_coverage.ipynb` (recovery vs coverage).

All paths are env-overridable (`PHENOPLIER_ENVBIN`, `PHENOPLIER_HOME`, `CLAMP_COV_DIR`,
`CLAMP_COV_SRC_SSH/ROOT`, `CLAMP_LOCAL_WORKLIST`, `CLAMP_SLOT{0,1}_CORES`).
