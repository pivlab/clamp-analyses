# ARCHS4 coverage study — GLS run (pico)

Runs `phenoplier shortcut gls` for the ARCHS4 **coverage study** models — how the
LV–trait signal grows as more of the compendium is used (coverage fraction rs).

Models (on pico):
`…/02_coverage_study/models/archs4/rs{1,5,10,25,50,75,100}/seed{1,2,3}/CLAMPfull_bp/CLAMPfull_bp.rds`
— 7 fractions × 3 seeds = **21**. `rs100/seed1` is byte-identical to the final
archs4 (already run) → **20 new**. (Distinct from the saturation study in
`../saturation_k1728/`, which varies K.) Every run uses
**`--trait-filter biomedical`** (2,366 of 4,049 traits), same as finals/saturation.

## How to reproduce (pico)

Two throttled SLURM arrays, memory tiered by fraction. `pico_coverage.sbatch` runs
the big models on `defq` (`-c 32`, 150 G) `%2`-wide so they backfill around other
work; `pico_coverage_tiny.sbatch` runs the small fractions in a wider, lower-memory
array. Both are idempotent + resume-safe (they skip finished stages), and the models
already live on pico's filesystem, so there is no transfer.

```bash
cd scripts/phenoplier/coverage
sbatch --array=0-6%2 pico_coverage.sbatch        # big models (rs50 s2/s3, rs75×3, rs100 s2/s3)
sbatch pico_coverage_tiny.sbatch                 # small fractions (see its header for the array spec)
```

phenoplier sizes its step-6/7 process pools from the SLURM cgroup (`-c N`), which
both caps and pins them — so on pico there is nothing extra to set. All paths are
env-overridable (`PHENOPLIER_HOME`, `CLAMP_COV_DIR`, the model source root, …); see
each script's header.

> **As-run note.** The 21 models were actually produced across pico + a workstation
> (`local`) + `alpine` in parallel for throughput (pico was busy on the saturation
> run at the time). The off-pico launchers — `run_local_worklist.sh`,
> `run_one_coverage.sh`, `alpine_coverage.sbatch`, and their thread-capping notes
> for un-cgrouped hosts — are preserved verbatim on the archival branch
> `phenoplier-asrun-launchers`
> (`81de511c6ee22d5fbaff1ef107ac613a1ce5bfb9`). This directory ships only the
> pico reproduction path; the results are identical regardless of which host ran a
> given model (see `../RUN_SUMMARY.md`).

## Outputs
- GLS summaries: `<workspace>/projects/cov_rs<f>_seed<s>_CLAMPfull_bp/results/gls/phenoplier/`
- HDF5 stores: `<coverage-dir>/stores/cov_rs<f>_seed<s>_CLAMPfull_bp.h5`
- Inspect: `nbs/03_model_biology/00_phenoplier/02_coverage.ipynb` (recovery vs coverage).
