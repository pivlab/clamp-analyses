# ARCHS4 saturation (k=1728) GLS run

Runs `phenoplier shortcut gls` for the ARCHS4 **saturation-study** models at fixed
**k = 1728** LVs across coverage fractions and seeds, to trace how many LV–trait
associations are recovered as more of the compendium is used.

Models (on pico):
`…/output/01_model_building/02_archs4/03_saturation_study/models/rs{1,5,10,25,50,75}/k1728/seed{1,2,3}/CLAMPfull_bp/CLAMPfull_bp.rds`

**17 models** (not 18 — `rs1/seed1` is missing on pico). `.rds` size scales with
coverage (rs1 ~130 MB → rs75 ~5.7 GB), and so does GLS wall-time (steps 1–5 scale
with the model's gene count).

## How it was run (pico)

Three ordered steps; the two heavy ones are **SLURM job arrays throttled to 6
concurrent** so we stay under **80 % of pico** (defq = 168 cores / 900 GB → cap
134 cores / 720 GB). At `-c 22`, 6 concurrent = **132 cores (78.6 %) / 600 GB (67 %)**.

```bash
cd scripts/phenoplier/saturation_k1728
REG=$(sbatch --parsable 00_register_models.sbatch)                       # register 17 once
GLS=$(sbatch --parsable --dependency=afterok:$REG --array=0-16%6 01_run_saturation_gls.sbatch)
sbatch --dependency=afterok:$GLS --array=0-16%6 02_build_stores.sbatch   # one HDF5 store per model
```

- Register once → run each model with `shortcut gls --model` (read-only on the
  shared `registry.toml`; also avoids re-extracting each up-to-5.7 GB `.rds` per job).
- Array indices are **big-model-first** (rs75 → rs1), so the three long rs75 jobs
  start in the first wave and the smaller ones backfill. `-c 22` covers step 3's 22
  chromosomes in one wave; phenoplier sizes its step-6/7 pools from the cgroup.
- Every run uses **`--trait-filter biomedical`** (2,366 of 4,049 traits), same as the
  final-models run, so the saturation curve is consistent across fractions.

Estimated cost: **35 SLURM allocations** (1 register + 17 GLS + 17 store tasks),
**~14–16 h** wall-clock for GLS (bounded by the 3 rs75 jobs) + ~30–60 min stores,
~2,300 core-hours.

Results land in `<workspace>/projects/sat_rs<f>_k1728_seed<s>_CLAMPfull_bp/`; HDF5
stores in `<workspace>/clamp_saturation_stores/`. Inspect them with
`nbs/03_model_biology/00_phenoplier/01_saturation_k1728.ipynb`.

## Customizing

Everything is env-overridable (defaults = pico): `PHENOPLIER_CONDA_SH`,
`PHENOPLIER_ENV`, `PHENOPLIER_HOME`, `CLAMP_SAT_ROOT`, `CLAMP_SAT_STORE_DIR`,
`CLAMP_TRAIT_FILTER`. To change the concurrency cap, edit the `%6` in `--array`
and/or `-c` (keep `-c × concurrency` ≤ 80 % of the node). See
`../final_models/README.md` for the equivalent single-node run and the
phenoplier-cli caveats.
