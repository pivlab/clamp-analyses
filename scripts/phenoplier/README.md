# phenoplier-cli integration — LV–trait association (GLS)

Runs `phenoplier shortcut gls` (GLS module–trait regression) against the CLAMP
models this repo produces, and rolls the results into cross-model
trait-association tables and HDF5 study stores.

phenoplier-cli needs `snakemake>=9` / Python ≥3.12, which is incompatible with
this repo's own `snakemake=8` / Python 3.11 (`envs/snakemake.yaml`). So it runs
in its **own conda env** (`phenoplier-cli-neo`) and each GLS run is a single
shell step, *not* part of this repo's Snakemake DAG's compute.

## Layout

| path | what |
|---|---|
| `setup_env.sh` | one-time: build the `phenoplier-cli-neo` env (installs phenoplier-cli `main` + rpy2/R) |
| `run_gls.sh` | run one model through `shortcut gls`, copy its summary back (driven by the Snakemake rules) |
| `store_build.sh` | composite one model's `.rds` + GLS results into a one-file HDF5 study store (#70) |
| `aggregate_traits.{py,sh}` | concatenate all per-model summaries into one long cross-model table |
| `final_models/` | the three **final production models** (register → GLS → store). See its README. |
| `../../workflow/rules/phenoplier.smk` | the rule module (`phenoplier_traits`, `phenoplier_final`, per-group targets) |
| `../../workflow/config/phenoplier.yaml` | config: env name, trait filter, cluster mapping, model paths |

## One-time setup (per machine/cluster)

```bash
bash scripts/phenoplier/setup_env.sh            # builds phenoplier-cli-neo
conda activate phenoplier-cli-neo
phenoplier workspace init
phenoplier workspace link /path/to/phenoplier_full_data
```

`setup_env.sh` tracks phenoplier-cli `main`. On an air-gapped cluster (e.g. pico
has no PyPI), install the env from wheels transferred from a networked host —
see `final_models/README.md`.

## Reproduce

**A — the three final models** (what PR #28 item 0 ran; the validated path):
see [`final_models/README.md`](final_models/README.md). In short, on pico:

```bash
cd scripts/phenoplier/final_models
REG=$(sbatch --parsable 00_register_models.sbatch)
for ds in archs4 gtex recount2; do
  sbatch --dependency=afterok:$REG 01_run_final_gls.sbatch $ds
done
sbatch 02_build_stores.sbatch
```

Inspect the results: `nbs/03_model_biology/00_phenoplier/00_final_models_inspection.ipynb`.

**B — the full model matrix** (ARCHS4 final/coverage/saturation, GTEx, plus the
finals) via the Snakemake rule module, from this repo's snakemake-8 env:

```bash
snakemake phenoplier_final           # the 3 final models
snakemake phenoplier_traits          # whole matrix + the aggregated long table
snakemake phenoplier_final_stores    # one HDF5 store per final model
```

> **Note:** the SLURM path of `run_gls.sh` (`--executor slurm --cluster …`)
> needs [phenoplier-cli#101](https://github.com/pivlab/phenoplier-cli/pull/101)
> (adds `--cluster` passthrough to `shortcut gls`). Until it merges, use the
> `local` target, or the `final_models/` launchers in (A), which run the local
> executor inside one sbatch allocation.

## Customizing cluster parameters

- **Which cluster the rules submit to** — `workflow/config/phenoplier.yaml: target`
  (or `--config phenoplier_target=<name>`). It maps to entries under `clusters:`:
  `local` (no SLURM), `server_cu` → phenoplier-cli's built-in `cu-pico` profile,
  `alpine` → `cu-alpine`. Account/partition/QOS live *inside* phenoplier-cli's
  profiles (`phenoplier/workflow/clusters.py`), not here.
- **Per-job resources** for the `phenoplier_gls` / `phenoplier_store` rules —
  their `resources:` / `threads:` in `phenoplier.smk`.
- **The final-models launchers** — pass `sbatch` flags per cluster and override
  paths via env vars (`PHENOPLIER_HOME`, `PHENOPLIER_CONDA_SH`, `PHENOPLIER_ENV`,
  `CLAMP_FINAL_MODELS`, `CLAMP_STORE_DIR`, …). `final_models/README.md` has pico
  and Alpine examples.
- **Trait filter** — `phenoplier.yaml: trait_filter` (`biomedical` | `all`).
  Filtering shrinks the Benjamini–Hochberg test count, so keep it fixed across a
  comparison (a coverage/saturation sweep, or these three models).
