# Setup

This project runs in two ways:

- **Ad hoc / cluster jobs** (`jobs/*.sh`): SLURM `sbatch` scripts that run `jupyter nbconvert --execute` on a single notebook directly inside the `clamp-analyses` (or `gpu-kmeans`) conda env.
- **Snakemake pipeline** (`workflow/Snakefile`): orchestrates the same scripts/notebooks as a DAG and tracks which outputs are up to date. This is the preferred way to (re)build a model end-to-end.

Both paths use the same two analysis conda envs; Snakemake additionally needs its own minimal env to run the orchestrator itself.

## 1. Analysis environments

Create `clamp-analyses` (CPU R/Python) and `gpu-kmeans` (GPU clustering). Full instructions are in the root [README.md](../../README.md#-install-dependencies) — do that first if you haven't.

> **Known issue:** `envs/clamp-analyses.yaml` and `config.R` currently contain unresolved git merge-conflict markers (`<<<<<<< Updated upstream` / `=======` / `>>>>>>> Stashed changes`). Resolve those before running `conda env create -f envs/clamp-analyses.yaml`, otherwise the YAML won't parse (and `config.R` won't source in R).

Snakemake rules reference these environments **by name** (e.g. `conda: "clamp-analyses"`), not by re-solving the `.yaml` files each run. This only works if a conda env with that exact name already exists — Snakemake will not create it for you when a name is given instead of a path to a `.yaml` file.

## 2. Snakemake itself

Keep Snakemake in its own env, separate from the analysis envs:

```bash
conda create -n snakemake -c bioconda -c conda-forge snakemake=8.4.11
conda activate snakemake
```

Any recent Snakemake 8.x should work; `8.4.11` is what's currently installed locally.

## 3. Running the pipeline

From the repo root, with the `snakemake` env active:

```bash
# dry run: see what would execute, without doing anything
snakemake -n

# build everything (each rule uses its own `conda:` env, not the snakemake env)
snakemake --use-conda --cores <N>

# build just one target, e.g. the GTEx models or the pseudobulk models
snakemake --use-conda --cores <N> full_models_gtex
snakemake --use-conda --cores <N> full_models_pseudobulk
```

`--use-conda` is required — without it, rules run inside whatever env `snakemake` itself was launched from and will fail on missing R/Python packages.

Some rules fetch external data over HTTP the first time they run (e.g. `download_gtex_raw`, `pathway_prior`). See [data/README.md](../../data/README.md) for what's downloaded automatically vs. what needs to be placed manually before running the pipeline.

## 4. Cluster execution

Rules declare `resources: mem_mb=..., runtime=...`, which is the interface Snakemake executor plugins (e.g. the SLURM plugin) read to submit each rule as its own cluster job. No such plugin or `--profile` is configured in this repo yet — today, `jobs/*.sh` submits `sbatch` scripts that call `jupyter nbconvert` directly rather than going through Snakemake's own SLURM integration. If that gets wired up later, document the profile/setup here.
