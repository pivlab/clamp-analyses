# CLAMP analyses

## Setup

CLAMP uses three Conda environments: `clamp-analyses` (core CPU modeling/analysis),
`gpu-kmeans` (GPU clustering), and `snakemake` (pipeline orchestrator).

```bash
./setup.sh
```

Installs Conda if needed, creates all three environments, and installs the pinned
[`chikinalab/CLAMP`](https://github.com/chikinalab/CLAMP) R package into
`clamp-analyses`. Safe to re-run. Verify with `nbs/00_setup/00_check_setup.ipynb`.

> [!WARNING]
> Requires CLAMP commit
> [`818e13ba55d66840e0710c3f1ac15f6d97e1dd8b`](https://github.com/chikinalab/CLAMP/commit/818e13ba55d66840e0710c3f1ac15f6d97e1dd8b) —
> do not update it independently.

<details>
<summary>Manual setup, step by step</summary>

1. Install Conda (Miniconda or Mambaforge).
2. (Optional, for `gpu-kmeans`) verify GPU drivers with `nvidia-smi`.
3. Create the environments:
   ```bash
   conda create --name clamp-analyses --file envs/clamp-analyses.lock
   conda run -n clamp-analyses python -m pip install -r envs/clamp-analyses.pip.lock
   conda create --name gpu-kmeans --file envs/gpu-kmeans.lock
   conda env create -n snakemake -f envs/snakemake.yaml
   ```
4. Install the CLAMP R package:
   ```bash
   conda activate clamp-analyses
   Rscript scripts/install_clamp.R
   Rscript scripts/install_clamp.R --check   # verify only
   ```
</details>

Each notebook states its required environment in its first cell. Run all `snakemake`
commands from the repo root with `conda activate snakemake` and `--use-conda`:

```bash
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile <target>
```

## Pseudobulk

Benchmarks CLAMP against other latent-variable/matrix-decomposition methods on pseudobulk
profiles built from single-cell/single-nucleus RNA-seq cohorts (brain, heart, PBMC, lung),
then projects the fitted latent variables back onto individual cells to check biological
interpretability. A donor-bulk extension collapses the same cells into one real bulk-like
library per donor, as an independent check that results aren't an artifact of pseudobulk
construction. Configuration: `workflow/config/pseudobulk.yaml`.

```bash
# End to end
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile biology_pseudobulk panels_pseudobulk

# Donor-bulk extension
snakemake --cores 8 --resources donor_bulk_io=1 --use-conda \
  --snakefile workflow/Snakefile donor_bulk_report donor_bulk_figure2
```

## GTEx

Fits CLAMP and comparison methods on GTEx v8 bulk tissue expression, then checks whether
the latent variables recover known tissue/subtissue structure and specific biology (e.g.
liver cell-type composition via xCell). Configuration: `workflow/config/gtex.yaml`.

```bash
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile biology_gtex
```

## ARCHS4

CLAMP fit across the full ARCHS4 human RNA-seq compendium (~605k samples). Its largest
fits need ~500 GB RAM and run on Slurm, not a workstation.

The coverage analysis retrains CLAMPfull with a pinned GO:BP prior and measures pathway
recovery (Reactome, canonical non-Reactome, CellMarker) as more studies enter training,
compared against GTEx and recount2. Subsampling is at the level of GEO series, not
individual samples, and one model is fit per coverage level. The three full-data GO:BP
models are published seed-free under `output/98_models/bp_models/`.

Configuration: `workflow/config/archs4.yaml` and `workflow/config/recount2.yaml`.

### Cluster setup

The heavy stages do not fit on a workstation. Run them on a cluster through a **site
profile**: a copy of the committed generic Slurm profile with your account and partition
filled in. `.gitignore` excludes `workflow/profiles/site-*/`, so site-specific settings
never get committed.

```bash
mkdir -p workflow/profiles/site-mycluster
cp workflow/profiles/slurm/config.yaml workflow/profiles/site-mycluster/config.yaml
```

Then add your scheduler details under `default-resources` in that copy:

```yaml
default-resources:
  - slurm_account="my_account"      # sacctmgr show assoc user=$USER format=account
  - slurm_partition="my_partition"  # sinfo -s
  - mem_mb=8000
  - runtime=5760
  - tasks=1
```

Run with `--profile workflow/profiles/site-mycluster`. The profile already carries the
per-rule memory and thread requests, so nothing else needs changing.

**Memory floors.** These are requests, not suggestions: the fits die without them.

| Stage | RAM | Notes |
| --- | --- | --- |
| `preprocess_archs4` | 200 GB | 30 cores; streams the 45 GB HDF5 |
| `svd_archs4` | 200 GB | 30 cores |
| `clampbase_archs4`, `clampfull_archs4` | 500 GB | multi-day |
| coverage fits | 32 GB (1%) to 320 GB (100%) | escalates to 800 GB on OOM retry |
| `preprocess_recount2`, `svd_recount2`, `clampbase_recount2` | 200 GB | |
| ORA | 64 GB | |
| aggregation, `bp_models`, report notebook | 16 GB | runs fine on a laptop |

### What is reproducible from scratch

Everything, given the raw inputs and enough cluster time. In practice the model fits are
adopted rather than recomputed, because they take days:

```bash
# Adopt model outputs that already exist on disk, so the expensive rules never fire.
# --touch fails loudly if anything is missing, so this doubles as a completeness check.
snakemake --touch archs4_precomputed
snakemake --touch recount2_precomputed

# Check what is actually present, without running or refitting anything. This
# also checks that each coverage cell's saved sample labels align with the
# filtered ARCHS4 universe used to construct its FBM.
snakemake --profile workflow/profiles/local coverage_validate
```

`coverage_validate` writes a read-only provenance and artifact report to
`output/03_model_biology/02_archs4/00_coverage/audit/` listing every coverage cell, model
and ORA directory, and what is missing from each. It never repairs or rebuilds anything.
Cells marked `raw_indexed_misaligned` were produced by the old, uncommitted
subsampling procedure and must be regenerated before they can support a fully
reproducible coverage curve.

### Running it

```bash
# Full-compendium model fits (Slurm only, multi-day)
snakemake --profile workflow/profiles/site-mycluster archs4_models

# GO:BP coverage sweep: one model per coverage level, 1% to 100% (Slurm, multi-day)
snakemake --profile workflow/profiles/site-mycluster coverage_bp_models

# Coverage report, once fits/ORA outputs exist (local, cheap)
snakemake --profile workflow/profiles/local archs4_coverage

# Publish the three GO:BP models (ARCHS4, GTEx, recount2) to output/98_models/bp_models/.
# Hardlinks the existing fits, so this costs no extra disk.
snakemake --profile workflow/profiles/local bp_models
```

## Running a single notebook

Each notebook-backed rule runs one specific notebook and writes the executed copy
back to it. Target it by rule name like any other rule:

```bash
snakemake --cores 1 --use-conda --snakefile workflow/Snakefile liver_disentangle_xcell_rf_true_labels_gtex
```

Snakemake skips a rule whose outputs are already newer than its inputs, so editing
a notebook's cells alone won't trigger a re-run. Force one with `-f`/`--forcerun`
(add `-R`/`--forceall` to also force everything downstream of it):

```bash
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile -f <target>
```

## Citation

## License

This project is licensed under the [CC-BY 4.0 License](http://creativecommons.org/licenses/by/4.0/).

## Acknowledgments

Supported by the National Human Genome Research Institute,
The Eunice Kennedy Shriver National Institute of Child Health and Human Development,
the National Science Foundation, and the National Eye Institute.
