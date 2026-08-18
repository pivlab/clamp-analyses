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
fits need ~500 GB RAM and run on Slurm, not a workstation. The coverage analysis retrains
with a GO:BP prior and measures pathway recovery (Reactome, canonical non-Reactome,
CellMarker) as more studies enter training, compared against GTEx and recount2.
Configuration: `workflow/config/archs4.yaml`.

```bash
# Full-compendium model fits (Slurm only, multi-day) -- see
# workflow/profiles/slurm/config.yaml for the generic profile to copy and fill in
snakemake --profile /path/to/user/slurm-profile archs4_models

# Coverage report, once fits/ORA outputs exist (local, cheap)
snakemake --profile workflow/profiles/local archs4_coverage
```

## Citation

## License

This project is licensed under the [CC-BY 4.0 License](http://creativecommons.org/licenses/by/4.0/).

## Acknowledgments

Supported by the National Human Genome Research Institute,
The Eunice Kennedy Shriver National Institute of Child Health and Human Development,
the National Science Foundation, and the National Eye Institute.
