# CLAMP analyses

## 🔧 Dependencies

CLAMP uses **three Conda environments** to separate core package development, large-scale analyses, and GPU-accelerated workflows:

| Environment | File | Purpose |
|--------------|------|----------|
| **`clamp-analyses.yaml`** | `envs/clamp-analyses.yaml` | Default environment for CPU-based modeling, priors, projections, and vignettes. |
| **`gpu-kmeans.yaml`** | `envs/gpu-kmeans.yaml` | Optional environment for GPU-accelerated clustering and benchmarking. |

This separation avoids dependency conflicts between R (Bioconductor) and GPU libraries (`rapids`, `cupy`, `cuml`).

### 🛠️ Install dependencies

Recommended steps to install the system-level and Conda tooling required to create the CLAMP environments.

1. Install a Conda distribution

- Install Miniconda or Mambaforge for your platform (Mambaforge is recommended for faster environment solves).

2. (Optional) Verify GPU drivers for RAPIDS/cuML workflows

- Ensure a compatible NVIDIA driver / CUDA version is installed before creating the GPU environment:

```bash
nvidia-smi
```

- Check RAPIDS compatibility matrix for the correct CUDA version (match driver/CUDA with RAPIDS/cuML requirements).

3. Create environments using conda

```bash
conda create --name clamp-analyses --file envs/clamp-analyses.lock
conda activate clamp-analyses

# Clone CLAMP the repo into REPO_PATH (adjust path as needed)
export REPO_PATH=~/path/to/CLAMP
mkdir -p "$(dirname "$REPO_PATH")"
git clone https://github.com/chikinalab/CLAMP.git "$REPO_PATH"

# Install and check CLAMP using devtools
Rscript -e "devtools::install_local('$REPO_PATH', force=TRUE, dependencies=FALSE)"
```

```bash
conda create --name gpu-kmeans --file envs/gpu-kmeans.lock
conda activate gpu-kmeans

# Clone CLAMP the repo into REPO_PATH (adjust path as needed)
export REPO_PATH=~/path/to/CLAMP
mkdir -p "$(dirname "$REPO_PATH")"
git clone https://github.com/chikinalab/CLAMP.git "$REPO_PATH"

# Install and check CLAMP using devtools
Rscript -e "devtools::install_local('$REPO_PATH', force=TRUE, dependencies=FALSE)"
```

## 📘 Notebook Headers

Each notebook explicitly states which environment to use in the first Markdown cell.

## 🧬 Pseudobulk single-cell benchmark pipeline

This Snakemake pipeline benchmarks CLAMP against other latent-variable/matrix-decomposition
methods on **pseudobulk** expression profiles built from single-cell/single-nucleus RNA-seq
cohorts, then projects the resulting CLAMP latent variables (LVs) back onto individual cells
to evaluate biological interpretability.

### Datasets

| Dataset | Tissue |
|---|---|
| `Brain_Mathys2023` | Brain |
| `Brain_Xiong2023` | Brain |
| `Heart_Datar2026` | Heart |
| `PBMC_1k1k` | PBMC |
| `PBMC_Perez2022` | PBMC |
| `Lung_Sikkema2023` | Lung (HLCA) |

Dataset ingestion metadata (raw file paths, sample/cell-type columns, cell-type label
mappings) lives in `workflow/config/pseudobulk.yaml` under `datasets:`.

### Methods benchmarked

| Method | Output | Approach |
|---|---|---|
| CLAMPbase | `models/CLAMPbase/B.csv` | CLAMP without pathway priors |
| CLAMPfull | `models/CLAMPfull/B.csv` | CLAMP with GO-BP pathway priors |
| PLIER | `models/PLIER/B.csv` | Pathway-level latent variable regression |
| PCA | `models/PCA/B.csv` | Principal component analysis |
| NMF | `models/NMF/B.csv` | Non-negative matrix factorization |
| ICA | `models/ICA/B.csv` | Independent component analysis |
| Flashier | `models/flashier/B.csv` | Empirical Bayes matrix factorization |
| MOFA-FLEX | `models/MOFA_FLEX_PRIOR/B.csv` | Multi-omics factor analysis, with priors |
| GSSig (GSS) | `models/GSS/B.csv` | GenomicSuperSignature |
| CoGAPS | `models/CoGAPS/B.csv` | Bayesian NMF (CoGAPS) |

### Prerequisites

- The `clamp-analyses` Conda environment (see Dependencies above) - pass `--use-conda` so
   Snakemake activates it automatically per rule.
- Run all commands from the repo root, so relative paths in the config resolve correctly.

### Pipeline stages

Runs in this order; each stage consumes the previous stage's outputs:

1. __Pathway prior__ (`pathway_prior`) - downloads/verifies the GO-BP gene set reference used by CLAMP, MOFA-FLEX, and the disentangle report.
2. __Build pseudobulk data__ (`pseudobulk_data`) - per dataset: builds a cohort manifest (`cohort_manifest`), then aggregates single-cell counts into pseudobulk expression profiles + ground-truth cell-type fractions (`build_pseudobulk`).
3. __Preprocess__ (`preprocess_pseudobulk`, runs automatically per dataset ahead of each model rule) - CPM normalization, low-expression/variance filtering, and rank (`k`) estimation.
4. __Fit models__ (`full_models_pseudobulk`) - fits all 10 methods above, for every dataset.
5. __QC report__ (`model_building_qc_pseudobulk`) - sanity-checks pseudobulk inputs and model outputs across every dataset/method combination.
6. __Grouped cross-validation__ (`grouped_cv_models_pseudobulk` → `grouped_cv_analysis_pseudobulk`) - 5-fold, sample-grouped CV of CLAMPfull (held out at the donor/sample level) to estimate out-of-fold cell-type-fraction prediction accuracy.
7. __Single-cell projection__ (`single_cell_projections`) - projects each individual cell onto the CLAMPfull latent variables learned from pseudobulk.
8. __Biology reports__ (`biology_pseudobulk`) - five notebooks: `benchmark_pseudobulk` (method comparison), `holdout_report_pseudobulk` (grouped-CV results), `disentangle_pseudobulk` (LV ↔ cell-type mapping + pathway enrichment), `single_cell_recovery_pseudobulk` (single-cell projection recovery), `hard_cell_types_pseudobulk` (cell types poorly captured by any LV).
9. __Panels__ (`panels_pseudobulk`) - figure 2 and supplementary figure 1, built from the biology report outputs.

### Running it

Run everything end to end:

```bash
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile --configfile workflow/config/pseudobulk.yaml biology_pseudobulk panels_pseudobulk
```

Or work through it stage by stage:

```bash
# Preview the DAG without running anything
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile -n biology_pseudobulk

# 1. Build pseudobulk expression matrices for all datasets
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile pseudobulk_data

# 2. Fit all 10 methods for all datasets (each dataset/method pair runs
#    independently, so this parallelizes well with a higher --cores value
#    or a cluster profile)
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile full_models_pseudobulk

# QC report (Python) - sanity-checks the outputs of steps 1-2
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile model_building_qc_pseudobulk

# 3. Grouped 5-fold cross-validation of CLAMPfull (sample-level holdout)
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile grouped_cv_models_pseudobulk
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile grouped_cv_analysis_pseudobulk

# 4. Project individual cells onto the CLAMPfull latent variables
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile single_cell_projections

# 5. Biology reports, run individually
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile benchmark_pseudobulk
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile holdout_report_pseudobulk
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile disentangle_pseudobulk
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile single_cell_recovery_pseudobulk
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile hard_cell_types_pseudobulk

# ...or all 5 biology reports in one go, via the umbrella rule at the
# bottom of pseudobulk.smk (it depends on all five .complete outputs)
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile biology_pseudobulk

# 6. Build the figure panels from the biology report outputs
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile panels_pseudobulk
```

Snakemake skips a rule if its declared outputs already exist and are newer than its
inputs (`Nothing to be done (all requested files are present and up to date).`). This
means editing a notebook's cells is *not* by itself enough to trigger a re-run - the
notebook file isn't a tracked input for most rules, only the upstream data files are.
To force a rule to run again regardless, add `-f`/`--forcerun`:

```bash
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile -f holdout_report_pseudobulk
```

Or force that rule and everything downstream of it with `-R`/`--forceall` on the target:

```bash
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile -R holdout_report_pseudobulk
```

### Outputs

- Model-building artifacts: `output/01_model_building/00_pseudobulk/<dataset>/{pseudobulk,preprocessing,models,grouped_cv,single_cell_projection}/`
- QC: `output/01_model_building/00_pseudobulk/qc/`, executed notebook (with plots) at `nbs/01_model_building/00_pseudobulk/00_model_building_qc.executed.ipynb`
- Biology reports: `output/03_model_biology/00_pseudobulk/{00_benchmark,01_holdout80,02_disentangle,03_b_matrix_singlecell,04_hard_cell_types}/`, executed notebooks (with plots) at `nbs/03_model_biology/00_pseudobulk/<name>.executed.ipynb`
- Panels: `output/99_panels/`

### Configuration

All pipeline parameters (preprocessing cutoffs, CLAMP/CoGAPS/MOFA-FLEX/grouped-CV
hyperparameters and seeds, projection chunk sizes, dataset ingestion metadata, and
pathway reference files) live in `workflow/config/pseudobulk.yaml`. Adding a new dataset
means adding an entry under `datasets:` there - no rule changes needed.

## Citation

## License

This project is licensed under the [CC-BY 4.0 License](http://creativecommons.org/licenses/by/4.0/).

## Acknowledgments

Supported by the National Human Genome Research Institute,  
The Eunice Kennedy Shriver National Institute of Child Health and Human Development,  
the National Science Foundation, and the National Eye Institute.
