# CLAMP analyses

## 🔧 Dependencies

CLAMP uses **three Conda environments** to separate core package development, large-scale analyses, GPU-accelerated workflows, and the Snakemake orchestrator itself:

| Environment | File | Purpose |
|--------------|------|----------|
| **`clamp-analyses`** | `envs/clamp-analyses.yaml` (`.lock`) | Default environment for CPU-based modeling, priors, projections, and vignettes. |
| **`gpu-kmeans`** | `envs/gpu-kmeans.yaml` (`.lock`) | Optional environment for GPU-accelerated clustering and benchmarking. |
| **`snakemake`** | `envs/snakemake.yaml` | Runs the pipeline (`--use-conda` activates the two envs above per rule). Kept separate so the orchestrator's Snakemake version doesn't drift with either compute env's dependency set. |

This separation avoids dependency conflicts between R (Bioconductor) and GPU libraries (`rapids`, `cupy`, `cuml`), and between the pipeline orchestrator and the code it runs.

### 🛠️ Install dependencies

```bash
./setup.sh
```

Installs Conda (Miniconda) if it isn't already on `PATH`, then creates the three
environments above from their lock/yaml files and installs the external
[`chikinalab/CLAMP`](https://github.com/chikinalab/CLAMP) R package into `clamp-analyses`.
Safe to re-run: existing environments are not recreated, and CLAMP is reinstalled only
when its recorded revision is missing or differs from the required pin (remove an
environment with `conda env remove -n <name>` first to recreate it from scratch).
Currently supports Linux x86_64 only, matching the platform-pinned `.lock` files.

> [!WARNING]
> These analyses currently require CLAMP commit
> [`818e13ba55d66840e0710c3f1ac15f6d97e1dd8b`](https://github.com/chikinalab/CLAMP/commit/818e13ba55d66840e0710c3f1ac15f6d97e1dd8b).
> Do not install CLAMP from its latest branch or update it independently: newer revisions
> are not yet guaranteed to be compatible with this workflow. A follow-up PR will update
> the analyses for the latest compatible CLAMP release and remove this temporary pin.

Verify the result with `nbs/00_setup/00_check_setup.ipynb`: a read-only diagnostic
notebook that checks Conda, all three environments, Snakemake, and the data files each
rule expects (see `data/README.md` for what's auto-downloaded vs. manual).

<details>
<summary>What <code>setup.sh</code> does, step by step (for manual setup)</summary>

1. Install a Conda distribution (Miniconda or Mambaforge, either works; `setup.sh`
   installs Miniconda only if neither is already present).
2. (Optional) Verify GPU drivers for RAPIDS/cuML workflows before creating `gpu-kmeans`:

   ```bash
   nvidia-smi
   ```

   Check the RAPIDS compatibility matrix for the correct CUDA version.
3. Create the environments:

   ```bash
   conda create --name clamp-analyses --file envs/clamp-analyses.lock
   conda run -n clamp-analyses python -m pip install -r envs/clamp-analyses.pip.lock
   conda create --name gpu-kmeans --file envs/gpu-kmeans.lock
   conda env create -n snakemake -f envs/snakemake.yaml
   ```
4. Install the CLAMP R package into `clamp-analyses`:

   ```bash
   conda activate clamp-analyses
   Rscript scripts/install_clamp.R

   # Read-only verification of the required revision
   Rscript scripts/install_clamp.R --check
   ```

</details>

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

### Pipeline stages

Runs in this order; each stage consumes the previous stage's outputs:

1. __Pathway prior__ (`pathway_prior`) - downloads/verifies the GO-BP gene set reference used by CLAMP, MOFA-FLEX, and the disentangle report.
2. __Build pseudobulk data__ (`pseudobulk_data`) - per dataset: builds a cohort manifest (`cohort_manifest`), then aggregates single-cell counts into pseudobulk expression profiles + ground-truth cell-type fractions (`build_pseudobulk`).
3. __Preprocess__ (`preprocess_pseudobulk`, runs automatically per dataset ahead of each model rule) - CPM normalization, low-expression/variance filtering, and rank (`k`) estimation.
4. __Fit models__ (`full_models_pseudobulk`) - fits all 10 methods above, for every dataset.
5. __QC report__ (`model_building_qc_pseudobulk`) - sanity-checks pseudobulk inputs and model outputs across every dataset/method combination.
6. __Grouped cross-validation__ (`grouped_cv_models_pseudobulk` → `grouped_cv_analysis_pseudobulk`) - 5-fold, sample-grouped CV of CLAMPfull (held out at the donor/sample level) to estimate out-of-fold cell-type-fraction prediction accuracy.
7. __Single-cell projection__ (`single_cell_projections`) - projects each individual cell onto the CLAMPfull latent variables learned from pseudobulk.
8. __Biology reports__ (`biology_pseudobulk`) - six notebooks: `benchmark_pseudobulk` (method comparison), `holdout_report_pseudobulk` (grouped-CV results), `disentangle_pseudobulk` (LV ↔ cell-type mapping + pathway enrichment), `single_cell_recovery_pseudobulk` (single-cell projection recovery), `hard_cell_types_pseudobulk` (cell types poorly captured by any LV), and `hard_pair_loadings_pseudobulk` (gene-loading separation for difficult cell-type pairs).
9. __Computational timing__ (`computational_timing_analysis`) - a separate analysis that reuses the six preprocessed pseudobulk matrices, fits each method with seeds 123, 456, and 789, and aggregates wall time across datasets. Its models never overwrite the production models.
10. __Panels__ (`panels_pseudobulk`) - Figure 2 and Supplementary Figures 1–3, built from the biology and timing report outputs.
11. __Donor-bulk extension__ (`donor_bulk_report`) - sums the same retained raw
    single-cell counts into one library per donor, fits six full CLAMPfull models plus
    donor-grouped CV and exact-cell mean-CPM controls, projects the original cells, and
    builds independent expression UMAPs. Aggregate targets are `donor_bulk_data`,
    `donor_bulk_models`, `donor_bulk_controls`, `donor_bulk_single_cell_projections`,
    `donor_bulk_umaps`, `donor_bulk_qc`, `donor_bulk_biology`, and
    `donor_bulk_figure2`.

See "▶️ Running Snakemake" below for how to run these.

### Outputs

- Model-building artifacts: `output/01_model_building/00_pseudobulk/<dataset>/{pseudobulk,preprocessing,models,grouped_cv,single_cell_projection}/`
- QC: `output/01_model_building/00_pseudobulk/qc/`, executed notebook (with plots) at `nbs/01_model_building/00_pseudobulk/00_model_building_qc.executed.ipynb`
- Computational timing: `output/02_model_performance/00_pseudobulk/00_computational_timing/` (isolated models, per-fit logs and timing records, pooled CSVs, and the executed report notebook with its embedded plot)
- Biology reports: `output/03_model_biology/00_pseudobulk/{00_benchmark,01_holdout80,02_disentangle,03_b_matrix_singlecell,04_hard_cell_types}/`, executed notebooks (with plots) at `nbs/03_model_biology/00_pseudobulk/<name>.executed.ipynb`
- Donor-bulk models and QC: `output/01_model_building/00_pseudobulk/donor_bulk/`;
  recovery tables: `output/03_model_biology/00_pseudobulk/06_donor_bulk_recovery/`.
- Panels: `output/99_panels/`

### Configuration

All pipeline parameters (preprocessing cutoffs, CLAMP/CoGAPS/MOFA-FLEX/grouped-CV
hyperparameters and seeds, projection chunk sizes, dataset ingestion metadata, and
pathway reference files) live in `workflow/config/pseudobulk.yaml`. Adding a new dataset
means adding an entry under `datasets:` there - no rule changes needed.

## 🧬 GTEx bulk tissue pipeline

This Snakemake pipeline fits CLAMP and comparison methods on bulk RNA-seq expression
across GTEx v8 tissues, then evaluates whether the resulting latent variables recover
known tissue/subtissue structure (clustering, tissue-predictive LV importance,
subtissue inference) and specific biology (e.g. liver cell-type composition via xCell).

### Tissues

All GTEx v8 sampled tissues go into model building (one whole-body panel, not
per-tissue models). Subtissue-recovery evaluation (`subtissue_lr_eval_gtex`) is
restricted to a 9-tissue subset with multiple annotated subtissues: Adipose Tissue,
Blood, Blood Vessel, Brain, Colon, Esophagus, Heart, Kidney, Skin
(`workflow/config/gtex.yaml` → `subtissue_inference.tissues`).

### Methods benchmarked

| Method | Output | Approach |
|---|---|---|
| CLAMPbase | `CLAMPbase/B.csv` | CLAMP without pathway priors |
| CLAMPfull | `CLAMPfull/B.csv` | CLAMP with GO-BP pathway priors |
| PLIER | `PLIER/B.csv` | Pathway-level latent variable regression |
| PCA | `PCA/gtex_pca_B.pkl` | Principal component analysis |
| NMF | `NMF/gtex_nmf_B.pkl` | Non-negative matrix factorization |
| ICA | `ICA/gtex_ica_B.pkl` | Independent component analysis |
| Flashier | `flashier/gtex_B.csv` | Empirical Bayes matrix factorization |
| MOFA-FLEX | `MOFA_FLEX_PRIOR/B_matrix.csv` | Multi-omics factor analysis, with priors |
| GSSig (GSS) | `GSS/gtex_B.csv` | GenomicSuperSignature |
| CoGAPS | `CoGAPS/gtex_B.csv` | Bayesian NMF (CoGAPS), excluded by default |

Paths are relative to `output/01_model_building/01_gtex/`.

### Pipeline stages

1. **Build models** (`full_models_gtex`) - downloads/filters the GTEx v8 bulk TPM matrix
   (`download_gtex_raw` → `clamp_gtex`), then fits CLAMP (base + pathway-prior) and
   comparison methods (PLIER, PCA/NMF/ICA, Flashier, MOFA-FLEX, GSSig; CoGAPS is excluded
   by default - see the note in `gtex.smk`).
2. **QC report** (`model_building_qc_gtex`).
3. **Clustering** (`kmeans_clustering_gtex` → `kmeans_clustering_report_gtex`) - GPU
   k-means ensemble across methods and gene subsampling fractions (`gpu-kmeans` env).
4. **LV importance** (`lv_importance_rf_true_labels_gtex` → report) - random forest +
   SHAP, tissue labels as ground truth.
5. **Biology reports** (`biology_gtex`) - LV↔tissue concordance (`02_b_matrix`), global
   alignment across methods (`03_global_alignment`), liver cell-type disentangling via
   xCell (`05_liver_disentangle_xcell`), and out-of-fold subtissue recovery
   (`04_subtissues`).

See "▶️ Running Snakemake" below for how to run these.

### Outputs

- Model-building artifacts: `output/01_model_building/01_gtex/`
- Biology reports: `output/03_model_biology/01_gtex/`, executed notebooks (with plots) at
  `nbs/03_model_biology/01_gtex/<name>.executed.ipynb`
- Figure panels combining GTEx and pseudobulk results: `output/99_panels/` (same
  `panels_pseudobulk` target as above)

### Configuration

All pipeline parameters (preprocessing cutoffs, per-method hyperparameters and seeds,
clustering settings, subtissue-inference settings) live in `workflow/config/gtex.yaml`.

## ▶️ Running Snakemake

Applies to both pipelines above.

### Prerequisites

- `conda activate snakemake`, plus `clamp-analyses` created (see Dependencies above);
  pass `--use-conda` so Snakemake activates `clamp-analyses`/`gpu-kmeans` per rule.
- Run all commands from the repo root, so relative paths in the config resolve correctly.

### Basic pattern

```bash
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile <target>
```

Preview the DAG without running anything by adding `-n`:

```bash
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile -n <target>
```

### Main targets

```bash
# Pseudobulk, end to end
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile biology_pseudobulk panels_pseudobulk

# Separate 180-fit timing analysis (6 datasets x 10 methods x 3 seeds)
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile computational_timing_analysis

# Donor-level bulk-like RNA-seq extension and its Figure 2 row
snakemake --cores 8 --resources donor_bulk_io=1 --use-conda \
  --snakefile workflow/Snakefile donor_bulk_report donor_bulk_figure2

# GTEx, end to end
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile biology_gtex
```

`donor_bulk_io=1` serializes the large raw-matrix aggregation and projection scans.
Keep this resource limit when running donor-bulk targets on a workstation.

Or target any individual rule named in the "Pipeline stages" lists above, e.g.
`full_models_pseudobulk`, `kmeans_clustering_gtex`; each dataset/method/rule runs
independently where possible, so this parallelizes well with a higher `--cores`.

### Running a single notebook

Each notebook-backed rule (`notebook:` directive) runs one specific notebook and writes
the executed copy to its `log.notebook` path. Target it by rule name like any other rule:

```bash
snakemake --cores 1 --use-conda --snakefile workflow/Snakefile liver_disentangle_xcell_rf_true_labels_gtex
```

This runs `nbs/03_model_biology/01_gtex/05_liver_disentangle_xcell.ipynb`, writing its CSV
outputs under `output/03_model_biology/01_gtex/04_liver_disentangle_xcell_rf_true_labels/`
and the executed notebook (with plots) to
`nbs/03_model_biology/01_gtex/05_liver_disentangle_xcell.executed.ipynb`.

To develop/debug interactively instead of just running it, use `--edit-notebook` with one
of that rule's output files; Snakemake infers the rule from it, injects its inputs/params,
and opens a live Jupyter session, saving your edits back into the source notebook on exit:

```bash
snakemake --cores 1 --use-conda --snakefile workflow/Snakefile --edit-notebook \
  output/03_model_biology/01_gtex/04_liver_disentangle_xcell_rf_true_labels/notebook.complete
```

### Forcing a re-run

Snakemake skips a rule if its declared outputs already exist and are newer than its
inputs. Editing a notebook's cells is *not* by itself enough to trigger a re-run, since
the notebook file isn't a tracked input for most rules, only the upstream data files are.
To force a rule to run again regardless, add `-f`/`--forcerun`; add `-R`/`--forceall`
instead to also force everything downstream of it:

```bash
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile -f <target>
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile -R <target>
```

## Citation

## License

This project is licensed under the [CC-BY 4.0 License](http://creativecommons.org/licenses/by/4.0/).

## Acknowledgments

Supported by the National Human Genome Research Institute,  
The Eunice Kennedy Shriver National Institute of Child Health and Human Development,  
the National Science Foundation, and the National Eye Institute.
