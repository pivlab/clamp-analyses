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

### Cell-type evaluation universe

All pseudobulk and donor-bulk recovery analyses use the same canonical set of
47 dataset-cell-type targets. The following low-prevalence annotations are excluded:

- `Brain_Mathys2023`: `Vas`
- `Brain_Xiong2023`: `Vas`
- `Lung_Sikkema2023`: `Hematopoietic stem cells`, `Mesothelium`, and `Smooth muscle`

The exclusions are defined once under `cell_type_analysis.excluded_targets` in
`workflow/config/cell_type_analysis.yaml` and apply to full-data LV matching, grouped-CV
calibration and evaluation, projected single-cell recovery, specificity analyses,
method comparisons, and figure/supplementary tables. The corresponding cells remain
in the pseudobulk and donor-bulk expression mixtures and in the unfiltered composition
tables. Thus the libraries continue to represent all retained cells and truth-table
rows continue to sum to one; the excluded annotations are simply not scored as
recovery targets and the retained target fractions are not renormalized.

### Single-cell recovery provenance

The projected-cell purity results from the pseudobulk and donor-bulk models are
separate benchmarks and must not be pooled or substituted for one another:

| Figure | Model fitted to | Pooled purity | Mean dataset purity | Mean purity ratio | Mean lift | Perez B-cell LV |
|---|---|---:|---:|---:|---:|---:|
| Supplementary Fig. 1e | Pseudobulk mixtures | 70.3% | 76.8% | 0.773 | 19.1 | LV10 |
| Supplementary Fig. 2a | Donor-summed bulk-like libraries | 70.6% | 77.5% | 0.782 | 19.6 | LV32 |

Supplementary Figure 1 reads the pseudobulk recovery outputs under
`output/03_model_biology/00_pseudobulk/03_b_matrix_singlecell/`.
Supplementary Figure 2 reads the donor-bulk recovery outputs under
`output/03_model_biology/00_pseudobulk/06_donor_bulk_recovery/`. Both panel
notebooks assert these pipeline identities and benchmark values before drawing.

For grouped-CV results, the live canonical path is
`output/01_model_building/00_pseudobulk/grouped_cv_analysis/`. The retired
`output/03_model_biology/00_pseudobulk/01_grouped_cv/` directory may contain stale
pre-exclusion results and must not be used by analyses or figures.

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
pathway reference files) live in `workflow/config/pseudobulk.yaml`. The shared
pseudobulk/donor-bulk evaluation exclusions live separately in
`workflow/config/cell_type_analysis.yaml`, so changing the evaluated target universe
does not invalidate aggregation or model fits. Adding a new dataset means adding an
entry under `datasets:` in the pseudobulk config - no rule changes needed.

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

## 🧬 ARCHS4 compendium pipeline

CLAMP fit across the whole ARCHS4 human RNA-seq compendium: **18,423 genes x 605,614
samples**, K = 1728. This is the model the coverage, saturation, drug-disease, projection
and CRISPR analyses all project into.

Unlike the pipelines above, ARCHS4 does not fit on a workstation. Its largest CLAMP fits
ask for ~500 GB of RAM, so they run through a generic high-memory Slurm profile. The
streaming ORA can run either on Slurm or locally. See "Running ARCHS4" below.

### Input

ARCHS4 `human_gene_v2.5.h5` (raw gene counts). Counts are TPM-normalized using a
per-symbol transcript-length table pulled once from Ensembl 107
(`data/archs4/gene_lengths.rds`); TPM rather than CPM because the compendium ships raw
counts and TPM puts it on the same footing as GTEx, which enters the repo already in TPM.

### Pipeline stages

1. **Preprocess** (`preprocess_archs4`) - stream the HDF5 in blocks, collapse duplicate
   gene symbols, TPM-normalize, clean, drop probable single-cell libraries, filter on mean
   and variance, z-score. Produces the file-backed matrix everything else reads.
2. **SVD** (`svd_archs4`) - randomized SVD and the CLAMP rank estimate, per seed.
3. **CLAMPbase** (`clampbase_archs4`) - unsupervised fit, per seed.
4. **CLAMPfull** (`clampfull_archs4`) - fit with the Hallmark + Reactome + GO:CC + C8
   prior, per seed. This is the shipped model.
5. **QC** (`model_building_qc_archs4`) - matrix shape, pathway recovery per seed, and
   seed-to-seed agreement.

At 100% of the compendium there is no subsampling, so stages 2 and 3 are computed once and
shared: the SVD and CLAMPbase fit are byte-identical across seeds, and only CLAMPfull
actually varies.

### Analyses

Each lives in its own directory under `nbs/03_model_biology/02_archs4/`, with the heavy
compute in `scripts/archs4/` and the notebook reduced to plotting:

| Directory | Target | What it asks |
|---|---|---|
| `00_coverage/` | `coverage_bp` | Retrain CLAMPfull with the GO:BP prior and measure Reactome, canonical non-Reactome, and CellMarker recovery as more **studies** enter training. Includes a full-data ARCHS4/GTEx/recount2 comparison with an explicit universe per dataset. |
| `01_saturation/` | `archs4_saturation` | How does coverage depend on model rank K, and does a larger compendium need a larger K? |
| `02_drug_diseases_associations/` | - | S-PrediXcan and LINCS projections, module- vs gene-based drug-disease prediction across ARCHS4/GTEx/recount2. |
| `03_projections/` | - | GTEx and external datasets projected into the ARCHS4 model. |
| `04_crispercas/` | - | CRISPR-Cas9 screen enrichment and LV deep dives. |

Subsampling is **by study**, not by random samples: samples within a study are correlated,
so drawing them at random overstates how much new biology each increment buys. The
random-sample track is retired to `output/_deprecated/`.

Coverage and saturation are scored by ORA only. Trait- and L1000-based coverage do not
exist yet.

### Outputs

- Model-building artifacts: `output/01_model_building/02_archs4/`
  - `00_preprocess/` filtered, z-scored FBM plus gene/sample metadata
  - `01_final_model/hall_coverage_rs100_seed_{1,2,3}/` the shipped models
  - `02_coverage_study/`, `03_saturation_study/` subsampling sweeps (by study)
  - `04_reconstruction_error/`, `05_drug_disease/`
- Biology reports: `output/03_model_biology/02_archs4/`
- Executed notebooks (with plots) at `nbs/**/02_archs4/<name>.executed.ipynb`

### Configuration

`workflow/config/archs4.yaml` - preprocessing cutoffs, CLAMP hyperparameters, the prior
GMT collections, ORA settings, and the final-model location.

## ▶️ Running ARCHS4

### Adopting the existing results

The ARCHS4 models already exist (2.2 TB). Adopt them instead of recomputing:

```bash
# One-time: move legacy output trees into the layout the rules declare.
# Dry run first; nothing is deleted, dropped tracks go to output/_deprecated/.
scripts/archs4/migrate_outputs.sh
scripts/archs4/migrate_outputs.sh --apply

# Mark every existing output as up to date.
snakemake --profile workflow/profiles/local --touch archs4_precomputed
```

`--touch` fails loudly on a missing file, so it doubles as a check that the migration put
everything where the rules expect it. Afterwards `snakemake -n archs4_models` should
report "Nothing to be done".

**Starting from scratch instead?** Skip both steps; the pipeline computes everything from
the raw HDF5. Expect multi-day jobs.

### Coverage rebuild resources

The coverage-only rebuild reuses its existing subsamples, SVDs, inferred ranks, and
CLAMPbase fits. Each row below is one CLAMPfull seed; the wall-time limits are deliberately
generous and do not represent expected elapsed time.

| ARCHS4 fraction | CPUs | Initial RAM | Wall time |
|---:|---:|---:|---:|
| 1% | 16 | 32 GB | 96 h |
| 5% | 16 | 64 GB | 144 h |
| 10% | 16 | 96 GB | 240 h |
| 25% | 16 | 224 GB | 400 h |
| 50% | 16 | 400 GB | 600 h |
| 75% | 16 | 400 GB | 800 h |
| 100% | 16 | 500 GB | 900 h |

GTEx and recount2 request 8 CPUs, 32 GB, and 96 h per seed. A model/database ORA
requests 1 CPU, 8 GB, and 96 h. Failed CLAMPfull jobs have two retries; their memory
requests increase to 1.5x and 2x the initial value, capped at 800 GB. Outputs are
published from a temporary directory only after their dimensions and finite values pass
validation. Each attempt records its exit code, elapsed time, peak RSS, requested
resources, seed, prior hash, and Slurm state when available; the final validation also
captures `sacct` accounting.

On a single 900 GB node, restrict aggregate fit memory to 800 GB and allow only one
50--100% fit at a time. Streaming ORA jobs can use spare resources, but should not delay a
pending CLAMPfull job. If the Slurm node is occupied, synchronize the completed model's
`Z.csv` and manifest and run its ORA through the local profile.

Reactome, canonical non-Reactome, and CellMarker are run as separate ORAs. Each uses the
genes in that dataset's `Z` as its explicit universe; BH correction and the eligible-term
denominator remain database-specific and are combined only in the reporting tables.

#### Reading the cross-dataset panel

Per-dataset universes mean the three compendia are not scored against the same
denominator, and the difference is large:

| Dataset | ORA universe | Eligible Reactome | Eligible canonical | Eligible CellMarker |
|---|---:|---:|---:|---:|
| ARCHS4 | 18,423 genes | 1,370 | 1,721 | 289 |
| GTEx | 21,613 genes | 1,381 | 1,735 | 289 |
| recount2 | 6,000 genes | 1,148 | 1,513 | 216 |

This is the intended design: an ORA universe must be the genes the model could actually
have recovered, and recount2 enters the repo already filtered to 6,000 genes. It does mean
the cross-dataset figure's absolute pathway counts are **not** a like-for-like comparison,
and that part of recount2's lower recovery is denominator, not model quality. Compare the
recovered percentages, which are per-dataset, before comparing the counts.

At 100% there is no subsampling, so the SVD and CLAMPbase fit are shared across seeds.
The CLAMPbase reference marks at that fraction therefore carry no seed-to-seed spread by
construction; only CLAMPfull varies.

### Running on Slurm

```bash
# Run this on the Slurm login node after providing account/partition values in
# a user-owned profile derived from workflow/profiles/slurm/config.yaml.
snakemake --profile /path/to/user/slurm-profile coverage_bp_archs4_models

# Fit the two smaller compendia locally (at most two fits at once).
snakemake --profile workflow/profiles/local coverage_bp_comparator_models
```

The tracked profile contains no hostname, SSH alias, account, partition, email address,
or personal path. Inspect actual usage with `sacct`, including `State`, `Elapsed`,
`MaxRSS`, `AllocCPUS`, and `ReqMem`, then refine site-local settings if necessary.

For compute nodes without internet access, pre-stage both environments and all four GMT/
XLSX inputs, set `use-conda: false` only in the site-owned profile, and prepend the
pre-staged analysis and Snakemake `bin` directories to the worker `PATH`. Some Slurm
executors still call `conda env export` for provenance even when environment activation is
disabled, so the worker must have an offline `conda` command capable of exporting the
declared environment.

After remote model jobs begin, `scripts/coverage/route_ora.py` can run on the local host.
It discovers validated remote models, checks free CPUs/RAM and pending CLAMPfull jobs,
records one backend assignment per model/database, and invokes the same Snakemake ORA
target remotely or locally. Its SSH destination, roots, and user-owned Slurm profile are
runtime arguments and are not stored in the repository. Use `--dry-run` first; use
`--once` for a single routing pass or omit it to poll continuously.

The router also counts local fit and ORA processes launched by other Snakemake instances,
reserves headroom for their requested memory, and checks current system-available RAM.
Thus independent local dispatchers still share the 96 GB coverage budget instead of each
assuming it owns the full machine.

Both profiles set `rerun-triggers: mtime`. Without it, editing an ARCHS4 rule's shell
command or params marks its outputs stale through Snakemake's provenance tracking, which
would queue multi-day 500 GB jobs to regenerate files that are already correct. If you do
edit these rules and see them queued unexpectedly, re-run the `--touch` above.

### Deriving a site Slurm profile

`workflow/profiles/slurm/config.yaml` is generic on purpose: no account, partition,
hostname or absolute path. A site profile is that file plus the local answers, kept
outside version control (`.gitignore` covers `workflow/profiles/site-*/`):

```bash
mkdir -p workflow/profiles/site-slurm
cp workflow/profiles/slurm/config.yaml workflow/profiles/site-slurm/config.yaml
# Then edit the copy to add the two site-specific deltas:
#   default-resources:  add slurm_partition (and slurm_account if not guessable)
#   use-conda: false    only when the compute nodes have no internet and the
#                       environments are pre-staged on PATH instead
```

Drop the `set-resources` block from the copy if the site profile only drives the coverage
rules, which carry their own per-fraction memory and runtime. Keep the copy on the machine
that submits jobs; the campaign scripts refer to it by path, never by content.

### Running the campaign unattended

The coverage campaign runs for days across two machines, so its drivers must outlive the
shell that starts them. `scripts/coverage/run_campaign.sh` starts the host half under
`setsid`, is safe to re-run at any time, and starts only what is not already running:

```bash
scripts/coverage/run_campaign.sh --status    # progress, cluster jobs, nothing started
scripts/coverage/run_campaign.sh --dry-run   # what it would start
scripts/coverage/run_campaign.sh             # start or resume the host half
```

It reads site values from `.campaign.env` (gitignored; override with `--site-env`):

```bash
CLAMP_SNAKEMAKE_BIN=/path/to/env/snakemake/bin   # host environments
CLAMP_ANALYSES_BIN=/path/to/env/clamp-analyses/bin
CLAMP_REMOTE=my-cluster-alias                    # SSH alias for the Slurm login node
CLAMP_REMOTE_ROOT=/path/to/checkout/on/cluster
CLAMP_REMOTE_PROFILE=workflow/profiles/site-slurm
CLAMP_REMOTE_PATH_DIRS="--remote-path-dir /pre/staged/snakemake/bin --remote-path-dir /pre/staged/analyses/bin"
CLAMP_TOTAL_MEM_MB=96000                         # host scheduling budget
CLAMP_FIT_SLOTS=2                                # drop to 1 under memory pressure
```

The script starts four components: the comparator model fits, the comparator ORAs (which
wait for those fits to validate), the ORA router, and the aggregation plus notebook (which
waits for every ORA). Both waiters derive their expected counts from
`workflow/config/archs4.yaml`, so adding a fraction, seed, database or comparator does not
leave them waiting on a stale number.

The cluster half is separate and is not started by this script. Launch the high-memory
fits inside a `tmux` session on the login node so they survive disconnection, and check
them with `ssh <alias> 'tmux ls; squeue -u $USER'`.

The umbrella `coverage_bp` target only resolves where the ARCHS4 subsamples, SVDs and
CLAMPbase fits live, which is the cluster. On a workstation it raises
`MissingInputException` for `fit_bp_archs4_coverage`, so drive the host half through the
narrower targets the script already uses (`coverage_bp_comparator_models`,
`coverage_bp_comparator_ora`, `aggregate_bp_coverage coverage_report_archs4`). Aggregation
needs only the published ORA directories, which the router synchronizes back, not the
models themselves.

If the host drivers are killed, re-running the script resumes them. Snakemake recomputes
nothing already finished, and the router reclaims any ORA left mid-flight: on startup it
discards output directories that were never published atomically and makes their targets
eligible again.

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

# ARCHS4 models (cluster only -- see "Running ARCHS4" above before invoking this,
# or you will queue multi-day 500 GB jobs instead of adopting existing results)
snakemake --profile workflow/profiles/slurm archs4_models

# ARCHS4 model-building QC and the coverage/saturation reports (local, cheap)
snakemake --profile workflow/profiles/local \
  model_building_qc_archs4 archs4_coverage archs4_saturation
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
