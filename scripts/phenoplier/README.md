# phenoplier-cli integration — LV–trait association (GLS)

Runs `phenoplier shortcut gls` (GLS module–trait regression, from
pivlab/phenoplier-cli) against the CLAMP models this repo's ARCHS4 workflow fits,
and publishes one GLS summary per model where `workflow/rules/archs4_traits.smk`
picks them up for the coverage / saturation / final-model trait-recovery reports.

phenoplier-cli needs `snakemake>=9` / Python ≥3.12, which is incompatible with
this repo's own `snakemake=8` / Python 3.11 (`envs/snakemake.yaml`). So it runs in
its **own conda env** (`phenoplier-cli-neo`); each GLS run is one shell job of
this repo's DAG (`run_gls.sh`) that drives phenoplier-cli's pipeline in one
workspace project per model.

**phenoplier-cli is pinned to v0.5.2** (`workflow/config/phenoplier.yaml: version`),
the release the 86 as-run models were produced with (`RUN_SUMMARY.md`). Every GLS
and store job checks `phenoplier --version` and fails on a mismatch, and
`setup_env.sh` installs that tag. To move to a newer release, bump the config key
and `setup_env.sh` together and re-run the whole model set — results from
different releases are not comparable.

## How it is wired

`workflow/rules/phenoplier.smk` is included by `workflow/Snakefile` (config:
`workflow/config/phenoplier.yaml`). It has one rule per model family, and each
takes its model from the rule file that builds that model — no model path is
written down here:

| rule | model, from | summary, into (`archs4.yaml: traits.*`) |
|---|---|---|
| `gls_bp_coverage_full` | `a4_cov_model_dir()` CLAMPfull_bp, after `validate_bp_coverage_model` | `coverage.clampfull_bp_dir/cov_rs<f>_seed<s>.tsv.gz` |
| `gls_bp_coverage_base` | `a4_cov_cell()/CLAMPbase.rds` | `coverage.clampbase_dir/cov_rs<f>_seed<s>.tsv.gz` |
| `gls_bp_saturation_full` | `a4_sat_model_dir()` CLAMPfull_bp, after `validate_bp_saturation_model` | `saturation.clampfull_bp_dir/sat_rs<f>_k<k>_seed<s>.tsv.gz` |
| `gls_bp_saturation_base` | `a4_sat_base_rds()` | `saturation.clampbase_dir/sat_rs<f>_k<k>_seed<s>.tsv.gz` |
| `gls_bp_final_full` | gtex / recount2: `final_models.clampfull.bp.<ds>`, after `publish_bp_model` | `finals.clampfull_bp_dir/final_<ds>.tsv.gz` |
| `gls_final_base` | gtex / recount2: `coverage.clampbase_source.<ds>` (`clamp_gtex`, `clampbase_recount2`) | `finals.clampbase_dir/final_<ds>.tsv.gz` |
| `link_bp_final_full`, `link_final_base` | archs4: the finals *are* the coverage rs100/seed1 cells, so their summaries are copied, not re-run | `finals.*_dir/final_archs4.tsv.gz` |
| `gls_canonical` | `A4_CAN_FINAL_ROOT/<ds>/CLAMPfull_canonical.rds` | beside the model: `<ds>/traits/canon_<ds>.tsv.gz` |
| `store_bp_final_full`, `store_final_base` | the finals above | `finals.*_dir/../stores/final_<ds>_<model>.h5` |

The filenames are the ones `scripts/archs4/traits/aggregate_*_traits.R` parse; the
directory, not the filename, says whether a summary is CLAMPfull_bp or CLAMPbase.
`aggregate_archs4_{coverage,saturation}_traits` list these files as inputs, so

```bash
snakemake --profile workflow/profiles/local archs4_traits
```

fits any missing model, runs GLS on it and renders the trait-recovery reports in
one DAG. Workspace project names (`cov_rs5_seed2_CLAMPfull_bp`, …) match the
as-run pico names, so a workspace that already holds a finished project is reused
rather than recomputed. `run_gls.sh` fingerprints the `.rds` it registered
(`<project>/clamp_source.sha256`); if a model is refit, the next run moves the old
project aside (`<project>.stale.<timestamp>`), re-registers the model and starts
over. Projects made before this wiring carry no fingerprint and are reused as-is
(with a warning) — delete the project and `phenoplier model remove clamp/<name>`
to force a recompute.

### Targets

| target | what |
|---|---|
| `phenoplier_coverage` / `phenoplier_saturation` / `phenoplier_finals` | GLS summaries for one family (CLAMPfull_bp + CLAMPbase) |
| `phenoplier_canonical` | GLS summaries for the three canonical-prior models |
| `phenoplier_final_stores` | one HDF5 composite study store per final model (`phenoplier store build`) |
| `phenoplier_traits` | every per-model summary (coverage, saturation, finals, canonical) |
| `phenoplier_long_table` | every LV x trait row of every summary, streamed into one gzipped CSV under `phenoplier.yaml: report_root` (tens of GB) |
| `archs4_traits` | the trait-recovery reports (archs4_traits.smk); pulls the GLS runs it needs |

Dry-run first — a GLS job is hours to a day per model:

```bash
snakemake -n --snakefile workflow/Snakefile phenoplier_finals
```

### Configuration (`workflow/config/phenoplier.yaml`)

- **`trait_filter`** — `biomedical` (2,366 of 4,049 phenotypes) or `all`. Filtering
  shrinks the Benjamini–Hochberg test count, so keep it fixed across a comparison.
- **`saturation_k_values`** — ranks to run; the report is K=1728 only.
- **`workers`** — phenoplier-cli's step-6/7 pools. Serial step 6 (the forked pool
  deadlocks otherwise) and a pinned step 7 (`step7 × blas_step7` = the job's
  threads). On a host without a cgroup use many single-threaded workers.
- **`target` / `clusters`** — which of phenoplier-cli's *own* cluster profiles
  `shortcut gls` submits to. `local` (the validated path) runs each model's whole
  pipeline inside the allocation this repo's profile gives the job; override per
  run with `--config phenoplier_target=server_cu`.
- **`resources.gls`** — threads / memory / wall clock this repo's scheduler requests
  per GLS job.

## Layout

| path | what |
|---|---|
| `setup_env.sh` | one-time: build the `phenoplier-cli-neo` env (phenoplier-cli **v0.5.2** + rpy2/R) |
| `run_gls.sh` | one model → one GLS summary (register + init, pin pools, run, gate, publish); driven by the rules |
| `store_build.sh` | one model → one HDF5 study store; driven by the rules |
| `verify_summary.py` | completeness gate: expected phenotype count, full phenotype x LV grid, p-values in (0, 1], degenerate-only NaN |
| `aggregate_traits.{py,sh}` | stream per-model summaries into one long cross-model table (chunked, gzipped) |
| `final_models/`, `saturation_k1728/`, `coverage/`, `clampbase/` | the as-run pico `sbatch` recipes (see `RUN_SUMMARY.md`) |
| `RUN_SUMMARY.md` | **what was actually run** (86 models across 3 CLAMP families) + reproduce-on-pico index |

## One-time setup (per machine/cluster)

```bash
bash scripts/phenoplier/setup_env.sh            # builds phenoplier-cli-neo
conda activate phenoplier-cli-neo
phenoplier workspace init
phenoplier workspace link /path/to/phenoplier_full_data
```

The workspace defaults to `$HOME/phenoplier`; set `PHENOPLIER_HOME` to use another
(e.g. `/pividori_lab/phenoplier_workspace` on pico). `setup_env.sh` installs
phenoplier-cli **v0.5.2** (the pinned release; `phenoplier --version` must report
it). On an air-gapped cluster install the env from wheels transferred from a
networked host — see `final_models/README.md`.

## Reproduce the as-run pico jobs

The 86 models in `RUN_SUMMARY.md` were produced with the `sbatch` recipes under
`final_models/`, `saturation_k1728/`, `coverage/` and `clampbase/` before this
module was wired into the Snakefile. They remain the record of what ran; the rules
above are the reproducible path from here on and go through the same hardened
flow (`run_gls.sh` ≡ `clampbase/_run_clampbase_model.sh` minus the store step).
