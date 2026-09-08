# Final-models GLS run (PR #28, item 0)

Runs `phenoplier shortcut gls` (LV–trait association) for the three CLAMP
**final production models** — one `CLAMPfull_bp` per dataset — and composites
each result with its CLAMP model into a one-file HDF5 study store.

| Dataset  | Model (`<final_models_root>/…`)  | genes × LVs |
|----------|----------------------------------|-------------|
| ARCHS4   | `archs4/CLAMPfull_bp.rds` (~8 GB) | 18,423 × 1,728 |
| GTEx     | `gtex/CLAMPfull_bp.rds`           | 21,613 × 578   |
| Recount2 | `recount2/CLAMPfull_bp.rds`       | 6,000 × 724    |

`final_models_root` is `workflow/config/phenoplier.yaml: final_models_root`
(default: Marc's coverage-bp output `.../output/98_final_models`). The three
models are also wired into `workflow/rules/phenoplier.smk` as the `final/`
group (targets `phenoplier_final`, `phenoplier_final_stores`).

## How it was run

Three ordered sbatch steps. Everything is env-overridable; the defaults are the
**pico** layout, where the production run happened. Partition / account / QOS /
time-limit are **not** baked into the scripts — pass them as `sbatch` flags (or
rely on the cluster default), so the same scripts work on either cluster.

```bash
cd scripts/phenoplier/final_models

# pico (defq is the default partition, no account, no wall-time limit):
REG=$(sbatch --parsable 00_register_models.sbatch)
for ds in archs4 gtex recount2; do
  sbatch --dependency=afterok:$REG 01_run_final_gls.sbatch $ds
done
# after the GLS jobs finish:
sbatch 02_build_stores.sbatch

# Alpine (override paths + site policy):
#   export PHENOPLIER_CONDA_SH=/pl/active/pivlab/projects/hzhang/miniforge3/etc/profile.d/conda.sh
#   export PHENOPLIER_HOME=/pl/active/pivlab/projects/hzhang/phenoplier_workspace
#   export CLAMP_FINAL_MODELS=/pl/active/pivlab/projects/hzhang/clamp_final_models   # staged copy
#   REG=$(sbatch --parsable -p acpu --account=amc-general --qos=cpu-normal --time=03:00:00 00_register_models.sbatch)
#   sbatch --dependency=afterok:$REG -p acpu --account=amc-general --qos=cpu-long --time=2-00:00:00 01_run_final_gls.sbatch archs4
#   ... (gtex/recount2 can use --qos=cpu-normal --time=20:00:00)
```

Results land in the workspace under
`projects/final_<ds>_CLAMPfull_bp/results/gls/phenoplier/` (combined
`gls-summary-phenomexcan.tsv.gz` + per-phenotype TSVs); the stores land in
`$CLAMP_STORE_DIR` (default `<workspace>/clamp_final_stores/`) as
`final_<ds>_CLAMPfull_bp.h5`.

### Production run (pico, 2026-08-21/22, phenoplier-cli main @ 3394d44, `--trait-filter biomedical`)

All three tested 2,366 of 4,049 phenotypes (1,683 excluded), 0 zero-p-values,
0 degenerate LVs, 0 NaN FDR:

| Dataset  | GLS wall-time | summary rows | FDR<0.05 |
|----------|---------------|--------------|----------|
| ARCHS4   | ~20h          | 4,088,448    | 6,140    |
| GTEx     | ~15h          | 1,367,548    | 2,337    |
| Recount2 | ~2.7h         | 1,712,984    | 2,975    |

## Design notes / phenoplier-cli caveats hit during this run

- **Register once, then `--model`.** `shortcut gls --input <rds>` re-extracts
  and registers the model, writing the shared `models/registry.toml`
  non-atomically; three in parallel can race. So step 1 registers all three
  sequentially and step 2 uses `shortcut gls --model <key>` (read-only on the
  registry).
- **Local executor, not `--cluster`.** `shortcut gls` has no
  `--cluster/--partition/--qos` passthrough to `workflow gls init`, so its slurm
  path can't pick up a site profile. Running the whole pipeline locally inside
  one allocation avoids that and matches the validated path. (The
  `phenoplier_gls` rule still passes `--cluster` for slurm targets; that path
  needs the upstream passthrough before it works on a real cluster.)
- **Step-6/7 pool sizing.** phenoplier sizes its process pools from
  `sched_getaffinity`, so an `sbatch -c N` cgroup both caps and pins them
  correctly. An un-cgrouped run (e.g. bare `nohup`) auto-detects the whole node
  and can exhaust it — the origin of the step-7 `BrokenProcessPool` once seen on
  pico's 168-core node.
- **rpy2 + R required.** Reading a CLAMP `.rds` (in `shortcut gls` and in
  `store build --clamp-rds`) goes through rpy2 → R. The env must be created from
  phenoplier-cli's `environment.yml` (which pulls `rpy2`/`r-base`).
- **`store build` needs `anndata`** (+ zarr; the store's element encodings use
  `anndata.io`). pico is air-gapped (no PyPI), so `anndata` was installed from
  wheels transferred from a networked host: `pip install --no-index --no-deps
  <wheels>`, transferring ONLY the missing packages (+ `typing_extensions>=4.16`,
  which anndata 0.13 needs) so the env's pinned numpy/pandas are left untouched.
- **Trait filter shrinks the FDR test count.** `--trait-filter biomedical` tests
  2,366 of 4,049 phenotypes, so these numbers are not comparable to unfiltered
  runs. Keep the filter fixed across a model set.
