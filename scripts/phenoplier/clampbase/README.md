# CLAMPbase — GLS run (pico)

Runs `phenoplier shortcut gls` for the **CLAMPbase** model family — the same
finals / saturation / coverage design as the `CLAMPfull_bp` runs elsewhere in
`scripts/phenoplier/`, but on the CLAMPbase models (`model-namespace clampbase`).
**41 models:** finals 3 (archs4/gtex/recount2) + saturation-k1728 17 + coverage 21.
Every run uses `--trait-filter biomedical` (2,366 of 4,049 traits).

## The hardened per-model runner

All three sbatch scripts call `_run_clampbase_model.sh <name> <rds>`, which runs
one model end-to-end (register → GLS → store) with the fixes this pipeline needs:

- **Serial step 6** (`n_workers_step6=1`, in-process) — the forked
  `ProcessPoolExecutor` used to deadlock (workers parked in `futex`, 0 % CPU).
- **Pinned step 7** (`n_workers_step7 × blas_threads_step7`; default 6×2 for a
  SLURM cgroup) — an un-pinned run over-detects the node's cores and spawns ~29
  threads/job, thrashing the box ~10×. No-cgroup hosts export
  `CLAMPBASE_NW_STEP7=12 CLAMPBASE_BLAS_STEP7=1 OMP_NUM_THREADS=1`.
- **Register once, then resume.** `shortcut gls --unlock` registers the model +
  writes the config *without* computing; on requeue it does **not** re-init, so
  the expensive step1–6 sentinels survive and the run resumes.
- **Two gates before a store is trusted:** (1) the combined summary must exist and
  be *complete* — `verify_summary.py` checks 2,366 unique phenotypes, no
  zero-p-values, NaN only where `lv_degenerate=True`; `summarize` will happily
  concatenate a *partial* step-7, so this catches an interrupted run. (2) the
  store is built to a `.tmp`, verified to open, then atomically renamed — a
  failed build never leaves a partial `.h5` that a later `-s` check would reuse.

## Reproduce (pico)

```bash
cd scripts/phenoplier/clampbase
# finals (tiered memory):
sbatch --mem=88G --array=0   finals_clampbase.sbatch     # archs4
sbatch --mem=48G --array=1-2 finals_clampbase.sbatch     # gtex, recount2
# saturation k1728 (17):
sbatch --array=0-16%14 sat_clampbase.sbatch
# coverage (21, tiered by fraction — see the script header for the tier→--mem map):
sbatch --mem=48G --array=0-8%14   cov_clampbase.sbatch   # rs1, rs5, rs10
sbatch --mem=56G --array=9-14%14  cov_clampbase.sbatch   # rs25, rs50
sbatch --mem=64G --array=15-17%14 cov_clampbase.sbatch   # rs75
sbatch --mem=88G --array=18-20%14 cov_clampbase.sbatch   # rs100
```

Each sbatch sets `CLAMPBASE_CENTRAL` to a per-study collection under the workspace
(`clampbase_{finals,saturation,coverage}/{stores,summaries}/`). `CLAMPBASE_NS`
overrides the model namespace (default `clampbase`).

> **As-run note.** These 41 were produced across pico + `local` + `lab` in
> parallel; the off-pico worklist launchers (`run_local_clampbase.sh`,
> `alpine_clampbase.sbatch`) are archived on branch `phenoplier-asrun-launchers`
> (`81de511c6ee22d5fbaff1ef107ac613a1ce5bfb9`). See `../RUN_SUMMARY.md`.
