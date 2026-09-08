# PhenoPLIER GLS runs — summary & reproduction

LV–trait association (`phenoplier shortcut gls`, the GLS module–trait regression from
pivlab/phenoplier-cli) run against the CLAMP models this repo produces. This file is
the index of **what was actually run** and **how to reproduce it on pico**.

- Scripts live beside this file (`scripts/phenoplier/{final_models,saturation_k1728,coverage,clampbase}/`).
- The rule module `workflow/rules/phenoplier.smk` is shipped **un-wired** — see its
  header and `README.md` (§ Status). Nothing here is on the Snakefile DAG yet.
- The GLS per-model summaries produced here are the same directories Marc's
  `workflow/rules/archs4_traits.smk` (`A4_TRAIT_CFG` coverage / saturation / finals
  dirs) consumes to build its trait-recovery reports. **This is the upstream that
  produces them.**

## What was run — 86 runs across 3 CLAMP model families

All runs use `--trait-filter biomedical` (2,366 of 4,049 phenotypes) so the
Benjamini–Hochberg test count is fixed across every comparison (coverage/saturation
sweeps and the finals). Executed August–September 2026 on phenoplier-cli v0.5.2.

| Family | Finals (archs4/gtex/recount2) | Saturation k=1728 | Coverage (rs×seed) | Total |
|---|---|---|---|---|
| **CLAMPfull_bp** | 3 | 17¹ | 21² | 41 |
| **CLAMPbase** | 3 | 17¹ | 21² | 41 |
| **CLAMPfull_canonical** | 3 (archs4 full, gtex, recount2) | — | — | 3 + 1³ = 4 |
| | | | | **86** |

¹ 17 not 18 — `rs1/seed1` has no model on pico.
² 7 fractions × 3 seeds; `rs100/seed1` is byte-identical to the finals archs4 → 20 new + 1 reused.
³ CLAMPfull_canonical archs4 was additionally run as an **89-LV CRISPR-Cas subset** (a
  colleague's urgent request) — same model as the full archs4 run, restricted to the 89
  significant LVs in `significant_lvs_canonical.csv`. Deliverable = GLS summary only.

## Where the results are on pico (verified)

Workspace root: `/pividori_lab/phenoplier_workspace`. HDF5 study stores (one composite
`.rds` + GLS-results file per model, `phenoplier store build`) + normalized GLS summaries:

| Family / study | stores | summaries | store dir | summary dir |
|---|---:|---:|---|---|
| CLAMPfull_bp finals | 3 | — | `clamp_final_stores/` | in `projects/final_<ds>_CLAMPfull_bp/results/gls/phenoplier/` |
| CLAMPfull_bp saturation | 17 | — | `clamp_saturation_stores/` | in each `projects/sat_*/…` |
| CLAMPfull_bp coverage | 21 | 21 | `clamp_coverage/stores/` | `clamp_coverage/summaries/` |
| CLAMPbase finals | 3 | 3 | `clampbase_finals/stores/` | `clampbase_finals/summaries/` |
| CLAMPbase saturation | 17 | 17 | `clampbase_saturation/stores/` | `clampbase_saturation/summaries/` |
| CLAMPbase coverage | 21 | 21 | `clampbase_coverage/stores/` | `clampbase_coverage/summaries/` |
| CLAMPfull_canonical | 3 | 3 + 1 | `clampcanonical/stores/` | `clampcanonical/summaries/` (+ `clampcanonical/archs4_prelim89/`) |

**85 stores + the 89-LV subset summary** (summary-only, no store). Canonical stores:
`canon_archs4_full.h5`, `canon_gtex_CLAMPfull_canonical.h5`, `canon_recount2_CLAMPfull_canonical.h5`.

## Reproduce on pico

One directory per study; each has its own README with the exact `sbatch` lines.
Assume the run happens entirely on pico (`defq`), phenoplier-cli-neo env active.

| Study | Directory |
|---|---|
| CLAMPfull_bp finals | [`final_models/`](final_models/README.md) — `00_register` → `01_run_final_gls` → `02_build_stores` |
| CLAMPfull_bp saturation | [`saturation_k1728/`](saturation_k1728/README.md) — 3 ordered arrays, `%6`-throttled |
| CLAMPfull_bp coverage | [`coverage/`](coverage/README.md) — `pico_coverage.sbatch` + `pico_coverage_tiny.sbatch` |
| CLAMPbase (all three) | [`clampbase/`](clampbase/README.md) — `_run_clampbase_model.sh` via `finals/sat/cov_clampbase.sbatch` |

**CLAMPfull_canonical** was run on the workstation; the recipe is reproducible on pico
by the same pattern. The archs4 89-LV subset is driven manually because phenoplier-cli
has **no LV-subset knob** (step6 always builds `1~N`; step7 enumerates every model LV):
run step1–5 (`snakemake --until step5_filter_by_distance`, LV-independent), then
`pipeline gene-correlation step6-generate-lv-matrices -l "<bare numbers>"` and
`pipeline trait-association-batch --lv-list "LV… LV…"` on just the 89 LVs, then
`pipeline summarize`. step6 is per-LV idempotent, so the subsequent full run reuses
step1–5 + the 89 matrices. **The 89-LV subset was verified bit-identical to the same
89 LVs in the full 1,728-LV run (max |Δp-value| = 0)** — proof the manual subset path
is correct.

## Verification gate (every model, before a store is trusted)

Checker: [`verify_summary.py`](verify_summary.py) (also enforced inline by the CLAMPbase
runner). A model passes only if its combined summary:

- exists and has **2,366 unique phenotypes** (a partial step-7 yields a smaller
  but valid-*looking* summary — `summarize` concatenates whatever exists);
- has **0 zero-p-values** (the SE-collapse artifact that pivlab/phenoplier-cli#85 fixes);
- has **NaN only where `lv_degenerate=True`** (an expected, flagged near-singular LV).

Stores are built atomically (to a `.tmp`, verified to open with h5py, then renamed).
All 86 runs pass. Two CLAMPfull_bp models carry exactly one flagged degenerate LV each
(cov_rs50_seed2 LV161, sat_rs5_k1728_seed2 LV866) — expected, not a failure.

## Pipeline hardening (applied to every run)

- **Serial step 6** (`n_workers_step6=1`, in-process) — the forked `ProcessPoolExecutor`
  otherwise deadlocks (workers parked in `futex`, 0 % CPU).
- **Pinned step 7** — un-pinned, phenoplier over-detects a node's cores and spawns ~29
  threads/job (load ~380, ~10× slowdown). SLURM cgroup: 6×2; no-cgroup hosts: 12×1
  single-threaded (`OMP_NUM_THREADS=1` etc.).
- **Register once, then resume** — `shortcut gls --unlock` registers + writes config
  without computing; requeues do not re-init, preserving step1–6 sentinels.
- **Store build gated** on a *complete* summary + built atomically.

## Provenance — as-run vs. this reproduction subset

For throughput the 86 models were actually produced in parallel across **pico + a
workstation (`local`) + `lab` + `alpine`** (pico was busy; off-pico hosts have no cgroup,
hence the single-threaded BLAS strategy). Results are independent of which host ran a
given model. This directory ships only the **pico** reproduction path (normalized
recipes — not a byte-for-byte record of every submitted job).

The exact multi-machine launchers (`run_local_worklist.sh`, `run_one_coverage.sh`,
`alpine_coverage.sbatch`, `run_local_clampbase.sh`, `alpine_clampbase.sbatch`, the
canonical helpers) and the live cross-machine run log (`RUN_PROGRESS.md`) are preserved
verbatim on branch **`phenoplier-asrun-launchers`**, commit
**`81de511c6ee22d5fbaff1ef107ac613a1ce5bfb9`**.
