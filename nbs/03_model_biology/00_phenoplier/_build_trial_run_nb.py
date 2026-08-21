"""Build nbs/03_model_biology/00_phenoplier/00_trial_run_archs4_final_seed1.ipynb.

The notebook documents the first end-to-end run of the phenoplier-cli
integration against a real production CLAMP model, and doubles as the
inspection tool for its output: run it (or re-run it) once the pipeline
finishes and every cell reports against the live project directory.

Rebuild with:
    python nbs/03_model_biology/00_phenoplier/_build_trial_run_nb.py
"""

import json
from pathlib import Path

NB_PATH = Path(__file__).resolve().parent / "00_trial_run_archs4_final_seed1.ipynb"


def md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text.strip().splitlines(keepends=True)}


def code(text):
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": text.strip().splitlines(keepends=True),
    }


cells = [
    md("""
# Trial run: ARCHS4 final model (seed 1) through phenoplier-cli

First end-to-end run of this PR's integration against a **real production
CLAMP model**, rather than the 10-LV smoke model used earlier. It exists to
answer one question: *does `phenoplier shortcut gls` work, unattended, on a
model of the size this repo actually produces?*

| | |
|---|---|
| Model | `output/01_model_building/02_archs4/01_final_model/seed1/CLAMPfull.rds` (7.5 GB, **1,728 LVs**) |
| Source checkout | `/pividori_lab/marc_projects/clamp-analyses-coverage-bp` (pico) |
| phenoplier-cli | **v0.5.1** + the fix from [phenoplier-cli#100](https://github.com/pivlab/phenoplier-cli/pull/100) |
| Cohort / panel / eQTL | `phenomexcan_rapid_gwas` / GTEX_V8 / MASHR |
| Trait filter | `biomedical` — 2,366 of 4,049 traits kept |
| Executor | local, 22 parallel jobs (pico: 336 cores, 1 TB RAM) |
| Project dir | `/pividori_lab/phenoplier_workspace/projects/archs4_final_seed1_CLAMPfull` |

The exact launcher is committed as
`scripts/phenoplier/trial_run_archs4_seed1.sh`, so this is reproducible
rather than a one-off shell invocation.
"""),
    md("""
## What the run found

Three defects surfaced before the pipeline reached real computation. Each is
worth recording, because each would have hit anyone running this integration
for the first time.

**1. The trait filter crashed on every real workspace** — *fixed upstream in
[phenoplier-cli#100](https://github.com/pivlab/phenoplier-cli/pull/100).*

```
AttributeError: 'Settings' object has no attribute 'TWAS'
RuntimeError: init failed (exit 1): ... workflow gls init ... --trait-filter biomedical
```

Outside dev/test, phenoplier's Dynaconf object starts with **no settings file
loaded**; each command loads the workspace's settings itself. `pipeline init`
never had to, because it takes every path as an argument — so the trait
filter, the first code there to read the UK Biobank data dictionary, hit an
unconfigured `Settings`. phenoplier-cli's own test suite could not catch it:
it runs with `ENV_FOR_DYNACONF=test`, where the templates are pre-loaded.

**2. `snakemake` not found** — *fixed in this PR's `run_gls.sh`.*

`phenoplier` shells out to a **bare** `snakemake`, and each Snakemake rule
shells out to a bare `phenoplier`. Invoking the CLI by absolute path is not
enough; the env's `bin` must be on `PATH`. `run_gls.sh` now exports it
explicitly rather than relying on `conda activate` having run hooks.

**3. The summary copy path was wrong** — *fixed in this PR's `run_gls.sh`.*

The script copied `results/gls/gls-summary-phenomexcan.tsv.gz`, but summaries
are written to `results/gls/**phenoplier**/` and named after the cohort
(`gls-summary-<cohort-slug>.tsv.gz`, with a `-phenomexcan` alias for the
default cohort). Every run would have completed the science and then failed
on the final copy. Confirmed against a finished project on pico.
"""),
    code('''
from pathlib import Path
import gzip
import pandas as pd

PROJECT = Path(
    "/pividori_lab/phenoplier_workspace/projects/archs4_final_seed1_CLAMPfull"
)
SENTINELS = PROJECT / "pipeline_results" / "sentinels"
SUMMARY_DIR = PROJECT / "results" / "gls" / "phenoplier"

def step_status():
    """Which pipeline steps have completed, from the on-disk sentinels."""
    if not SENTINELS.exists():
        return pd.DataFrame({"step": [], "done": []})
    done = {p.stem.replace(".done", "") for p in SENTINELS.glob("*.done")}
    rows = []
    for step in ["step1", "step2", "step3_complete", "step4", "step5", "step6", "step7"]:
        rows.append({"step": step, "done": step in done})
    n_chr = len([d for d in done if d.startswith("step3_chr")])
    rows.append({"step": "step3 (per-chromosome)", "done": f"{n_chr}/22"})
    return pd.DataFrame(rows)

step_status()
'''),
    md("""
## Configuration actually used

Read back from the generated `pipeline_config.yaml` rather than restated, so
this reflects what ran.
"""),
    code('''
import yaml

cfg = yaml.safe_load((PROJECT / "pipeline_config.yaml").read_text())
{
    k: cfg[k]
    for k in ["cohort", "reference_panel", "eqtl_model", "num_lvs",
              "lv_percentile", "trait_filter"]
    if k in cfg
} | {"n_phenotype_files": len(cfg.get("phenotype_files", []))}
'''),
    md("""
## The trait filter, on real data

`--trait-filter biomedical` is applied at `init`, so excluded traits are never
computed: **4,049 discovered → 2,366 kept, 1,683 excluded**. Every exclusion is
logged with its category and the rule that dropped it — the log is copied next
to the results by `run_gls.sh`, so the provenance travels with the numbers.
"""),
    code('''
excluded = pd.read_csv(PROJECT / "trait_filter_excluded.tsv", sep="\\t")
print(f"{len(excluded):,} traits excluded")
display(excluded["exclusion_reason"].value_counts().rename("traits").to_frame())
excluded.head(5)
'''),
    code('''
# Spot-check the borderline decisions the standard makes (see phenoplier-cli
# docs/development/trait-filtering-standard.md): treatment codes go, heritable
# behavioural phenotypes stay.
kept_codes = set(Path(f).name for f in cfg.get("phenotype_files", []))
print("medication codes excluded:",
      int(excluded["exclusion_reason"].str.startswith("treatment").sum()))
print("job codes excluded:",
      int(excluded["category"].isin(["Employment history", "Employment"]).sum()))
print("phenotype files that will actually be tested:", len(kept_codes))
'''),
    md("""
## Results

Empty until step 7 + summarize finish. `pvalue`/`fdr` may be `NaN` for
degenerate LVs — that is the SE-collapse guard from phenoplier-cli#85 doing
its job, and those fits are excluded from the BH-FDR count rather than
counted as significant.
"""),
    code('''
summary_path = SUMMARY_DIR / "gls-summary-phenomexcan.tsv.gz"
if not summary_path.exists():
    print(f"Not finished yet — no summary at {summary_path}")
    print("Steps still pending:")
    display(step_status())
else:
    gls = pd.read_csv(summary_path, sep="\\t")
    print(f"{len(gls):,} (LV, trait) rows | "
          f"{gls['lv'].nunique():,} LVs x {gls['phenotype'].nunique():,} traits")
    print(f"significant (FDR < 0.05): {(gls['fdr'] < 0.05).sum():,}")
    if "lv_degenerate" in gls.columns:
        deg = gls["lv_degenerate"].fillna(False).astype(bool)
        print(f"degenerate fits (excluded from FDR): {int(deg.sum()):,} "
              f"across {gls.loc[deg, 'lv'].nunique()} LV(s)")
    display(gls.sort_values("fdr").head(10))
'''),
    code('''
# Sanity check: the trait axis must match the filter, not the raw 4,049.
if summary_path.exists():
    n_traits = gls["phenotype"].nunique()
    print(f"traits in results: {n_traits:,} (expected 2,366 with --trait-filter biomedical)")
    assert n_traits <= 2366, "more traits than the filter should have kept"
    excluded_codes = set(excluded["phenotype"].astype(str))
    leaked = excluded_codes & set(gls["phenotype"].astype(str))
    print(f"excluded traits that leaked into results: {len(leaked)} (expected 0)")
'''),
    md("""
## Reproducing / resuming

```bash
# on pico, from this repo
bash scripts/phenoplier/trial_run_archs4_seed1.sh
```

The pipeline is idempotent: re-running skips every step whose sentinel exists
and resumes at the first missing one. If a run was interrupted, clear the
stale Snakemake lock first:

```bash
cd /pividori_lab/phenoplier_workspace/projects/archs4_final_seed1_CLAMPfull
snakemake --configfile pipeline_config.yaml --profile profiles/gls_local --unlock
```

**Timing note (from Marc's earlier smoke test, confirmed here):** steps 1–5
run against the full reference panel and do **not** scale down with LV count,
so a "small" model is not a fast test. Only steps 6–7 scale with LVs — and
with the trait filter, step 7 now does 41.6% less work.
"""),
]

nb = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

NB_PATH.parent.mkdir(parents=True, exist_ok=True)
NB_PATH.write_text(json.dumps(nb, indent=1) + "\n")
print(f"wrote {NB_PATH}")
