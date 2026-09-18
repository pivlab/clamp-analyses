"""Build nbs/03_model_biology/00_phenoplier/00_final_models_inspection.ipynb.

Inspects and visualizes the LV–trait association (GLS) results for the three
CLAMP **final production models** (archs4 / gtex / recount2 CLAMPfull_bp), run
via `phenoplier shortcut gls` (see scripts/phenoplier/final_models/). It reads
the per-model `gls-summary-phenomexcan.tsv.gz` from the phenoplier workspace and
reports counts, p-value distributions, significant associations and the traits
each model recovers.

Rebuild with:
    python nbs/03_model_biology/00_phenoplier/_build_final_models_inspection.py

Execute (needs the results + pandas/matplotlib) e.g. on pico:
    jupyter nbconvert --to notebook --execute --inplace \\
        nbs/03_model_biology/00_phenoplier/00_final_models_inspection.ipynb
"""

import json
from pathlib import Path

NB_PATH = Path(__file__).resolve().parent / "00_final_models_inspection.ipynb"


_counter = iter(range(1000))


def md(text):
    return {
        "cell_type": "markdown",
        "id": f"cell-{next(_counter):02d}",
        "metadata": {},
        "source": text.strip().splitlines(keepends=True),
    }


def code(text):
    return {
        "cell_type": "code",
        "id": f"cell-{next(_counter):02d}",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": text.strip().splitlines(keepends=True),
    }


cells = [
    md("""
# Final CLAMP models — LV–trait association (GLS) results

Inspection of `phenoplier shortcut gls` output for the three **final production
models**, one `CLAMPfull_bp` per dataset, run against the `phenomexcan_rapid_gwas`
cohort. The exact launcher is committed under `scripts/phenoplier/final_models/`.

> **The GLS step ran on the biomedical trait subset only.** Every run used
> `--trait-filter biomedical`, which **drops the 1,683 non-biomedical phenotypes**
> (occupational/job codes, medication codes, 24-hour diet-recall, household/admin
> and environmental-exposure fields) *before* computing — so results below cover
> the **2,366 of 4,049** biomedical traits, and the excluded traits are never
> tested. This also shrinks the Benjamini–Hochberg test count, so these numbers
> are not comparable to an unfiltered run (keep the filter fixed across a model
> set). See `phenoplier-cli` `docs/development/trait-filtering-standard.md`.

| Dataset | genes × LVs |
|---|---|
| ARCHS4 | 18,423 × 1,728 |
| GTEx | 21,613 × 578 |
| Recount2 | 6,000 × 724 |

Results live in the phenoplier workspace under
`projects/final_<ds>_CLAMPfull_bp/results/gls/phenoplier/`. Point `WS` below at
that workspace (defaults to the pico path; override with `CLAMP_PHENOPLIER_WS`).
"""),
    code('''
%matplotlib inline
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

WS = Path(os.environ.get("CLAMP_PHENOPLIER_WS", "/pividori_lab/phenoplier_workspace"))
DATASETS = ["archs4", "gtex", "recount2"]

def summary_path(ds):
    return (WS / "projects" / f"final_{ds}_CLAMPfull_bp"
            / "results" / "gls" / "phenoplier" / "gls-summary-phenomexcan.tsv.gz")

gls = {}
for ds in DATASETS:
    p = summary_path(ds)
    if p.exists():
        gls[ds] = pd.read_csv(p, sep="\\t", low_memory=False)
        print(f"{ds:9s} {len(gls[ds]):>10,} rows  <-  {p}")
    else:
        print(f"{ds:9s} MISSING  ({p})")
'''),
    md("""
## Counts per model

`lv_degenerate` marks fits the phenoplier-cli#85 SE-collapse guard flagged; their
`fdr` is `NaN` and never counts as significant.
"""),
    code('''
rows = []
for ds, d in gls.items():
    deg = d["lv_degenerate"].fillna(False).astype(bool) if "lv_degenerate" in d else pd.Series(False, index=d.index)
    rows.append({
        "dataset": ds,
        "phenotypes": d["phenotype"].nunique(),
        "LVs": d["lv"].nunique(),
        "assoc (rows)": len(d),
        "min pvalue": d["pvalue"].min(),
        "zero pvalues": int((d["pvalue"] == 0).sum()),
        "degenerate LVs": int(d.loc[deg, "lv"].nunique()),
        "sig FDR<0.05": int((d["fdr"] < 0.05).sum()),
        "NaN fdr": int(d["fdr"].isna().sum()),
    })
counts = pd.DataFrame(rows).set_index("dataset")
counts
'''),
    md("""
## P-value distributions

Under the null the p-values are ~uniform; a spike near 0 is real signal. A spike
of *exact* zeros would indicate the SE-collapse artifact (phenoplier-cli#85) — we
expect none here.
"""),
    code('''
fig, axes = plt.subplots(1, len(gls), figsize=(5 * len(gls), 3.4), squeeze=False)
for ax, (ds, d) in zip(axes[0], gls.items()):
    ax.hist(d["pvalue"].dropna(), bins=50, color="#4C72B0", edgecolor="white")
    ax.axhline(len(d) / 50, color="grey", ls="--", lw=1, label="uniform (null)")
    ax.set_title(f"{ds}  (n={len(d):,})")
    ax.set_xlabel("GLS p-value")
    ax.set_ylabel("count")
    ax.legend(fontsize=8)
fig.tight_layout()
plt.show()
'''),
    md("""
## Significant LV–trait associations (FDR < 0.05)

Note: `--trait-filter biomedical` shrinks the Benjamini–Hochberg test count, so
these are comparable across the three models here but **not** to unfiltered runs.
"""),
    code('''
sig = counts["sig FDR<0.05"]
fig, ax = plt.subplots(figsize=(5, 3.2))
bars = ax.bar(sig.index, sig.values, color=["#4C72B0", "#DD8452", "#55A868"])
ax.bar_label(bars, fmt="{:,.0f}")
ax.set_ylabel("significant (FDR < 0.05)")
ax.set_title("LV–trait associations recovered per model")
fig.tight_layout()
plt.show()
'''),
    md("""
## Top associations per model
"""),
    code('''
cols = [c for c in ["lv", "phenotype", "phenotype_desc", "pvalue", "fdr"] if c in next(iter(gls.values())).columns]
for ds, d in gls.items():
    print(f"=== {ds}: top 12 by FDR ===")
    display(d.sort_values("fdr")[cols].head(12).reset_index(drop=True))
'''),
    md("""
## Traits recovered

How many distinct LV modules each trait is significantly associated with, per
model — the traits with the broadest module signal. Uses phenotype descriptions
so the same trait lines up across models even though the LV sets differ.
"""),
    code('''
def sig_per_trait(d):
    s = d[d["fdr"] < 0.05]
    key = "phenotype_desc" if "phenotype_desc" in d.columns else "phenotype"
    return s.groupby(key)["lv"].nunique().rename("n_sig_LVs")

per = pd.concat({ds: sig_per_trait(d) for ds, d in gls.items()}, axis=1).fillna(0).astype(int)
per["total"] = per.sum(axis=1)
print(f"{len(per):,} traits significant in >=1 model")
per.sort_values("total", ascending=False).head(20)
'''),
    code('''
# Overlap of traits with >=1 significant module association across the models.
sig_traits = {}
for ds, d in gls.items():
    key = "phenotype_desc" if "phenotype_desc" in d.columns else "phenotype"
    sig_traits[ds] = set(d.loc[d["fdr"] < 0.05, key])
for ds, s in sig_traits.items():
    print(f"{ds:9s} {len(s):>5,} traits with a significant module")
if len(sig_traits) == 3:
    common = set.intersection(*sig_traits.values())
    print(f"\\nshared by all three: {len(common):,}")
'''),
    md("""
## Reproducing

See `scripts/phenoplier/final_models/README.md`. In short (pico):

```bash
cd scripts/phenoplier/final_models
REG=$(sbatch --parsable 00_register_models.sbatch)
for ds in archs4 gtex recount2; do
  sbatch --dependency=afterok:$REG 01_run_final_gls.sbatch $ds
done
sbatch 02_build_stores.sbatch      # one HDF5 study store per model
```
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
