"""Build nbs/03_model_biology/00_phenoplier/02_coverage.ipynb.

Tracks the ARCHS4 coverage-study GLS run and inspects the recovery curve: how
many LV-trait associations (and distinct traits) are recovered as the coverage
fraction rs grows, per seed. Reads per-model `gls-summary` files collected into
one directory (results were produced across local + alpine, so they are gathered
into `COV_SUMMARY_DIR`, one file per model named `cov_rs<f>_seed<s>.tsv.gz`).

Rebuild:
    python nbs/03_model_biology/00_phenoplier/_build_coverage_nb.py
Execute (needs the collected summaries):
    jupyter nbconvert --to notebook --execute --inplace \\
        nbs/03_model_biology/00_phenoplier/02_coverage.ipynb
"""

import json
from pathlib import Path

NB_PATH = Path(__file__).resolve().parent / "02_coverage.ipynb"
_counter = iter(range(1000))


def md(text):
    return {"cell_type": "markdown", "id": f"c{next(_counter):02d}",
            "metadata": {}, "source": text.strip().splitlines(keepends=True)}


def code(text):
    return {"cell_type": "code", "id": f"c{next(_counter):02d}",
            "execution_count": None, "metadata": {}, "outputs": [],
            "source": text.strip().splitlines(keepends=True)}


cells = [
    md("""
# ARCHS4 coverage study — LV–trait association recovery vs coverage

`phenoplier shortcut gls` was run for the ARCHS4 coverage models across coverage
fractions `rs ∈ {1,5,10,25,50,75,100}` and seeds `{1,2,3}` (21 models; rs100/seed1
reused from the final archs4). This notebook traces how many LV–trait associations
are recovered as more of the compendium is used.

> **Biomedical subset only.** Every run used `--trait-filter biomedical`, so all
> counts are over the **2,366 of 4,049** biomedical traits — consistent across
> fractions, not comparable to unfiltered runs. Launcher + multi-machine notes:
> `scripts/phenoplier/coverage/`.

Per-model `gls-summary` files are collected into `COV_SUMMARY_DIR` (default below),
one per model named `cov_rs<f>_seed<s>.tsv.gz`.
"""),
    code('''
%matplotlib inline
import os, re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

COV = Path(os.environ.get("COV_SUMMARY_DIR", "/media/data/clamp_coverage/summaries"))
pat = re.compile(r"cov_rs(\\d+)_seed(\\d+)")

rows = []
for f in sorted(COV.glob("cov_rs*_seed*.tsv.gz")):
    m = pat.search(f.name)
    rs, seed = int(m.group(1)), int(m.group(2))
    d = pd.read_csv(f, sep="\\t", low_memory=False)
    sig = d["fdr"] < 0.05
    deg = d["lv_degenerate"].fillna(False).astype(bool) if "lv_degenerate" in d else pd.Series(False, index=d.index)
    rows.append({
        "rs": rs, "seed": seed,
        "n_phenotypes": d["phenotype"].nunique(), "n_LVs": d["lv"].nunique(),
        "n_sig_assoc": int(sig.sum()),
        "n_traits_recovered": int(d.loc[sig, "phenotype"].nunique()),
        "n_LVs_with_sig": int(d.loc[sig, "lv"].nunique()),
        "zero_pvalues": int((d["pvalue"] == 0).sum()),
        "degenerate_LVs": int(d.loc[deg, "lv"].nunique()), "NaN_fdr": int(d["fdr"].isna().sum()),
    })
res = pd.DataFrame(rows).sort_values(["rs", "seed"]).reset_index(drop=True)
print(f"{len(res)} coverage models loaded")
res
'''),
    md("## Coverage curves — recovery vs coverage fraction"),
    code('''
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for ax, col, title in [
    (axes[0], "n_sig_assoc", "significant LV–trait associations (FDR<0.05)"),
    (axes[1], "n_traits_recovered", "distinct traits recovered (>=1 sig LV)"),
]:
    for seed, g in res.groupby("seed"):
        g = g.sort_values("rs")
        ax.plot(g["rs"], g[col], marker="o", label=f"seed {seed}")
    mean = res.groupby("rs")[col].mean()
    ax.plot(mean.index, mean.values, color="black", lw=2, ls="--", label="mean")
    ax.set_xlabel("coverage fraction rs (%)"); ax.set_ylabel(col)
    ax.set_title(title); ax.legend(fontsize=8)
fig.suptitle("ARCHS4 coverage (biomedical traits)")
fig.tight_layout(); plt.show()
'''),
    md("## Per-model counts (sanity: 2,366 phenotypes, 0 zero-pvalues, 0 degenerate)"),
    code('''
piv = res.pivot_table(index="rs", columns="seed", values="n_sig_assoc")
display(piv)
bad = res[(res.n_phenotypes != 2366) | (res.zero_pvalues != 0) | (res.NaN_fdr != 0)]
print("rows failing sanity (expect none):", len(bad))
res
'''),
    md("""
## Reproduce
See `scripts/phenoplier/coverage/README.md` (local worklist 2-wide with core
pinning + alpine `%6` backfill, `--trait-filter biomedical`).
"""),
]

nb = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
    },
    "nbformat": 4, "nbformat_minor": 5,
}
NB_PATH.write_text(json.dumps(nb, indent=1) + "\n")
print(f"wrote {NB_PATH}")
