"""Build nbs/03_model_biology/00_phenoplier/01_saturation_k1728.ipynb.

Tracks the ARCHS4 saturation (k=1728) GLS run and inspects the recovery curve:
how many LV-trait associations (and distinct traits) are recovered as the
coverage fraction rs increases, per seed. Reads each model's
`gls-summary-phenomexcan.tsv.gz` from the phenoplier workspace
(one project per rs x seed). See scripts/phenoplier/saturation_k1728/.

Rebuild:
    python nbs/03_model_biology/00_phenoplier/_build_saturation_k1728_nb.py
Execute (needs the results; run where the workspace is, e.g. pico):
    jupyter nbconvert --to notebook --execute --inplace \\
        nbs/03_model_biology/00_phenoplier/01_saturation_k1728.ipynb
"""

import json
from pathlib import Path

NB_PATH = Path(__file__).resolve().parent / "01_saturation_k1728.ipynb"
_counter = iter(range(1000))


def md(text):
    return {"cell_type": "markdown", "id": f"s{next(_counter):02d}",
            "metadata": {}, "source": text.strip().splitlines(keepends=True)}


def code(text):
    return {"cell_type": "code", "id": f"s{next(_counter):02d}",
            "execution_count": None, "metadata": {}, "outputs": [],
            "source": text.strip().splitlines(keepends=True)}


cells = [
    md("""
# ARCHS4 saturation (k=1728) — LV–trait association recovery vs coverage

`phenoplier shortcut gls` was run for the ARCHS4 saturation models at fixed
**k = 1728** LVs across coverage fractions `rs ∈ {1,5,10,25,50,75}` and seeds
`{1,2,3}` — **17 models** (`rs1/seed1` is missing). This notebook traces how many
LV–trait associations are recovered as more of the compendium is used.

> **Biomedical subset only.** Every run used `--trait-filter biomedical`, so the
> 1,683 non-biomedical phenotypes are dropped before computing and all counts
> below are over the **2,366 of 4,049** biomedical traits — consistent across
> fractions, but not comparable to unfiltered runs. Launcher:
> `scripts/phenoplier/saturation_k1728/`.

Point `WS` at the workspace holding `projects/sat_rs<f>_k1728_seed<s>_CLAMPfull_bp/`
(defaults to the pico path; override with `CLAMP_PHENOPLIER_WS`).
"""),
    code('''
%matplotlib inline
import os, re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

WS = Path(os.environ.get("CLAMP_PHENOPLIER_WS", "/pividori_lab/phenoplier_workspace"))
projects = sorted((WS / "projects").glob("sat_rs*_k1728_seed*_CLAMPfull_bp"))
pat = re.compile(r"sat_rs(\\d+)_k1728_seed(\\d+)_CLAMPfull_bp")

rows = []
for p in projects:
    m = pat.search(p.name)
    rs, seed = int(m.group(1)), int(m.group(2))
    f = p / "results" / "gls" / "phenoplier" / "gls-summary-phenomexcan.tsv.gz"
    if not f.exists():
        print(f"  (missing summary) rs{rs} seed{seed}")
        continue
    d = pd.read_csv(f, sep="\\t", low_memory=False)
    sig = d["fdr"] < 0.05
    deg = d["lv_degenerate"].fillna(False).astype(bool) if "lv_degenerate" in d else pd.Series(False, index=d.index)
    rows.append({
        "rs": rs, "seed": seed,
        "n_phenotypes": d["phenotype"].nunique(), "n_LVs": d["lv"].nunique(),
        "n_sig_assoc": int(sig.sum()),
        "n_traits_recovered": int(d.loc[sig, "phenotype"].nunique()),
        "n_LVs_with_sig": int(d.loc[sig, "lv"].nunique()),
        "min_pvalue": d["pvalue"].min(), "zero_pvalues": int((d["pvalue"] == 0).sum()),
        "degenerate_LVs": int(d.loc[deg, "lv"].nunique()), "NaN_fdr": int(d["fdr"].isna().sum()),
    })
res = pd.DataFrame(rows).sort_values(["rs", "seed"]).reset_index(drop=True)
print(f"{len(res)} models loaded")
res
'''),
    md("## Saturation curves — recovery vs coverage fraction"),
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
fig.suptitle("ARCHS4 saturation, k=1728 (biomedical traits)")
fig.tight_layout(); plt.show()
'''),
    md("""
## Per-model counts

`zero_pvalues` and `degenerate_LVs` should be 0 (SE-collapse guard,
phenoplier-cli#85); `n_phenotypes` should be 2,366 everywhere (the biomedical
subset).
"""),
    code('''
piv = res.pivot_table(index="rs", columns="seed", values="n_sig_assoc")
print("significant associations by rs x seed:")
display(piv)
# sanity
bad = res[(res.n_phenotypes != 2366) | (res.zero_pvalues != 0) | (res.NaN_fdr != 0)]
print("rows failing sanity (expect none):", len(bad))
res[["rs","seed","n_phenotypes","n_LVs","n_sig_assoc","n_traits_recovered",
     "zero_pvalues","degenerate_LVs"]]
'''),
    md("""
## Reproduce

See `scripts/phenoplier/saturation_k1728/README.md` (register once → GLS array
`%6` → store array `%6`, `--trait-filter biomedical`, ≤80% of pico).
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
