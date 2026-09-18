#!/usr/bin/env python
"""Completeness/validity gate for a GLS combined summary.

Exit 0 only if the summary is COMPLETE and clean:
  - exactly EXPECT_PHENOS (default 2366) unique phenotypes,
  - a full phenotype x LV grid: no duplicate (phenotype, lv) pair, and every
    phenotype tested against the same LV set,
  - every non-NaN p-value finite and in (0, 1] (a zero p-value is the SE
    collapse artifact pivlab/phenoplier-cli#85 fixes),
  - NaN p-values confined to rows flagged lv_degenerate=True.

Summary existence alone is NOT enough: `summarize` concatenates whatever
per-phenotype files exist, so an interrupted step 7 yields a smaller-but-valid-
looking summary. The phenotype count and the grid check catch that.

Usage: verify_summary.py <summary.tsv.gz> [expected_phenos]
"""
import sys

import numpy as np
import pandas as pd

f = sys.argv[1]
expect = int(sys.argv[2]) if len(sys.argv) > 2 else 2366

try:
    d = pd.read_csv(f, sep="\t", low_memory=False)
except Exception as e:  # noqa: BLE001
    print(f"FAIL read: {e}")
    sys.exit(1)

problems = []
for col in ("phenotype", "lv", "pvalue"):
    if col not in d.columns:
        problems.append(f"missing column {col}")
if problems:
    print(f"FAIL {f}: " + "; ".join(problems))
    sys.exit(1)

nph = d["phenotype"].nunique()
nlv = d["lv"].nunique()
if nph != expect:
    problems.append(f"phenotypes {nph} != {expect}")

dup = int(d.duplicated(["phenotype", "lv"]).sum())
if dup:
    problems.append(f"{dup} duplicate (phenotype, lv) rows")
per_pheno = d.groupby("phenotype")["lv"].nunique()
short = int((per_pheno != nlv).sum())
if short:
    problems.append(f"{short} phenotypes missing LVs (grid {nph} x {nlv} = {nph * nlv}, rows {len(d)})")

p = pd.to_numeric(d["pvalue"], errors="coerce")
nan_p = p.isna()
nan_total = int(nan_p.sum())
if "lv_degenerate" in d.columns:
    degenerate = d["lv_degenerate"].fillna(False).astype(bool)
else:
    degenerate = pd.Series(False, index=d.index)
nan_bad = int((nan_p & ~degenerate).sum())
if nan_bad:
    problems.append(f"{nan_bad} NaN p-values on non-degenerate LVs")

zero_p = int((p == 0).sum())
if zero_p:
    problems.append(f"{zero_p} zero p-values")
out_of_range = int((~nan_p & (~np.isfinite(p) | (p < 0) | (p > 1))).sum())
if out_of_range:
    problems.append(f"{out_of_range} p-values outside (0, 1] or non-finite")

status = "OK" if not problems else "FAIL"
print(
    f"{status} {f}: phenos={nph}/{expect} lvs={nlv} rows={len(d)} zero_p={zero_p} "
    f"nan_p={nan_total} nan_nondegenerate={nan_bad}"
    + ("" if not problems else " :: " + "; ".join(problems))
)
sys.exit(0 if not problems else 1)
