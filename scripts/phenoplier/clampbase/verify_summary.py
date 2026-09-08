#!/usr/bin/env python
"""Completeness/validity gate for a GLS combined summary.

Exit 0 only if the summary is COMPLETE and clean:
  - exactly EXPECT_PHENOS (default 2366) unique phenotypes,
  - no zero p-values,
  - any NaN p-values are confined to rows flagged lv_degenerate=True.

Summary existence alone is NOT enough: `summarize` concatenates whatever
per-phenotype files exist, so an interrupted step 7 yields a smaller-but-valid-
looking summary. This checks the phenotype count.

Usage: verify_summary.py <summary.tsv.gz> [expected_phenos]
"""
import sys
import pandas as pd

f = sys.argv[1]
expect = int(sys.argv[2]) if len(sys.argv) > 2 else 2366

try:
    d = pd.read_csv(f, sep="\t", low_memory=False)
except Exception as e:  # noqa: BLE001
    print(f"FAIL read: {e}")
    sys.exit(1)

nph = d["phenotype"].nunique()
zero_p = int((d["pvalue"] == 0).sum())
nan_p = d["pvalue"].isna()
nan_total = int(nan_p.sum())
if "lv_degenerate" in d.columns:
    nan_bad = int((nan_p & ~d["lv_degenerate"].fillna(False)).sum())
else:
    nan_bad = nan_total

ok = (nph == expect) and (zero_p == 0) and (nan_bad == 0)
status = "OK" if ok else "FAIL"
print(
    f"{status} {f}: phenos={nph}/{expect} zero_p={zero_p} "
    f"nan_p={nan_total} nan_nondegenerate={nan_bad} lvs={d['lv'].nunique()}"
)
sys.exit(0 if ok else 1)
