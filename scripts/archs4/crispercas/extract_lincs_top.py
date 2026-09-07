#!/usr/bin/env python3
# Extract top-ranked LINCS L1000 perturbagens for each significant LV.
import argparse

import pandas as pd


DESCRIPTION = "Extract top-ranked LINCS L1000 perturbagens for each significant LV."


def main():
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("--lincs", required=True)
    parser.add_argument("--lv-sig", required=True)
    parser.add_argument("--top-pct", type=float, default=0.01)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    lincs = pd.read_pickle(args.lincs)
    lv_sig = pd.read_csv(args.lv_sig)["lv"].tolist()

    missing = [lv for lv in lv_sig if lv not in lincs.index]
    if missing:
        raise SystemExit(f"LVs missing from LINCS projection: {missing}")

    n_top = max(1, int(len(lincs.columns) * args.top_pct))
    rows = []
    for lv in lv_sig:
        ranked = lincs.loc[lv].sort_values(ascending=False)
        top = ranked.head(n_top)
        rows.append(pd.DataFrame({
            "lv": lv,
            "perturbagen": top.index,
            "score": top.values,
            "rank": range(1, len(top) + 1),
        }))

    out = pd.concat(rows, ignore_index=True)
    out.to_csv(args.out, index=False)
    print(f"Wrote {args.out}: {len(lv_sig)} LVs x top {n_top} of {len(lincs.columns)} perturbagens")


if __name__ == "__main__":
    main()
