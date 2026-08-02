"""
GPU-accelerated one-sided hypergeometric (Fisher's exact, "greater") ORA against
MSigDB-style GMT gene sets, vectorized over all (LV, pathway) pairs at once via cupy.

Reproduces the same statistical convention as the clusterProfiler/decoupleR ORA
notebooks in this repo: top-1%-per-LV hit set, gene universe = model genes that
overlap the GMT, pathway size filter [min_size, max_size], BH adjustment within
each LV, then minimum adjusted p-value per pathway across LVs.
"""

import numpy as np
import pandas as pd
import cupy as cp
from cupyx.scipy.special import gammaln


def read_gmt(gmt_path):
    """Tab-delimited GMT parser: field[0]=term name, field[2:]=genes."""
    gene_sets = {}
    with open(gmt_path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            gene_sets[fields[0]] = fields[2:]
    return gene_sets


def build_pathway_matrix(universe_genes, library, min_size=10, max_size=50000):
    """
    Returns (P, term_names): P is an (N, n_pathways_filtered) cupy 0/1 array
    aligned to universe_genes (rows) and term_names (columns), restricted to
    pathways whose overlap with universe_genes is within [min_size, max_size].
    """
    gene_index = {g: i for i, g in enumerate(universe_genes)}
    n_genes = len(universe_genes)

    term_names = []
    cols = []
    for term, genes in library.items():
        idx = [gene_index[g] for g in genes if g in gene_index]
        if min_size <= len(idx) <= max_size:
            term_names.append(term)
            col = np.zeros(n_genes, dtype=np.float32)
            col[idx] = 1.0
            cols.append(col)

    if len(cols) == 0:
        raise ValueError("No pathways passed the min_size/max_size filter")

    P = cp.asarray(np.stack(cols, axis=1))
    return P, term_names


def build_hit_matrix(Z, universe_genes, pct=0.01):
    """
    Returns (H, n_top): H is an (n_LV, N) cupy 0/1 array, top-pct genes per LV
    column by descending loading, restricted to universe_genes (in that order).
    n_top is constant across all LVs (depends only on len(universe_genes)).
    """
    Zu = Z.loc[universe_genes]
    n_genes = len(universe_genes)
    n_top = int(np.ceil(pct * n_genes))

    n_lv = Zu.shape[1]
    H = np.zeros((n_lv, n_genes), dtype=np.float32)
    values = Zu.to_numpy()
    for l in range(n_lv):
        top_idx = np.argpartition(-values[:, l], n_top - 1)[:n_top]
        H[l, top_idx] = 1.0

    return cp.asarray(H), n_top


def hypergeom_ora(H, P, N, n_top):
    """
    Core kernel. H: (n_LV, N) hit matrix, P: (N, n_pathways) pathway matrix,
    N: universe size, n_top: hit-set size (constant across LVs).

    Returns pvals: (n_LV, n_pathways) cupy float64 array, one-sided P(X >= x).
    """
    n_pathways = P.shape[1]

    # overlap counts: single GPU matmul
    X = cp.asarray(cp.rint(H @ P), dtype=cp.int64)  # (n_LV, n_pathways)

    pw_size = cp.asarray(cp.rint(P.sum(axis=0)), dtype=cp.int64)  # (n_pathways,)

    x_range = cp.arange(n_top + 1, dtype=cp.float64)  # (n_top+1,)
    K = pw_size.astype(cp.float64)[:, None]  # (n_pathways, 1)

    const = gammaln(n_top + 1.0) + gammaln(N - n_top + 1.0) - gammaln(N + 1.0)

    logpmf = (
        gammaln(K + 1.0) - gammaln(x_range[None, :] + 1.0) - gammaln(K - x_range[None, :] + 1.0)
        + gammaln(N - K + 1.0) - gammaln(n_top - x_range[None, :] + 1.0) - gammaln(N - K - n_top + x_range[None, :] + 1.0)
        + const
    )  # (n_pathways, n_top+1)

    # mask combinatorially invalid (pathway, x) cells to -inf before exponentiating
    invalid = (x_range[None, :] > K) | (x_range[None, :] < cp.maximum(0.0, n_top - N + K))
    logpmf = cp.where(invalid, -cp.inf, logpmf)

    # reverse cumulative tail sum P(X>=x), stable via row-max shift
    row_max = logpmf.max(axis=1, keepdims=True)
    exp_shifted = cp.exp(logpmf - row_max)
    exp_shifted = cp.where(cp.isneginf(logpmf), 0.0, exp_shifted)
    reverse_cumsum = cp.cumsum(exp_shifted[:, ::-1], axis=1)[:, ::-1]
    tail_logp = cp.log(reverse_cumsum) + row_max  # (n_pathways, n_top+1)

    # gather: pvals[l, t] = exp(tail_logp[t, X[l, t]])
    col_idx = cp.arange(n_pathways)[None, :]
    pvals = cp.exp(tail_logp[col_idx, X])  # (n_LV, n_pathways)

    return cp.clip(pvals, 0.0, 1.0)


def bh_adjust_rows(pvals):
    """
    Row-wise (per-LV) Benjamini-Hochberg adjustment, matching R's
    p.adjust(method="BH") exactly (rank-based scaling + tail cummin monotonicity).
    pvals: (n_rows, m) cupy array. Returns padj of the same shape.
    """
    m = pvals.shape[1]
    order = cp.argsort(pvals, axis=1)  # ascending p-value order, per row
    sorted_p = cp.take_along_axis(pvals, order, axis=1)

    scale = m / cp.arange(1, m + 1, dtype=cp.float64)
    raw_adj = sorted_p * scale[None, :]

    # enforce non-decreasing from the tail (R's rev(cummin(rev(...))))
    # cupy has no minimum.accumulate; cummin is cheap relative to the gammaln/matmul
    # steps, so do it on CPU with numpy (which does support it) and transfer back.
    raw_adj_np = cp.asnumpy(raw_adj)
    adj_sorted_np = np.minimum.accumulate(raw_adj_np[:, ::-1], axis=1)[:, ::-1]
    adj_sorted = cp.clip(cp.asarray(adj_sorted_np), 0.0, 1.0)

    # scatter back to original column order: rank[i] = position of column i in `order`
    rank = cp.argsort(order, axis=1)
    padj = cp.take_along_axis(adj_sorted, rank, axis=1)
    return padj


def combine_min_across_lvs(padj_per_lv, term_names):
    """Column-wise min across LVs -> pandas Series indexed by term_names."""
    terms_padj = cp.asnumpy(padj_per_lv.min(axis=0))
    return pd.Series(terms_padj, index=term_names)


def run_gpu_ora_for_model(z_path, library, min_size=10, max_size=50000, pct=0.01):
    """
    End-to-end per-model driver. Returns a dict with terms_padj (pd.Series),
    n_samples, n_lvs, n_top_genes, n_total_msigdb -- same schema as the R
    run_ora_for_model() closures used by the decoupler_ora notebooks.
    """
    Z = pd.read_csv(z_path, index_col=0)

    # universe = all model genes (matches clusterProfiler::enricher's
    # `universe <- rownames(Z)` convention) -- NOT intersected with the GMT a
    # priori; pathway genes are intersected with this universe when building P.
    universe_genes = list(Z.index)
    N = len(universe_genes)

    P, term_names = build_pathway_matrix(universe_genes, library, min_size=min_size, max_size=max_size)
    n_total_msigdb = len(term_names)

    H, n_top = build_hit_matrix(Z, universe_genes, pct=pct)

    pvals = hypergeom_ora(H, P, N, n_top)
    padj_per_lv = bh_adjust_rows(pvals)
    terms_padj = combine_min_across_lvs(padj_per_lv, term_names)

    # n_samples = number of columns in B.csv's header. B.csv can be multi-GB with a
    # single very wide header line (one column per sample), so read just that line
    # directly rather than pd.read_csv(nrows=0), which still tokenizes the full
    # header through the C parser and is dramatically slower for wide headers.
    b_path = "/".join(str(z_path).split("/")[:-1]) + "/B.csv"
    try:
        with open(b_path) as f:
            header = f.readline()
        n_samples = header.rstrip("\n").count(",")
    except FileNotFoundError:
        n_samples = None

    return {
        "n_samples": n_samples,
        "n_lvs": Z.shape[1],
        "n_top_genes": n_top,
        "n_total_msigdb": n_total_msigdb,
        "terms_padj": terms_padj,
    }
