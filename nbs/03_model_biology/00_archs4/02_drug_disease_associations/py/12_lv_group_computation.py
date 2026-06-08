# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all,-execution,-papermill,-trusted,-slideshow
#     notebook_metadata_filter: -jupytext.text_representation.jupytext_version
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: clamp-analyses
#     language: python
#     name: python3
# ---

# %% [markdown]
# 💡 **Environment:** `clamp-analyses`  
#

# %% [markdown]
# # LV group computation
#
# Computes drug-disease prediction performance for:
# 1. All methods (ARCHS4, GTEx, recount2, gene-based) from aggregated predictions.
# 2. ARCHS4 LV groups (Pathway+Trait / Pathway only / Trait only / Neither) by running the full prediction pipeline on LV subsets.

# %% [markdown]
# # Setup
#

# %%
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.metrics import roc_auc_score, average_precision_score
from pyprojroot import here

sys.path.insert(0, str(here('libs')))
from drug_disease_utils import map_traits_to_doid

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

DATA_DIR   = here('data/drug_disease_associations')
OUTPUT_DIR = here('output/03_model_biology/00_archs4/02_drug_disease_associations/12_lv_group_computation')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# %% [markdown]
# # Data loading
#

# %%
gold = pd.read_pickle(DATA_DIR / 'gold_standard.pkl')
positive_rate = gold['true_class'].mean()
print(f'Gold standard: {len(gold)} pairs, positive rate={positive_rate:.3f}')

preds_avg = pd.read_pickle(
    here('output/03_model_biology/00_archs4/02_drug_disease_associations/10_prediction_performance') / 'lincs' / 'predictions' / 'predictions_results_aggregated.pkl')
preds_avg['method'] = preds_avg['method'].cat.rename_categories({
    'module_based_archs4':   'ARCHS4 LV-based',
    'module_based_gtex':     'GTEx LV-based',
    'module_based_recount2': 'recount2 LV-based',
})
print(f'Predictions: {preds_avg.shape}')
print('Methods:', preds_avg['method'].cat.categories.tolist())

ukb_efo = pd.read_csv(DATA_DIR / 'phenomexcan_traits_fullcode_to_efo.tsv',
                      sep='\t', index_col='ukb_fullcode')
ukb_efo.index = [idx.replace('-', '_', 1) for idx in ukb_efo.index]
efo_xrefs     = pd.read_csv(DATA_DIR / 'term_id_xrefs.tsv.gz', sep='\t')
do_xrefs      = pd.read_csv(DATA_DIR / 'xrefs-prop-slim.tsv',  sep='\t')
preferred_doids = set(gold['trait'])

lincs_proj = pd.read_pickle(
    here('output/03_model_biology/00_archs4/02_drug_disease_associations/01_lincs_projection_archs4') / 'lincs' / 'lincs-projection.pkl')
print(f'LINCS proj: {lincs_proj.shape}')

spr_files = sorted(
    (here('output/03_model_biology/00_archs4/02_drug_disease_associations/00_spredixcan_projection_archs4') / 'spredixcan' / 'proj').glob('*-projection-archs4.pkl'))
print(f'S-PrediXcan tissue files: {len(spr_files)}')


# %% [markdown]
# # LV group definitions
#

# %%
CLAMP_MODEL_FILE = here('output/01_model_building/04_archs4/06_bp_coverage_rshall/06_bp_coverage_hall_rs_100/hall_coverage_rs100_seed_1/CLAMPfull_hall.rds')
assert CLAMP_MODEL_FILE.exists()

readRDS = ro.r['readRDS']
clamp   = readRDS(str(CLAMP_MODEL_FILE))
with localconverter(ro.default_converter + pandas2ri.converter):
    model_sum = ro.conversion.rpy2py(clamp.rx2('summary'))
print(f'Pathway summary: {model_sum.shape}')

traits_df = pd.read_csv(
    here() / 'data' / 'archs4' / 'traits' /
    'hall_coverage_rs100_seed_1_CLAMPfull_hall' /
    'gls-summary-phenomexcan.tsv.gz', sep='\t')

all_lvs_set     = set(lincs_proj.index)
sig_path_strict = set(model_sum[(model_sum.AUC > 0.7) & (model_sum.FDR < 0.05)]['LV'].unique()) & all_lvs_set
sig_path_loose  = set(model_sum[(model_sum.AUC > 0.6) & (model_sum.FDR < 0.1) ]['LV'].unique()) & all_lvs_set
sig_trait       = set(traits_df[traits_df.fdr < 0.05]['lv'].unique()) & all_lvs_set

def make_groups(sp):
    return {
        'Pathway + Trait': sorted(sp & sig_trait),
        'Pathway only':    sorted(sp - sig_trait),
        'Trait only':      sorted(sig_trait - sp),
        'Neither':         sorted(all_lvs_set - sp - sig_trait),
    }

GROUP_ORDER = ['Pathway + Trait', 'Pathway only', 'Trait only', 'Neither']
THRESHOLDS  = {
    'strict (AUC>0.7, FDR<0.05)': make_groups(sig_path_strict),
    'loose (AUC>0.6, FDR<0.1)':   make_groups(sig_path_loose),
}
for label, groups in THRESHOLDS.items():
    print(f'{label}:')
    for g, lvs in groups.items():
        print(f'  {g}: {len(lvs)} LVs')


# %% [markdown]
# # Computation
#
# ## Bootstrap helper
#

# %%
np.random.seed(42)
N_BOOT = 300

def bootstrap_metrics(y_true, y_score, n_boot=N_BOOT):
    y_true  = np.array(y_true)
    y_score = np.array(y_score)
    n = len(y_true)
    aurocs, auprcs = [], []
    for _ in range(n_boot):
        idx = np.random.choice(n, n, replace=True)
        if len(np.unique(y_true[idx])) < 2:
            continue
        aurocs.append(roc_auc_score(y_true[idx], y_score[idx]))
        auprcs.append(average_precision_score(y_true[idx], y_score[idx]))
    if not aurocs:
        pt = roc_auc_score(y_true, y_score)
        pa = average_precision_score(y_true, y_score)
        return np.array([pt]), np.array([pa])
    return np.array(aurocs), np.array(auprcs)



# %% [markdown]
# ## Method performance
#

# %%
METHOD_ORDER = ['ARCHS4 LV-based', 'GTEx LV-based', 'recount2 LV-based', 'gene_based']
method_perf  = {}
for method in METHOD_ORDER:
    df_m = preds_avg[preds_avg['method'] == method]
    aurocs_b, auprcs_b = bootstrap_metrics(df_m['true_class'], df_m['score'])
    method_perf[method] = dict(
        auroc    = aurocs_b.mean(),
        auroc_lo = np.percentile(aurocs_b, 2.5),
        auroc_hi = np.percentile(aurocs_b, 97.5),
        auprc    = auprcs_b.mean(),
        auprc_lo = np.percentile(auprcs_b, 2.5),
        auprc_hi = np.percentile(auprcs_b, 97.5),
    )

df_r = preds_avg[preds_avg['method'] == 'gene_based'].copy()
df_r['score'] = np.random.permutation(df_r['score'].values)
aurocs_r, auprcs_r = bootstrap_metrics(df_r['true_class'], df_r['score'])
method_perf['Random'] = dict(
    auroc    = aurocs_r.mean(),
    auroc_lo = np.percentile(aurocs_r, 2.5),
    auroc_hi = np.percentile(aurocs_r, 97.5),
    auprc    = auprcs_r.mean(),
    auprc_lo = np.percentile(auprcs_r, 2.5),
    auprc_hi = np.percentile(auprcs_r, 97.5),
)
METHOD_ORDER_FULL = METHOD_ORDER + ['Random']

for m, v in method_perf.items():
    print(f'{m:25s}  AUROC={v["auroc"]:.3f}  AUPRC={v["auprc"]:.3f}')


# %% [markdown]
# ## LV group performance
#

# %%
def compute_group_perf(lv_subset, gold_df, lincs, spr_files,
                       ukb_efo, efo_xrefs, do_xrefs, preferred_doids):
    if len(lv_subset) == 0:
        return None
    lincs_sub = lincs.loc[lv_subset]
    pair_max  = {}
    for spr_file in spr_files:
        tissue_proj = pd.read_pickle(spr_file)
        tissue_sub  = tissue_proj.loc[lv_subset]
        scores      = -1.0 * lincs_sub.T.dot(tissue_sub)
        scores_doid = map_traits_to_doid(scores, preferred_doids, ukb_efo, efo_xrefs, do_xrefs)
        valid_doids = [d for d in gold_df['trait'].unique() if d in scores_doid.columns]
        valid_drugs = [d for d in gold_df['drug'].unique()  if d in scores_doid.index]
        if not valid_doids or not valid_drugs:
            continue
        sub = scores_doid.loc[valid_drugs, valid_doids]
        for drug in sub.index:
            for doid in sub.columns:
                key = (drug, doid)
                val = sub.loc[drug, doid]
                if key not in pair_max or val > pair_max[key]:
                    pair_max[key] = val
    if not pair_max:
        return None
    rows = [{'drug': r['drug'], 'trait': r['trait'],
             'score': pair_max[(r['drug'], r['trait'])], 'true_class': r['true_class']}
            for _, r in gold_df.iterrows() if (r['drug'], r['trait']) in pair_max]
    return pd.DataFrame(rows) if rows else None


grp_results = {}
for thresh_label, groups in THRESHOLDS.items():
    grp_results[thresh_label] = {}
    for grp_name in GROUP_ORDER:
        lvs = groups[grp_name]
        print(f'[{thresh_label}] {grp_name} ({len(lvs)} LVs) ... ', end='', flush=True)
        df = compute_group_perf(lvs, gold, lincs_proj, spr_files,
                                ukb_efo, efo_xrefs, do_xrefs, preferred_doids)
        if df is not None and len(np.unique(df['true_class'])) > 1:
            aur  = roc_auc_score(df['true_class'], df['score'])
            aupr = average_precision_score(df['true_class'], df['score'])
            aurocs_b, auprcs_b = bootstrap_metrics(df['true_class'], df['score'])
            grp_results[thresh_label][grp_name] = {
                'auroc':    aur,
                'auroc_lo': np.percentile(aurocs_b, 2.5),
                'auroc_hi': np.percentile(aurocs_b, 97.5),
                'auprc':    aupr,
                'auprc_lo': np.percentile(auprcs_b, 2.5),
                'auprc_hi': np.percentile(auprcs_b, 97.5),
                'n':        len(df),
                'n_pos':    int(df['true_class'].sum()),
            }
            print(f'AUROC={aur:.3f}  AUPRC={aupr:.3f}')
        else:
            print('skip')

best_thresh = max(grp_results,
    key=lambda k: np.mean([v['auroc'] for v in grp_results[k].values() if 'auroc' in v]))
print(f'\nBest threshold: {best_thresh}')


# %% [markdown]
# # Save results
#

# %%
import pickle
results = {
    'method_perf':      method_perf,
    'method_order_full': METHOD_ORDER_FULL,
    'grp_results':      grp_results,
    'grp_order':        GROUP_ORDER,
    'best_thresh':      best_thresh,
    'positive_rate':    positive_rate,
}
out_file = OUTPUT_DIR / 'results.pkl'
with open(out_file, 'wb') as f:
    pickle.dump(results, f)
print(f'Saved: {out_file}')

