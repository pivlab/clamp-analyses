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
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# 💡 **Environment:** `clamp-analyses`  
#

# %% [markdown]
# # Description
#
# Reads the drug-disease prediction HDF5 files from the notebooks 06–09 (gene-based and module-based) and computes final performance measures (AUROC, average precision) against the PharmacotherapyDB gold standard.
#
# Aggregation:
# 1. Group by (trait, drug, method, tissue) and average ranks across thresholds.
# 2. Group by (trait, drug, method) and take the max across tissues.

# %% [markdown]
# # Modules loading

# %%
# %load_ext autoreload
# %autoreload 2

# %%
from pathlib import Path
from collections import defaultdict
from IPython.display import display

import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.metrics import roc_auc_score, average_precision_score

from pyprojroot import here

# %% [markdown]
# # Settings

# %%
N_TISSUES = 49
N_THRESHOLDS = 5

METHOD_RENAME = {
    'Gene-based':             'gene_based',
    'Module-based (ARCHS4)':  'module_based_archs4',
    'Module-based':           'module_based_archs4',
    'Module-based (GTEx)':    'module_based_gtex',
    'Module-based (recount2)':'module_based_recount2',
}

METHOD_THRESHOLDS = {
    'gene_based': [-1.0, 50.0, 100.0, 250.0, 500.0],
    'module_based_archs4': [-1.0, 5.0, 10.0, 25.0, 50.0],
    'module_based_gtex': [-1.0, 5.0, 10.0, 25.0, 50.0],
    'module_based_recount2': [-1.0, 5.0, 10.0, 25.0, 50.0],
}
EXPECTED_METHODS = tuple(METHOD_THRESHOLDS)
N_PREDICTION_FILES_TOTAL = N_TISSUES * sum(len(v) for v in METHOD_THRESHOLDS.values())

# %%
DATA_DIR = here('data/drug_disease_associations')
display(DATA_DIR)
assert DATA_DIR.exists()

# Collect prediction HDF5 files from all 4 prediction notebooks
PREDICTIONS_DIRS = [
    here('output/03_model_biology/00_archs4/02_drug_disease_associations/06_prediction_single_gene_based') / 'lincs' / 'predictions' / 'dotprod_neg',
    here('output/03_model_biology/00_archs4/02_drug_disease_associations/07_prediction_module_based_archs4') / 'lincs' / 'predictions' / 'dotprod_neg',
    here('output/03_model_biology/00_archs4/02_drug_disease_associations/08_prediction_module_based_gtex') / 'lincs' / 'predictions' / 'dotprod_neg',
    here('output/03_model_biology/00_archs4/02_drug_disease_associations/09_prediction_module_based_recount2') / 'lincs' / 'predictions' / 'dotprod_neg',
]
for d in PREDICTIONS_DIRS:
    display(d)
    assert d.exists(), d

OUTPUT_DIR = here('output/03_model_biology/00_archs4/02_drug_disease_associations/10_prediction_performance') / 'lincs'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
(OUTPUT_DIR / 'predictions').mkdir(parents=True, exist_ok=True)
display(OUTPUT_DIR)

# %% [markdown]
# # Load PharmacotherapyDB gold standard

# %%
gold_standard = pd.read_pickle(DATA_DIR / 'gold_standard.pkl')
display(gold_standard.shape)
display(gold_standard.head())
display(gold_standard['true_class'].value_counts())

# %% [markdown]
# # Load drug-disease predictions

# %%
current_prediction_files = []
for d in PREDICTIONS_DIRS:
    current_prediction_files.extend(sorted(d.glob('*.h5')))
current_prediction_files.sort()
display(len(current_prediction_files))

# %%
current_prediction_files[:5]


# %%
def _get_tissue(data_value):
    """
    Extracts tissue name from the metadata 'data' field.
    Handles raw data and dataset-specific projected suffixes.
    """
    prefix = 'spredixcan-mashr-zscores-'
    assert data_value.startswith(prefix), data_value
    tissue = data_value[len(prefix):]

    for suffix in (
        '-projection-archs4',
        '-projection-gtex',
        '-projection-recount2',
        '-projection',
        '-data',
    ):
        if tissue.endswith(suffix):
            return tissue[:-len(suffix)]

    raise ValueError(f'Cannot extract tissue from metadata data value: {data_value}')


# %%
# Load all prediction files, rank scores, merge with gold standard
predictions = []
skipped_files = []

for f in tqdm(current_prediction_files, ncols=100):
    metadata = pd.read_hdf(f, key='metadata')
    method_name = METHOD_RENAME.get(metadata['method'].values[0], metadata['method'].values[0])
    if method_name not in METHOD_THRESHOLDS:
        skipped_files.append((f.name, method_name))
        continue

    # Load DOID-mapped predictions and rank within the full DOID distribution
    prediction_data = pd.read_hdf(f, key='prediction')
    prediction_data['score'] = prediction_data['score'].rank()

    # Filter to gold-standard pairs
    prediction_data = pd.merge(
        prediction_data, gold_standard, on=['trait', 'drug'], how='inner'
    )
    prediction_data['trait'] = prediction_data['trait'].astype('category')
    prediction_data['drug'] = prediction_data['drug'].astype('category')

    prediction_data = prediction_data.assign(method=method_name)
    prediction_data['method'] = pd.Categorical(
        prediction_data['method'], categories=EXPECTED_METHODS, ordered=True
    )

    prediction_data = prediction_data.assign(n_top_genes=metadata['n_top_genes'].values[0])

    data_value = metadata['data'].values[0]
    prediction_data = prediction_data.assign(data=data_value)
    prediction_data['data'] = prediction_data['data'].astype('category')

    # Extract tissue name from data field
    prediction_data = prediction_data.assign(tissue=_get_tissue(data_value))

    predictions.append(prediction_data)

display(f'Skipped files: {len(skipped_files)}')
if skipped_files:
    display(skipped_files[:10])

# %%
predictions = pd.concat(predictions, ignore_index=True)

# %%
display(predictions.shape)
display(predictions.head())

# %% [markdown]
# ## Validation checks

# %%
assert not predictions.isna().any().any()

_method_counts = predictions['method'].value_counts().reindex(EXPECTED_METHODS)
display(_method_counts)

N_PREDICTIONS = predictions[['drug', 'trait']].drop_duplicates().shape[0]
display(f'Unique drug-disease pairs: {N_PREDICTIONS}')

for method_name, thresholds in METHOD_THRESHOLDS.items():
    expected = N_TISSUES * len(thresholds) * N_PREDICTIONS
    actual = int(_method_counts.loc[method_name])
    assert actual == expected, f'{method_name}: expected {expected}, got {actual}'

# %%
_tmp = predictions.groupby(['method', 'n_top_genes'], observed=True).size().rename('n')
display(_tmp)

for method_name, thresholds in METHOD_THRESHOLDS.items():
    method_counts = _tmp.loc[method_name]
    actual_thresholds = sorted(float(x) for x in method_counts.index)
    expected_thresholds = sorted(thresholds)
    assert actual_thresholds == expected_thresholds, (
        f'{method_name}: expected thresholds {expected_thresholds}, got {actual_thresholds}'
    )
    assert np.all(method_counts.values == N_TISSUES * N_PREDICTIONS), method_name

# %% [markdown]
# ## Save raw predictions

# %%
output_file = OUTPUT_DIR / 'predictions' / 'predictions_results.pkl'
display(output_file)
predictions.to_pickle(output_file)


# %% [markdown]
# # Aggregate predictions
#
# 1. Average ranks across n_top_genes thresholds (per trait, drug, method, tissue).
# 2. Take maximum across tissues (per trait, drug, method).
#
# This matches the PhenoPlier aggregation exactly.

# %%
def _reduce_mean(x):
    return pd.Series({
        'score': x['score'].mean(),
        'true_class': x['true_class'].unique()[0],
    })


def _reduce_max(x):
    return pd.Series({
        'score': x['score'].max(),
        'true_class': x['true_class'].unique()[0],
    })


# %%
predictions_avg = (
    # Step 1: average across n_top_genes thresholds
    predictions
    .groupby(['trait', 'drug', 'method', 'tissue'], observed=True)
    .apply(_reduce_mean, include_groups=False)
    .dropna()
    # Step 2: take maximum across tissues
    .groupby(['trait', 'drug', 'method'], observed=True)
    .apply(_reduce_max, include_groups=False)
    .dropna()
    .sort_index()
    .reset_index()
)

# %%
display(predictions_avg.shape)
display(predictions_avg.head())

assert predictions_avg.shape[0] == len(EXPECTED_METHODS) * N_PREDICTIONS
assert predictions_avg.dropna().shape == predictions_avg.shape

# %% [markdown]
# ## Save aggregated predictions

# %%
output_file = OUTPUT_DIR / 'predictions' / 'predictions_results_aggregated.pkl'
display(output_file)
predictions_avg.to_pickle(output_file)

# %% [markdown]
# # ROC performance

# %%
# AUROC per method and n_top_genes threshold
predictions.groupby(['method', 'tissue', 'n_top_genes'], observed=True).apply(
    lambda x: roc_auc_score(x['true_class'], x['score']), include_groups=False
).groupby(['method', 'n_top_genes'], observed=True).describe()

# %%
# Final AUROC using aggregated predictions
auroc_final = predictions_avg.groupby('method', observed=True).apply(
    lambda x: roc_auc_score(x['true_class'], x['score']), include_groups=False
).rename('AUROC')
display(auroc_final)

# %% [markdown]
# These are the final performance measures using AUROC.

# %% [markdown]
# # Precision-Recall performance

# %%
# Average precision per method and n_top_genes threshold
predictions.groupby(['method', 'tissue', 'n_top_genes'], observed=True).apply(
    lambda x: average_precision_score(x['true_class'], x['score']), include_groups=False
).groupby(['method', 'n_top_genes'], observed=True).describe()

# %%
# Final average precision using aggregated predictions
ap_final = predictions_avg.groupby('method', observed=True).apply(
    lambda x: average_precision_score(x['true_class'], x['score']), include_groups=False
).rename('AvgPrecision')
display(ap_final)

# %% [markdown]
# These are the final performance measures using average precision.

# %% [markdown]
# # Summary

# %%
summary = pd.concat([auroc_final, ap_final], axis=1)
display(summary)
