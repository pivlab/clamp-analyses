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
# Predicts drug-disease associations using the **gene-based** approach: raw gene-level z-scores from S-PrediXcan (disease, 49 tissues) and LINCS L1000 (drug).
#
# Based on `phenoplier/nbs/30_drug_disease_associations/100-lincs/011-prediction-single_gene_based.ipynb`
#
# For each of 49 tissues and 5 gene-count thresholds (all, 50, 100, 250, 500), runs:
# $$\text{score} = -1 \times \mathbf{drug}^T \mathbf{disease}$$
# on the intersection of genes in LINCS and S-PrediXcan.

# %% [markdown]
# # Modules loading

# %%
# %load_ext autoreload
# %autoreload 2

# %%
from pathlib import Path
from IPython.display import display

import numpy as np
import pandas as pd

from pyprojroot import here

# %% [markdown]
# # Settings

# %%
PREDICTION_METHOD = 'gene_based'

# %%
DATA_DIR = here('data/drug_disease_associations')
display(DATA_DIR)
assert DATA_DIR.exists()

NB_NAME = '06_prediction_single_gene_based'
OUTPUT_DIR = here('output/03_model_biology/00_archs4/02_drug_disease_associations/' + NB_NAME)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Inputs from upstream notebooks
LINCS_RAW_FILE = DATA_DIR / 'lincs-data.pkl'
display(LINCS_RAW_FILE)
assert LINCS_RAW_FILE.exists()

SPREDIXCAN_RAW_DIR = here('output/03_model_biology/00_archs4/02_drug_disease_associations/00_spredixcan_projection_archs4') / 'spredixcan' / 'raw'
display(SPREDIXCAN_RAW_DIR)
assert SPREDIXCAN_RAW_DIR.exists()

OUTPUT_PREDICTIONS_DIR = OUTPUT_DIR / 'lincs' / 'predictions'
OUTPUT_PREDICTIONS_DIR.mkdir(parents=True, exist_ok=True)
display(OUTPUT_PREDICTIONS_DIR)


# %% [markdown]
# # Helper functions
#

# %%
import sys
sys.path.insert(0, str(here('libs')))
from drug_disease_utils import map_traits_to_doid, _zero_nontop_genes, predict_dotprod_neg

# %% [markdown]
# # Load PharmacotherapyDB gold standard

# %%
gold_standard = pd.read_pickle(DATA_DIR / 'gold_standard.pkl')
display(gold_standard.shape)
display(gold_standard.head())

doids_in_gold_standard = set(gold_standard['trait'])
print(f'Unique DOIDs in gold standard: {len(doids_in_gold_standard)}')

# %% [markdown]
# # Load trait → DOID mapping files

# %%
ukb_efo = pd.read_csv(
    DATA_DIR / 'phenomexcan_traits_fullcode_to_efo.tsv',
    sep='\t',
    index_col='ukb_fullcode',
)
# PhenoPlier stores trait full codes with hyphens (e.g. "I70-Diagnoses_...") but
# our S-PrediXcan data uses underscores throughout (e.g. "I70_Diagnoses_...").
# Normalize the index
ukb_efo.index = [idx.replace('-', '_', 1) for idx in ukb_efo.index]

efo_xrefs = pd.read_csv(DATA_DIR / 'term_id_xrefs.tsv.gz', sep='\t')
do_xrefs = pd.read_csv(DATA_DIR / 'xrefs-prop-slim.tsv', sep='\t')

# %% [markdown]
# # Load LINCS raw data

# %%
input_file = LINCS_RAW_FILE
display(input_file)
lincs_data = pd.read_pickle(input_file)
display(lincs_data.shape)

# %% [markdown]
# # Load S-PrediXcan per-tissue files

# %%
spredixcan_file_list = sorted(
    f for f in SPREDIXCAN_RAW_DIR.glob('*.pkl') if f.name.startswith('spredixcan-')
)
display(len(spredixcan_file_list))
assert len(spredixcan_file_list) == 49

# %% [markdown]
# # Predict drug-disease associations

# %%
N_TOP_GENES_LIST = [None, 50, 100, 250, 500]

for spredixcan_file in spredixcan_file_list:
    print(spredixcan_file.name)

    # Load tissue-specific S-PrediXcan data
    tissue_data = pd.read_pickle(spredixcan_file)

    # Intersect genes with LINCS
    common_genes = tissue_data.index.intersection(lincs_data.index)
    tissue_data_common = tissue_data.loc[common_genes]
    lincs_common = lincs_data.loc[common_genes]

    print(f'  shape: {tissue_data_common.shape} ({len(common_genes)} common genes)')

    for ntc in N_TOP_GENES_LIST:
        predict_dotprod_neg(
            lincs_common,
            spredixcan_file,
            tissue_data_common,
            OUTPUT_PREDICTIONS_DIR,
            PREDICTION_METHOD,
            doids_in_gold_standard,
            ukb_efo,
            efo_xrefs,
            do_xrefs,
            n_top_conditions=ntc,
            use_abs=True,
        )

    print('')
