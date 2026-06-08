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
# Predicts drug-disease associations using the recount2 **module-based** CLAMP model.
#
# Reads LINCS and S-PrediXcan projections from:
# - `04_spredixcan_projection_recount2` (S-PrediXcan, LVs × traits)
# - `05_lincs_projection_recount2` (LINCS, LVs × drugs)
#
# For each of 49 tissues and 5 LV-count thresholds (all, 5, 10, 25, 50), scores:
# $$\text{score} = -1 \times \mathbf{drug}^T \mathbf{disease}$$
#

# %% [markdown]
# # Modules loading
#

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
#

# %%
PREDICTION_METHOD = 'module_based_recount2'
MODEL_KEY = PREDICTION_METHOD.removeprefix('module_based_')

# %%
DATA_DIR = here('data/drug_disease_associations')
display(DATA_DIR)
assert DATA_DIR.exists()

NB_NAME = '09_prediction_module_based_recount2'
OUTPUT_DIR = here('output/03_model_biology/00_archs4/02_drug_disease_associations/' + NB_NAME)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Inputs from upstream projection notebooks
LINCS_PROJ_FILE = here('output/03_model_biology/00_archs4/02_drug_disease_associations/05_lincs_projection_recount2') / 'lincs' / 'lincs-projection.pkl'
display(LINCS_PROJ_FILE)
assert LINCS_PROJ_FILE.exists()

SPREDIXCAN_PROJ_DIR = here('output/03_model_biology/00_archs4/02_drug_disease_associations/04_spredixcan_projection_recount2') / 'spredixcan'
display(SPREDIXCAN_PROJ_DIR)
assert SPREDIXCAN_PROJ_DIR.exists()

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
#

# %%
gold_standard = pd.read_pickle(DATA_DIR / 'gold_standard.pkl')
display(gold_standard.shape)
display(gold_standard.head())

doids_in_gold_standard = set(gold_standard['trait'])
print(f'Unique DOIDs in gold standard: {len(doids_in_gold_standard)}')

# %% [markdown]
# # Load trait → DOID mapping files
#

# %%
ukb_efo = pd.read_csv(
    DATA_DIR / 'phenomexcan_traits_fullcode_to_efo.tsv',
    sep='\t',
    index_col='ukb_fullcode',
)
# PhenoPlier stores trait full codes with hyphens (e.g. "I70-Diagnoses_...") but
# our S-PrediXcan data uses underscores throughout (e.g. "I70_Diagnoses_...").
# Normalize the index so lookups work correctly.
ukb_efo.index = [idx.replace('-', '_', 1) for idx in ukb_efo.index]

efo_xrefs = pd.read_csv(DATA_DIR / 'term_id_xrefs.tsv.gz', sep='\t')
do_xrefs = pd.read_csv(DATA_DIR / 'xrefs-prop-slim.tsv', sep='\t')

# %% [markdown]
# # Load LINCS projection
#

# %%
lincs_projection = pd.read_pickle(LINCS_PROJ_FILE)
print(f'LINCS projection shape: {lincs_projection.shape}')
assert not lincs_projection.isna().any().any()
display(lincs_projection.head())


# %% [markdown]
# # Load S-PrediXcan projected files
#

# %%
spredixcan_file_list = sorted(
    f for f in SPREDIXCAN_PROJ_DIR.glob(f'spredixcan-*-projection-{MODEL_KEY}.pkl')
)
display(len(spredixcan_file_list))
assert len(spredixcan_file_list) == 49

display(pd.read_pickle(spredixcan_file_list[0]).head())


# %% [markdown]
# # Predict drug-disease associations
#

# %%
# LV thresholds: None = all LVs (same as PhenoPlier module-based thresholds)
N_TOP_LVS_LIST = [None, 5, 10, 25, 50]

for spredixcan_file in spredixcan_file_list:
    print(spredixcan_file.name)

    tissue_proj = pd.read_pickle(spredixcan_file)
    print(f'  shape: {tissue_proj.shape}')
    assert tissue_proj.index.equals(lincs_projection.index)

    for ntc in N_TOP_LVS_LIST:
        predict_dotprod_neg(
            lincs_projection,
            spredixcan_file,
            tissue_proj,
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

