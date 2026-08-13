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
# Projects LINCS L1000 consensus drug signatures into the GTEx CLAMP latent space.
#
# **Input**: `data/drug_disease_associations/lincs-data.pkl` (7120 Ensembl genes × 1170 drugs).
#
# **Output**: `03_lincs_projection_gtex/lincs/lincs-projection.pkl` (LVs × drugs).
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

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

from pyprojroot import here


# %% [markdown]
# # Settings
#

# %%
CLAMP_MODEL_FILE = here('output/01_model_building/02_gtex/10_CLAMP_hall/CLAMPfull_hall.rds')
display(CLAMP_MODEL_FILE)
assert CLAMP_MODEL_FILE.exists()


# %%
# Input: processed LINCS data (Ensembl IDs, 7120 PhenomeXcan genes)
DATA_DIR = here('data/drug_disease_associations')
LINCS_INPUT_FILE = DATA_DIR / 'lincs-data.pkl'
display(LINCS_INPUT_FILE)
assert LINCS_INPUT_FILE.exists()

# Output
OUTPUT_DIR = here('output/03_model_biology/00_archs4/02_drug_disease_associations/03_lincs_projection_gtex') / 'lincs'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
display(OUTPUT_DIR)


# %% [markdown]
# # Projection helpers
#

# %%
def prepare_clamp_projector(clamp_model_file):
    CLAMP = importr('CLAMP')
    readRDS = ro.r['readRDS']
    clamp = readRDS(str(clamp_model_file))
    print('CLAMP model loaded')

    gene_symbols = list(ro.r['rownames'](clamp.rx2('Z')))
    lv_names = list(ro.r['colnames'](clamp.rx2('Z')))
    print(f'CLAMP genes: {len(gene_symbols)}, LVs: {len(lv_names)}')

    # Map CLAMP gene symbols (HGNC) to Ensembl IDs used by LINCS and S-PrediXcan.
    clusterProfiler = importr('clusterProfiler')
    bitr_result = clusterProfiler.bitr(
        ro.StrVector(gene_symbols),
        fromType='SYMBOL',
        toType='ENSEMBL',
        OrgDb='org.Hs.eg.db',
    )

    with localconverter(ro.default_converter + pandas2ri.converter):
        mapping_df = ro.conversion.rpy2py(bitr_result)

    print(f'Raw mapping shape: {mapping_df.shape}')
    display(mapping_df.head())

    # Keep only 1:1 unambiguous symbol <-> Ensembl mappings.
    dup_symbols = mapping_df['SYMBOL'].duplicated(keep=False)
    dup_ensembl = mapping_df['ENSEMBL'].duplicated(keep=False)
    mapping_1to1 = mapping_df[~dup_symbols & ~dup_ensembl].set_index('SYMBOL')
    print(f'1:1 mappings: {mapping_1to1.shape[0]} / {len(gene_symbols)} CLAMP genes')

    mapped_symbols = mapping_1to1.index.tolist()
    mapped_ensembl = mapping_1to1['ENSEMBL'].tolist()

    subset_Z = ro.r('function(clamp, genes) { clamp$Z <- as.matrix(clamp$Z[genes, ]); clamp }')
    clamp_sub = subset_Z(clamp, ro.StrVector(mapped_symbols))
    print(f'Subsetted CLAMP Z: {len(mapped_symbols)} genes x {len(lv_names)} LVs')

    return CLAMP, clamp_sub, mapped_ensembl, lv_names


def project_to_clamp(data, CLAMP, clamp_sub, mapped_ensembl, lv_names):
    aligned = data.reindex(mapped_ensembl).fillna(0.0).values

    r_mat = ro.r['matrix'](
        ro.FloatVector(aligned.flatten('F')),
        nrow=aligned.shape[0],
        ncol=aligned.shape[1],
    )

    proj_r = CLAMP.projectCLAMP(clamp_sub, newdata=r_mat)

    with localconverter(ro.default_converter + pandas2ri.converter):
        proj_values = ro.conversion.rpy2py(proj_r)

    return pd.DataFrame(proj_values, index=lv_names, columns=data.columns)



# %% [markdown]
# # Load LINCS data
#

# %%
input_file = LINCS_INPUT_FILE
display(input_file)
lincs_data = pd.read_pickle(input_file)
display(lincs_data.shape)
display(lincs_data.head())
assert lincs_data.index.is_unique
assert lincs_data.columns.is_unique
assert not lincs_data.isna().any().any()


# %% [markdown]
# # Prepare CLAMP projector
#

# %%
CLAMP, clamp_sub, mapped_ensembl, lv_names = prepare_clamp_projector(CLAMP_MODEL_FILE)


# %% [markdown]
# # Project LINCS
#

# %%
lincs_projection_file = OUTPUT_DIR / 'lincs-projection.pkl'
display(lincs_projection_file)

lincs_projection = project_to_clamp(
    lincs_data, CLAMP, clamp_sub, mapped_ensembl, lv_names
)
print(f'LINCS projection shape: {lincs_projection.shape}')
assert not lincs_projection.isna().any().any()
lincs_projection.to_pickle(lincs_projection_file)
print('Saved.')

display(lincs_projection.head())

