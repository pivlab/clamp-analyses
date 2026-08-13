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
# Projects LINCS L1000 consensus drug signatures into the ARCHS4 CLAMP latent space.
#
# **Input**: `data/drug_disease_associations/lincs-data.pkl` (7120 Ensembl genes × 1170 drugs).
#
# **Outputs**:
# - `01_lincs_projection_archs4/lincs/lincs-data.pkl` (genes × drugs, copy)
# - `01_lincs_projection_archs4/lincs/lincs-projection.pkl` (LVs × drugs)
#

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

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

from pyprojroot import here

# %% [markdown]
# # Settings

# %%
MODEL_KEY = 'archs4'

# Input: processed LINCS data
DATA_DIR = here('data/drug_disease_associations')
LINCS_INPUT_FILE = DATA_DIR / 'lincs-data.pkl'
display(LINCS_INPUT_FILE)
assert LINCS_INPUT_FILE.exists()

CLAMP_MODEL_FILE = here('output/01_model_building/04_archs4/06_bp_coverage_rshall/06_bp_coverage_hall_rs_100/hall_coverage_rs100_seed_1/CLAMPfull_hall.rds')
display(CLAMP_MODEL_FILE)
assert CLAMP_MODEL_FILE.exists()


# %%
OUTPUT_DIR = here('output/03_model_biology/00_archs4/02_drug_disease_associations/01_lincs_projection_archs4') / 'lincs'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
display(OUTPUT_DIR)


# %% [markdown]
# # Load LINCS data

# %%
lincs_data = pd.read_pickle(LINCS_INPUT_FILE)
display(lincs_data.shape)
display(lincs_data.head())
assert lincs_data.index.is_unique
assert lincs_data.columns.is_unique
assert not lincs_data.isna().any().any()

# %%
# Verify all index entries are Ensembl IDs (15 chars)
_tmp = pd.Series(lincs_data.index.map(len)).value_counts()
display(_tmp)
assert _tmp.shape[0] == 1

# %% [markdown]
# # Save raw LINCS data to output directory

# %%
output_raw_file = OUTPUT_DIR / 'lincs-data.pkl'
display(output_raw_file)
lincs_data.to_pickle(output_raw_file)
print('Saved.')

# %% [markdown]
# # Load CLAMP model and prepare gene mapping

# %%
CLAMP = importr('CLAMP')
readRDS = ro.r['readRDS']
clamp = readRDS(str(CLAMP_MODEL_FILE))
print('CLAMP model loaded')

# %%
gene_symbols = list(ro.r['rownames'](clamp.rx2('Z')))
lv_names = list(ro.r['colnames'](clamp.rx2('Z')))
print(f'CLAMP genes: {len(gene_symbols)}, LVs: {len(lv_names)}')

# %%
# Map CLAMP gene symbols (HGNC) → Ensembl IDs
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

# %%
# Keep only 1:1 unambiguous symbol ↔ Ensembl mappings
dup_symbols = mapping_df['SYMBOL'].duplicated(keep=False)
dup_ensembl = mapping_df['ENSEMBL'].duplicated(keep=False)
mapping_1to1 = mapping_df[~dup_symbols & ~dup_ensembl].set_index('SYMBOL')
print(f'1:1 mappings: {mapping_1to1.shape[0]} / {len(gene_symbols)} CLAMP genes')

mapped_symbols = mapping_1to1.index.tolist()
mapped_ensembl = mapping_1to1['ENSEMBL'].tolist()

# %%
# Build CLAMP sub-object with Z restricted to 1:1-mapped genes.
subset_Z = ro.r('function(clamp, genes) { clamp$Z <- as.matrix(clamp$Z[genes, ]); clamp }')
clamp_sub = subset_Z(clamp, ro.StrVector(mapped_symbols))
print(f'Subsetted CLAMP Z: {len(mapped_symbols)} genes x {len(lv_names)} LVs')

# %% [markdown]
# # Project LINCS into CLAMP

# %%
# Align LINCS to CLAMP gene order (mapped Ensembl IDs), fill missing genes with 0
aligned = lincs_data.reindex(mapped_ensembl).fillna(0.0).values  # (n_genes, n_drugs)

r_mat = ro.r['matrix'](
    ro.FloatVector(aligned.flatten('F')),
    nrow=aligned.shape[0],
    ncol=aligned.shape[1],
)

proj_r = CLAMP.projectCLAMP(clamp_sub, newdata=r_mat)

with localconverter(ro.default_converter + pandas2ri.converter):
    proj_values = ro.conversion.rpy2py(proj_r)

lincs_projection = pd.DataFrame(proj_values, index=lv_names, columns=lincs_data.columns)
print(f'LINCS projection shape: {lincs_projection.shape}')
display(lincs_projection.head())

# %%
assert not lincs_projection.isna().any().any()

# %% [markdown]
# # Save

# %%
output_proj_file = OUTPUT_DIR / 'lincs-projection.pkl'
display(output_proj_file)
lincs_projection.to_pickle(output_proj_file)
print('Saved.')
