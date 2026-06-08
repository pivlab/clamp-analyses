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

# %% [markdown]
# # Description
#
# Projects S-PrediXcan MASHR gene-trait z-scores (49 tissues) into the ARCHS4 CLAMP latent space.
#
# **Input**: `data/phenomexcan/gene_assoc/spredixcan/pkl/` (49 per-tissue pkl files, genes × traits).
#
# **Outputs**:
# - `00_spredixcan_projection_archs4/spredixcan/raw/{stem}-data.pkl` (genes × traits)
# - `00_spredixcan_projection_archs4/spredixcan/proj/{stem}-projection.pkl` (LVs × traits)
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

SPREDIXCAN_FOLDER = here('data/phenomexcan/gene_assoc/spredixcan/pkl')
display(SPREDIXCAN_FOLDER)
assert SPREDIXCAN_FOLDER.exists()

CLAMP_MODEL_FILE = here('output/01_model_building/04_archs4/06_bp_coverage_rshall/06_bp_coverage_hall_rs_100/hall_coverage_rs100_seed_1/CLAMPfull_hall.rds')
display(CLAMP_MODEL_FILE)
assert CLAMP_MODEL_FILE.exists()

# %%
DRUG_DISEASE_DIR = here('output/03_model_biology/00_archs4/02_drug_disease_associations/00_spredixcan_projection_archs4')
DRUG_DISEASE_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_RAW_DIR = DRUG_DISEASE_DIR / 'spredixcan' / 'raw'
OUTPUT_RAW_DIR.mkdir(parents=True, exist_ok=True)
display(OUTPUT_RAW_DIR)

OUTPUT_PROJ_DIR = DRUG_DISEASE_DIR / 'spredixcan' / 'proj'
OUTPUT_PROJ_DIR.mkdir(parents=True, exist_ok=True)
display(OUTPUT_PROJ_DIR)


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
# S-PrediXcan files use Ensembl IDs as index.
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
# # Load and project S-PrediXcan tissue files

# %%
input_file_list = sorted(SPREDIXCAN_FOLDER.glob('spredixcan-mashr-zscores-*.pkl'))
_tmp = len(input_file_list)
display(_tmp)
assert _tmp == 49, f'Expected 49 tissue files, found {_tmp}'

# %%
for input_file in input_file_list:
    print(input_file.name)

    data = pd.read_pickle(input_file)
    print(f'  shape: {data.shape}')

    n_dups = data.index.duplicated().sum()
    if n_dups > 0:
        print(f'  dropping {n_dups} duplicate gene rows')
        data = data[~data.index.duplicated(keep='first')]

    assert data.index.is_unique
    assert data.columns.is_unique

    data = data.dropna(how='any')
    print(f'  shape (no NaN, no dups): {data.shape}')
    assert not data.isna().any().any()

    output_raw = OUTPUT_RAW_DIR / f'{input_file.stem}-data.pkl'
    print(f'  saving raw to: {output_raw}')
    data.to_pickle(output_raw)

    print('  projecting through CLAMP...')

    aligned = data.reindex(mapped_ensembl).fillna(0.0).values

    r_mat = ro.r['matrix'](
        ro.FloatVector(aligned.flatten('F')),
        nrow=aligned.shape[0],
        ncol=aligned.shape[1],
    )

    proj_r = CLAMP.projectCLAMP(clamp_sub, newdata=r_mat)

    with localconverter(ro.default_converter + pandas2ri.converter):
        proj_values = ro.conversion.rpy2py(proj_r)

    projection = pd.DataFrame(proj_values, index=lv_names, columns=data.columns)
    print(f'    projection shape: {projection.shape}')

    output_proj = OUTPUT_PROJ_DIR / f'{input_file.stem}-projection-{MODEL_KEY}.pkl'
    print(f'    saving projection to: {output_proj}')
    projection.to_pickle(output_proj)
