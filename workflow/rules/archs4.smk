import os

# Shared ARCHS4 workflow context. Final models are published under
# output/98_final_models and consumed by their dedicated analyses.
A4_CFG = config["archs4"]
A4_PROD = A4_CFG["paths"]["production"]
A4_BIO = A4_CFG["paths"]["biology"]
A4_MODEL_NB = os.path.join(REPO_ROOT, A4_CFG["paths"]["model_notebooks"])
A4_BIO_NB = os.path.join(REPO_ROOT, A4_CFG["paths"]["biology_notebooks"])
A4_PREP = f"{A4_PROD}/00_preprocess"
