# CLAMP analyses


## 🔧 Dependencies

CLAMP uses **three Conda environments** to separate core package development, large-scale analyses, and GPU-accelerated workflows:

| Environment | File | Purpose |
|--------------|------|----------|
| **`clamp-analyses.yaml`** | `envs/clamp-analyses.yaml` | Default environment for CPU-based modeling, priors, projections, and vignettes. |
| **`gpu-kmeans.yaml`** | `envs/gpu-kmeans.yaml` | Optional environment for GPU-accelerated clustering and benchmarking. |

This separation avoids dependency conflicts between R (Bioconductor) and GPU libraries (`rapids`, `cupy`, `cuml`).

### 🛠️ Install dependencies

Recommended steps to install the system-level and Conda tooling required to create the CLAMP environments.

1. Install a Conda distribution
- Install Miniconda or Mambaforge for your platform (Mambaforge is recommended for faster environment solves).

2. (Optional) Verify GPU drivers for RAPIDS/cuML workflows
- Ensure a compatible NVIDIA driver / CUDA version is installed before creating the GPU environment:
```bash
nvidia-smi
```
- Check RAPIDS compatibility matrix for the correct CUDA version (match driver/CUDA with RAPIDS/cuML requirements).

3. Create environments using conda

```bash
conda env create -f envs/clamp-analyses.yaml
conda activate clamp-analyses

# Clone CLAMP the repo into REPO_PATH (adjust path as needed)
export REPO_PATH=~/path/to/CLAMP
mkdir -p "$(dirname "$REPO_PATH")"
git clone https://github.com/chikinalab/CLAMP.git "$REPO_PATH"

# Install and check CLAMP using devtools
Rscript -e "devtools::install_local('$REPO_PATH', force=TRUE, dependencies=FALSE)"
Rscript -e "library(CLAMP); cat('CLAMP version:', packageVersion('CLAMP'), '\n')"
```

```bash
conda env create -f envs/envs/gpu-kmeans.yaml
conda activate gpu-kmeans.yaml

# Clone CLAMP the repo into REPO_PATH (adjust path as needed)
export REPO_PATH=~/path/to/CLAMP
mkdir -p "$(dirname "$REPO_PATH")"
git clone https://github.com/chikinalab/CLAMP.git "$REPO_PATH"

# Install and check CLAMP using devtools
Rscript -e "devtools::install_local('$REPO_PATH', force=TRUE, dependencies=FALSE)"
```

## 📘 Notebook Headers

Each notebook explicitly states which environment to use in the first Markdown cell.

## Citation

## License

This project is licensed under the [CC-BY 4.0 License](http://creativecommons.org/licenses/by/4.0/).

## Acknowledgments

Supported by the National Human Genome Research Institute,  
The Eunice Kennedy Shriver National Institute of Child Health and Human Development,  
the National Science Foundation, and the National Eye Institute.
