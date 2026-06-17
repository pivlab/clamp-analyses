# CLAMP analyses

This repository contains the analysis code used for CLAMP benchmarking, projections, biological validation, and downstream analyses.

## 🔧 Dependencies

This repository uses separate Conda environments to avoid dependency conflicts between CPU-based analyses, R/Bioconductor packages, and GPU-accelerated workflows.

| Environment | File | Purpose |
|------------|------|---------|
| `clamp-analyses` | `envs/clamp-analyses.lock` | Frozen environment for reproducing the main CPU-based analyses. |
| `gpu-kmeans` | `envs/gpu-kmeans.lock` | Frozen environment for GPU-accelerated clustering and benchmarking. |

The environments are provided as frozen Conda spec files. This allows collaborators to recreate the same package versions used for the analyses.

## 🛠️ Install dependencies

### 1. Install Conda

Install Miniconda, Mambaforge, or Miniforge for your platform.

### 2. Create the main analysis environment

Create the main frozen analysis environment:

```bash
conda create -n clamp-analyses --file envs/clamp-analyses.lock
conda activate clamp-analyses
```

### 3. Install the pinned CLAMP version

This analysis repository used this pinned CLAMP version while waiting for the Bioconductor submission.

Install CLAMP with the provided script:

```bash
Rscript scripts/install_clamp.R
```

The script installs CLAMP from this pinned commit:

```text
4a6a32006624b942c847becd71f73baf7dedfed6
```

Check that CLAMP loads correctly:

```bash
Rscript -e "library(CLAMP); packageVersion('CLAMP'); sessionInfo()"
```

## Optional: GPU environment

Some clustering and benchmarking workflows can use GPU acceleration through RAPIDS/cuML.

Before creating the GPU environment, verify that a compatible NVIDIA driver and CUDA version are available:

```bash
nvidia-smi
```

Create the frozen GPU environment:

```bash
conda create -n gpu-kmeans --file envs/gpu-kmeans.lock
conda activate gpu-kmeans
```

Install the same pinned CLAMP version if needed:

```bash
Rscript scripts/install_clamp.R
```

## Updating frozen environments

The frozen environment files should be regenerated whenever packages are added or updated.

### Update the main analysis environment lock file

```bash
conda activate clamp-analyses
conda list --explicit > envs/clamp-analyses.lock
conda deactivate
```

### Update the GPU environment lock file

```bash
conda activate gpu-kmeans
conda list --explicit > envs/gpu-kmeans.lock
conda deactivate
```

## 📘 Notebook headers

Each notebook should state which Conda environment is required in the first Markdown cell.

Example:

```markdown
Environment: `clamp-analyses`
```

or:

```markdown
Environment: `gpu-kmeans`
```

## Citation

Citation information will be added once available.

## License

This project is licensed under the [CC-BY 4.0 License](http://creativecommons.org/licenses/by/4.0/).

## Acknowledgments

Supported by the National Human Genome Research Institute,  
the Eunice Kennedy Shriver National Institute of Child Health and Human Development,  
the National Science Foundation, and the National Eye Institute.