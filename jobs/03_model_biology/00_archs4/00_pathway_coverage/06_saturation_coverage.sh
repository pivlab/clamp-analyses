#!/bin/bash
#SBATCH --job-name=saturation_coverage
#SBATCH --output=log_saturation_coverage.%j.log
#SBATCH --error=log_saturation_coverage.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marc.subiranagranes@cuanschutz.edu

MINIFORGE="/pividori_lab/software/miniforge3"
source "${MINIFORGE}/etc/profile.d/conda.sh"
source /pividori_lab/marc_projects/conda_env/clamp-analyses/bin/activate

n_jobs=4
export NUMBA_NUM_THREADS=$n_jobs
export MKL_NUM_THREADS=$n_jobs
export OPENBLAS_NUM_THREADS=$n_jobs
export NUMEXPR_NUM_THREADS=$n_jobs
export OMP_NUM_THREADS=$n_jobs

start=$(date +%s)
jupyter nbconvert --to notebook --execute nbs/03_model_biology/00_archs4/00_pathway_coverage/06_saturation_coverage.ipynb --output 06_saturation_coverage_executed.ipynb
end=$(date +%s)
echo "Elapsed seconds: $((end - start))"
echo "End: $(date)"
