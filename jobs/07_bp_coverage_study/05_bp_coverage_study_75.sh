#!/bin/bash
#SBATCH --job-name=archs4_study_cov_75
#SBATCH --output=log_ARCHS4_study_coverage_75.%j.log
#SBATCH --error=log_ARCHS4_study_coverage_75.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=100GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marc.subiranagranes@cuanschutz.edu

MINIFORGE="/pividori_lab/software/miniforge3"
source "${MINIFORGE}/etc/profile.d/conda.sh"
source /pividori_lab/marc_projects/conda_env/clamp-analyses/bin/activate

n_jobs=16
export NUMBA_NUM_THREADS=$n_jobs
export MKL_NUM_THREADS=$n_jobs
export OPENBLAS_NUM_THREADS=$n_jobs
export NUMEXPR_NUM_THREADS=$n_jobs
export OMP_NUM_THREADS=$n_jobs

start=$(date +%s)
jupyter nbconvert --to notebook --execute \
  nbs/01_model_building/04_archs4/07_bp_coverage_study/05_bp_coverage_study_75.ipynb \
  --output 05_bp_coverage_study_75_executed.ipynb
end=$(date +%s)
echo "Elapsed seconds: $((end - start))"
echo "End: $(date)"
