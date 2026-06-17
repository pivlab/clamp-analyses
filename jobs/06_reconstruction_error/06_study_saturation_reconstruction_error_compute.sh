#!/bin/bash
#SBATCH --job-name=archs4_reconstruction_error_saturation_study
#SBATCH --output=log_ARCHS4_reconstruction_error_saturation_study.%j.log
#SBATCH --error=log_ARCHS4_reconstruction_error_saturation_study.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=64GB
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
  nbs/01_model_building/06_reconstruction_error/06_study_saturation_reconstruction_error_compute.ipynb \
  --output 06_study_saturation_reconstruction_error_compute_executed.ipynb
end=$(date +%s)
echo "Elapsed seconds: $((end - start))"
echo "End: $(date)"
