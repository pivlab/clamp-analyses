#!/bin/bash
#SBATCH --job-name=archs4_sat_study_01
#SBATCH --output=log_ARCHS4_saturation_study_01.%j.log
#SBATCH --error=log_ARCHS4_saturation_study_01.%j.err
#SBATCH --time=96:00:00
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
  nbs/01_model_building/04_archs4/08_saturation_study/00_saturation_study_rs_01.ipynb \
  --output 00_saturation_study_rs_01_executed.ipynb
end=$(date +%s)
echo "Elapsed seconds: $((end - start))"
echo "End: $(date)"
