#!/usr/bin/env bash
set -euo pipefail
cd /home/marc/Documents/pivlab/clamp-analyses
export PATH="$PWD/.conda/envs/clamp-analyses/bin:$PATH"
LOG=logs/run_03_rf_true_labels_bg.log
echo "RUN_START $(date)" > "$LOG"
for nb in \
  nbs/03_model_biology/01_gtex/03_rf_true_labels/00_LV_importance_true_labels.ipynb \
  nbs/03_model_biology/01_gtex/03_rf_true_labels/01_LV_importance_true_labels_biology.ipynb \
  nbs/03_model_biology/01_gtex/03_rf_true_labels/02_LV_importance_true_labels_tissue.ipynb \
  nbs/03_model_biology/01_gtex/03_rf_true_labels/03_global_alignment.ipynb \
  nbs/03_model_biology/01_gtex/03_rf_true_labels/04_tissue_cellmarker_validation.ipynb \
  nbs/03_model_biology/01_gtex/03_rf_true_labels/05_tissue_validation_multi_resource.ipynb
 do
  echo "===== START $nb $(date) =====" | tee -a "$LOG"
  jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 "$nb" >> "$LOG" 2>&1
  status=$?
  echo "===== END $nb status=$status $(date) =====" | tee -a "$LOG"
  if [ "$status" -ne 0 ]; then exit "$status"; fi
 done
echo "RUN_DONE $(date)" | tee -a "$LOG"
