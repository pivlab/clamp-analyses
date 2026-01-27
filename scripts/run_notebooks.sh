#!/bin/bash

#===========================================
# CONFIGURATION
# Format: "conda_env | notebook_path"
#===========================================

JOBS=(
    "clamp-analyses | nbs/01_model_building/02_gtex/09_cogaps.ipynb"
    "clamp-analyses | nbs/01_model_building/02_gtex/00_gtex.ipynb"
    "clamp-analyses | nbs/01_model_building/02_gtex/08_PLIER.ipynb"
    "clamp-analyses | nbs/01_model_building/02_gtex/01_gtex_priors.ipynb"
    "clamp-analyses | nbs/01_model_building/02_gtex/07_wgcna.ipynb"
)

#===========================================
# EXECUTION
#===========================================

echo "Starting notebook execution sequence..."
echo "========================================="

for job in "${JOBS[@]}"; do
    env=$(echo "$job" | cut -d'|' -f1 | xargs)
    nb=$(echo "$job" | cut -d'|' -f2 | xargs)
    
    echo ""
    echo "Running: $nb"
    echo "Environment: $env"
    echo "Started at: $(date)"
    echo "-----------------------------------------"
    
    if ! conda run -n "$env" --no-capture-output jupyter execute "$nb"; then
        echo ""
        echo "FAILED: $nb"
        echo "Stopping execution."
        exit 1
    fi
    
    echo "Finished: $nb"
    echo "Completed at: $(date)"
    echo "========================================="
done

echo ""
echo "All notebooks completed successfully!"