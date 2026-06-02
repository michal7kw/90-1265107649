#!/bin/bash

#===============================================================================
# SCRIPT: 5_DEG_analysis.sh
# PURPOSE: Run Differential Gene Expression Analysis
#
# DESCRIPTION:
# Wrapper script to run DEG analysis using DESeq2 with fold-change filtering.
#
# USAGE:
# sbatch 5_DEG_analysis.sh
#===============================================================================

#SBATCH --job-name=5_DEG_ribotag
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/5_DEG.err"
#SBATCH --output="./logs/5_DEG.out"

# Set up conda environment with R and required packages
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/diffbind_analysis

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"

echo "=== Starting DEG Analysis ==="
echo "Timestamp: $(date)"

cd ${BASE_DIR}/scripts

# Run R script
Rscript 5_DEG_analysis.R

echo "=== DEG Analysis Complete ==="
echo "Timestamp: $(date)"
echo "Results in: ${BASE_DIR}/results/04_DEG/"
