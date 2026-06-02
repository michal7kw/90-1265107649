#!/bin/bash
#SBATCH --job-name=6b_viz_enhanced
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/6b_viz_enhanced.err"
#SBATCH --output="./logs/6b_viz_enhanced.out"

#===============================================================================
# Wrapper script for 6b_visualize_splicing_enhanced.R
# Creates enhanced visualizations for differential splicing analysis
#===============================================================================

echo "============================================"
echo "Enhanced Splicing Visualization"
echo "Start time: $(date)"
echo "============================================"

# Activate conda environment with R
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/rnaseq-quant

# Change to scripts directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649/scripts

# Run R script
Rscript 6b_visualize_splicing_enhanced.R

echo ""
echo "============================================"
echo "End time: $(date)"
echo "============================================"
