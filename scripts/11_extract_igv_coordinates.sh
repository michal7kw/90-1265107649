#!/bin/bash
#SBATCH --job-name=11_igv_coords
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/11_igv_coords.err"
#SBATCH --output="./logs/11_igv_coords.out"

#===============================================================================
# Wrapper script for 11_extract_igv_coordinates.R
# Extracts significant splicing events with IGV-ready coordinates
#===============================================================================

echo "============================================"
echo "IGV Coordinate Extraction"
echo "Start time: $(date)"
echo "============================================"

# Activate conda environment with R
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/rnaseq-quant

# Change to scripts directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649/scripts

# Run R script
Rscript 11_extract_igv_coordinates.R

echo ""
echo "============================================"
echo "End time: $(date)"
echo "============================================"
