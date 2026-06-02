#!/bin/bash

#===============================================================================
# SCRIPT: 4_multiqc.sh
# PURPOSE: Aggregate QC Reports with MultiQC
#
# DESCRIPTION:
# Generates a comprehensive MultiQC report combining FastQC and STAR alignment
# statistics for all samples.
#
# USAGE:
# sbatch 4_multiqc.sh
#===============================================================================

#SBATCH --job-name=4_multiqc_ribotag
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/4_multiqc.err"
#SBATCH --output="./logs/4_multiqc.out"

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/multiqc

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"

echo "=== Running MultiQC ==="
echo "Timestamp: $(date)"

cd ${BASE_DIR}

# Run MultiQC on all results
multiqc \
    --force \
    --outdir ${BASE_DIR}/results \
    --filename multiqc_report \
    --title "RiboTag Analysis - Project 90-1265107649" \
    ${BASE_DIR}/results/01_fastqc \
    ${BASE_DIR}/results/02_aligned

echo "=== MultiQC complete ==="
echo "Report: ${BASE_DIR}/results/multiqc_report.html"
echo "Timestamp: $(date)"
